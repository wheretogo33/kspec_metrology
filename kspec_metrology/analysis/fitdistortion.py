import numpy as np
from scipy.optimize import curve_fit
from kspec_metrology.analysis.utils import (transform_polynomial, transform, camera2focal_coeff,
                                            model_mask, COEFF_NAMES)
from kspec_metrology.logging.log import get_logger


def fitdistortion(x, y, fid_flag
                  , xobs, yobs
                  , imatch, theta_guess, xoff_guess, yoff_guess
                  , model='full', m_init=None):
    """fiducial 로 camera -> focal 변환 계수를 맞추고 모든 점을 focal 좌표로 옮긴다.

    model   : 풀 왜곡 항 (utils.DISTORTION_MODELS 의 이름이나 계수 이름 목록).
              'full'(기본) 은 15 개 항 전부, 실험실처럼 단순한 광학계에서는 'radial'
              이나 'similarity' 로 줄인다. 풀지 않는 항은 0 으로 고정된다.
    m_init  : 배율 초기값. None 이면 utils.CAMERA2FOCAL_M.

    반환: xfocal, yfocal, dx, dy, coeff (길이 15 의 전체 계수)
    """
    log = get_logger()

    xobs_match = xobs[imatch]
    yobs_match = yobs[imatch]

    xfit = np.concatenate((xobs_match[fid_flag], xobs_match[fid_flag]))
    yfit = np.concatenate((yobs_match[fid_flag], yobs_match[fid_flag]))
    ccd_fit = np.concatenate((x[fid_flag], y[fid_flag]))

    mask = model_mask(model)
    nfree, nfid = int(mask.sum()), int(np.count_nonzero(fid_flag))
    if 2*nfid < nfree:
        raise ValueError(f"distortion model {model!r} has {nfree} free coefficients but only "
                         f"{nfid} fiducials ({2*nfid} equations); use a simpler model")

    p_full = np.copy(camera2focal_coeff)
    p_full[4:][~mask[4:]] = 0.                 # 풀지 않는 왜곡 항은 0 으로 고정
    p_full[1] = theta_guess
    p_full[2] = xoff_guess
    p_full[3] = yoff_guess
    if m_init is not None:
        p_full[0] = m_init

    def f(xin, *pfree):
        p = p_full.copy()
        p[mask] = pfree
        return transform_polynomial(xin, *p)

    popt, _ = curve_fit(f, (xfit, yfit), ccd_fit, p0=p_full[mask], maxfev=20000)
    inv_popt_obs = p_full.copy()
    inv_popt_obs[mask] = popt

    xfocal_obs, yfocal_obs = transform(xobs_match, yobs_match, inv_popt_obs)

    dx = xfocal_obs-x
    dy = yfocal_obs-y

    res = np.hypot(dx[fid_flag], dy[fid_flag])*1e3
    log.info("Distortion fit (%s, %d coefficients, %d fiducials): fiducial residual "
             "rms %.1f um, max %.1f um; magnification %.4f",
             model if isinstance(model, str) else 'custom', nfree, nfid,
             np.sqrt(np.mean(res**2)), res.max(), inv_popt_obs[0])

    # 각도 변환은 arm 길이가 필요하므로 mtlcal에서 한다.
    return xfocal_obs, yfocal_obs, dx, dy, inv_popt_obs
