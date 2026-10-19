"""실험실 시험용 도우미: 목표 파일 만들기, fiber 표 줄이기, 단계별 확인 그림.

실험실에서는 목표가 관측 타일이 아니라 사람이 정한 (x, y) 나 (alpha, beta) 로
주어진다. 이 모듈이 그것을 mtlcal 이 읽는 target 파일(object.info 형식) 로
바꾸고, 꽂힌 positioner 만 남긴 fiber 표를 만들고, 각 단계의 결과를 그림으로
보여 준다. 분석 자체는 analysis/ 의 함수를 그대로 쓴다.

    make_fiber_table   꽂힌 홀만 positioner 로 둔 Fiber_Configuration 표
    ab2xy, xy2ab       positioner 각도 <-> 초점면 좌표 (패키지 각도 규약 그대로)
    make_target_file   'ID x y' 또는 'ID alpha beta' 텍스트 -> target json
    inspect_trial      한 trial 의 이미지를 단계별로 다시 분석해 그림으로 (결과는 안 바꿈)
    show_image, show_detections, show_matching, show_offsets   확인 그림
"""

import json
import os

import numpy as np
from astropy.io import ascii, fits

from kspec_metrology.analysis import fiber_config as cfg
from kspec_metrology.analysis.utils import find_angle_double_method2_select_elbow
from kspec_metrology.analysis.findpeak import findpeak
from kspec_metrology.analysis.matchfiber import matchfiber
from kspec_metrology.analysis.fitdistortion import fitdistortion
from kspec_metrology.naming import image_path

PIX_MM = 3.76e-3          # 카메라 화소 크기 [mm]
rtod = 180./np.pi


# ============================================================================ 표와 목표
def make_fiber_table(positioners, out_path, src_table=None, fiducials=None):
    """src_table 에서 positioners 에 적힌 홀만 FiducialFlag 0 으로, fiducials (None 이면
    원래 표의 fiducial 전부) 만 1 로 두고 나머지는 -9 로 바꾼 표를 out_path 에 쓴다
    (config73 표와 같은 고정폭 CDS 형식). 쓴 뒤 KSPEC_FIBER_TABLE 을 이 경로로 두고
    패키지를 import 한다 (표는 import 때 읽는다).
    반환: (positioner 수, fiducial 수)"""
    src_table = src_table or cfg.DEFAULT_FIBER_TABLE_PATH
    pos = [str(p).strip() for p in positioners]
    tab = ascii.read(src_table)
    c = tab.colnames
    ids = [str(v).strip() for v in tab[c[0]]]
    flags = {i: int(f) for i, f in zip(ids, tab[c[3]])}
    unknown = [p for p in pos if p not in flags]
    if unknown:
        raise ValueError(f"holes not in {src_table}: {unknown}")
    fid = set(i for i, f in flags.items() if f == 1) if fiducials is None else set(map(str, fiducials))
    clash = sorted(fid & set(pos))
    if clash:
        raise ValueError(f"listed both as positioner and fiducial: {clash}")
    has_ij = len(c) >= 6
    name = os.path.basename(out_path)
    hdr = ["Title: K-SPEC Fiber Configuration (lab test)",
           f"Authors: made by kspec_metrology.lab.make_fiber_table from {os.path.basename(src_table)}",
           "Table: Fiber Configuration",
           "=" * 80,
           f"Byte-by-byte Description of file: {name}",
           "-" * 80,
           " Bytes Format Units  Label     Explanations",
           "-" * 80,
           " 1- 3  A3     ---    ID           Hole ID",
           " 5-13  F9.4   mm     X            X coordinate",
           "15-23  F9.4   mm     Y            Y coordinate",
           "25-26  I2     ---    FiducialFlag [-9/1] Fiducial flag (1)"]
    if has_ij:
        hdr += ["28-29  I2     ---    Row          Row number (i)",
                "31-32  I2     ---    Column       Column number (j)"]
    hdr += ["-" * 80,
            "Note (1) : Flag as follows :",
            "    0 = Fiber positioner",
            "    1 = Fiducial",
            "    -9 = Neither (positioner not installed)",
            "-" * 80]
    npos = nfid = 0
    rows = []
    for k, hole in enumerate(ids):
        new = 0 if hole in pos else (1 if hole in fid else -9)
        npos += new == 0; nfid += new == 1
        line = f"{hole:<3} {float(tab[c[1]][k]):9.4f} {float(tab[c[2]][k]):9.4f} {new:>2}"
        if has_ij:
            line += f" {int(tab[c[4]][k]):>2} {int(tab[c[5]][k]):>2}"
        rows.append(line)
    with open(out_path, 'w') as ff:
        ff.write('\n'.join(hdr + rows) + '\n')
    return npos, nfid


def ab2xy(xorigin, yorigin, alpha_deg, beta_deg, arm1, arm2):
    """positioner 출력 각도 (alpha, beta) [deg, json 과 같은 부호 규약] -> fiber 끝 (x, y) [mm].
    mtlcal 은 계산각 (theta, phi) 를 부호를 뒤집어 내보내므로 (출력 = -계산각) 여기서 되돌린다."""
    theta = -np.asarray(alpha_deg, float)/rtod
    phi = -np.asarray(beta_deg, float)/rtod
    if cfg.CLOCKWISE:
        theta, phi = -theta, -phi
    x = xorigin + arm1*np.cos(theta) + arm2*np.cos(theta + phi)
    y = yorigin + arm1*np.sin(theta) + arm2*np.sin(theta + phi)
    return x, y


def xy2ab(xorigin, yorigin, x, y, arm1, arm2):
    """(x, y) [mm] -> 출력 각도 (alpha, beta) [deg]. ab2xy 의 역이고 mtlcal 과 같은 계산이다."""
    x, y = np.atleast_1d(x), np.atleast_1d(y)
    xo, yo = np.broadcast_to(xorigin, x.shape), np.broadcast_to(yorigin, x.shape)
    a1, a2 = np.broadcast_to(arm1, x.shape), np.broadcast_to(arm2, x.shape)
    th = np.empty(x.shape); ph = np.empty(x.shape)
    for i in range(x.size):
        th.flat[i], ph.flat[i] = find_angle_double_method2_select_elbow(
            xo.flat[i], yo.flat[i], x.flat[i], y.flat[i], a1.flat[i], a2.flat[i],
            elbow=cfg.ELBOW, clockwise=cfg.CLOCKWISE)
    return -th*rtod, -ph*rtod


def read_spec(path):
    """'ID  v1  v2' 세 칸짜리 텍스트를 {ID: (v1, v2)} 로 읽는다. '#' 뒤는 주석."""
    out = {}
    with open(path) as ff:
        for n, line in enumerate(ff, 1):
            line = line.split('#')[0].strip()
            if not line:
                continue
            p = line.replace(',', ' ').split()
            if len(p) != 3:
                raise ValueError(f"{path}:{n}: 'ID v1 v2' 세 칸이어야 한다: {line!r}")
            if p[0] in out:
                raise ValueError(f"{path}:{n}: {p[0]} 가 중복이다")
            out[p[0]] = (float(p[1]), float(p[2]))
    return out


def make_target_file(spec_path, out_path, kind='xy', tile='lab', use_measured=None,
                     zenith_angle=0.0):
    """목표 텍스트 파일 -> mtlcal 이 읽는 target json (object.info 형식).

    spec_path : 'ID x y' (kind='xy', mm) 또는 'ID alpha beta' (kind='ab', deg; json 과
                같은 부호 규약) 세 칸 텍스트. 표의 positioner 전부가 있어야 한다.
    각도는 positioner 회전 중심(실측값이 켜져 있으면 실측)과 arm 길이로 (x, y) 로 바꾼다.
    반환: (xp, yp) 표 순서."""
    from kspec_metrology.analysis.mtlcal import load_configuration
    spec = read_spec(spec_path)
    ids = list(cfg.fiber_ids())
    missing = [i for i in ids if i not in spec]
    extra = [i for i in spec if i not in ids]
    if missing:
        raise ValueError(f"{spec_path}: positioner {missing} 의 목표가 없다 "
                         f"(표 {cfg.FIBER_TABLE_PATH} 의 positioner 전부가 있어야 한다)")
    if extra:
        raise ValueError(f"{spec_path}: 표에 positioner 로 없는 ID: {extra}")
    v1 = np.array([spec[i][0] for i in ids]); v2 = np.array([spec[i][1] for i in ids])
    if kind == 'xy':
        xp, yp = v1, v2
    elif kind == 'ab':
        # 회전 중심: load_configuration 이 실측/설계를 가려서 준다 (target 은 아직 없으니 더미)
        tmp = out_path + '.tmp'
        json.dump(dict(tile_id=tile, xp=[0.]*len(ids), yp=[0.]*len(ids)), open(tmp, 'w'))
        try:
            xo, yo, _, _, _, nfib = load_configuration(tmp, apply_zenith_offset=False,
                                                       use_measured=use_measured)
        finally:
            os.remove(tmp)
        a1, a2 = cfg.arm_lengths(use_measured)
        xp, yp = ab2xy(xo[:nfib], yo[:nfib], v1, v2, a1, a2)
    else:
        raise ValueError("kind must be 'xy' or 'ab'")
    data = dict(tile_id=tile, inst='', func='lab', ra=[], dec=[],
                xp=[float(v) for v in xp], yp=[float(v) for v in yp],
                zenith_angle=float(zenith_angle))
    with open(out_path, 'w') as ff:
        json.dump(data, ff, indent=1)
    return xp, yp


# ============================================================================ 그림
def _binned(im, n):
    ny, nx = im.shape
    im = im[:ny - ny % n, :nx - nx % n]
    return im.reshape(ny // n, n, nx // n, n).mean(axis=(1, 3))


def load_frames(data_dir, head, nexposure=1):
    """findpeak 과 같은 방향([::-1,:]) 으로 읽어 평균한 이미지."""
    im = None
    for i in range(nexposure):
        f = fits.getdata(image_path(data_dir, head, i)).astype(float)[::-1, :] / nexposure
        im = f if im is None else im + f
    return im


def chip2pix(xchip, ychip, shape):
    """findpeak 의 칩 좌표 [mm] -> 픽셀 (u, v)."""
    ny, nx = shape
    return xchip / PIX_MM + (nx - 1) / 2., ychip / PIX_MM + (ny - 1) / 2.


def show_image(im, n=8, ax=None, title=None, vmax_pct=99.5):
    """n x n 으로 줄인 전체 영상."""
    import matplotlib.pyplot as plt
    small = _binned(im, n)
    if ax is None:
        _, ax = plt.subplots(figsize=(13, 10))
    lo, hi = np.percentile(small, [1, vmax_pct])
    ax.imshow(small, cmap='gray', vmin=lo, vmax=hi,
              extent=(-.5, small.shape[1]*n - .5, small.shape[0]*n - .5, -.5))
    ax.set_xlabel('u [px]'); ax.set_ylabel('v [px]')
    if title:
        ax.set_title(title)
    return ax


def show_detections(im, xobs, yobs, labels=None, flag=None, n=8, ax=None, title=None, r=60):
    """검출된 spot (findpeak 의 칩 좌표) 을 영상 위에 표시. flag 가 True 인 것은 주황."""
    import matplotlib.pyplot as plt
    ax = show_image(im, n, ax, title)
    u, v = chip2pix(np.asarray(xobs), np.asarray(yobs), im.shape)
    flag = np.zeros(u.size, bool) if flag is None else np.asarray(flag, bool)
    for k in range(u.size):
        c = plt.Circle((u[k], v[k]), r, fc='none', ec='tab:orange' if flag[k] else 'c', lw=0.8)
        ax.add_artist(c)
        if labels is not None:
            ax.text(u[k], v[k] - 1.6*r, str(labels[k]), color='y', fontsize=7, ha='center')
    return ax


def show_matching(x, y, xfocal, yfocal, fid_flag, ids=None, ax=None, ref_um=200.):
    """초점면에서 목표(설계) 와 측정 위치: 화살표 = 측정 - 목표 (과장), fiducial 은 주황."""
    import matplotlib.pyplot as plt
    if ax is None:
        _, ax = plt.subplots(figsize=(9, 8.5))
    dx, dy = xfocal - x, yfocal - y
    k = 12. / (ref_um*1e-3)
    q = None
    for f, c, lab in ((~fid_flag, 'k', 'fiber'), (fid_flag, 'tab:orange', 'fiducial')):
        if f.any():
            q = ax.quiver(x[f], y[f], dx[f]*k, dy[f]*k, angles='xy', scale_units='xy', scale=1, color=c, width=0.004)
        ax.plot(x[f], y[f], 'o', mfc='none', color=c, ms=5, label=lab)
    if q is not None:
        ax.quiverkey(q, 0.85, 1.03, 12., f"{ref_um:.0f} um", labelpos='E')
    if ids is not None:
        for xi, yi, i in zip(x, y, ids):
            ax.text(xi, yi - 3, str(i), fontsize=6, ha='center', va='top')
    ax.set_aspect('equal'); ax.set_xlabel('x [mm]'); ax.set_ylabel('y [mm]'); ax.legend(loc='lower left')
    ax.set_title('measured - target (arrows exaggerated)')
    return ax


def show_offsets(offset_file, ax=None):
    """mtlcal 이 저장한 _offset.txt 를 막대로."""
    import matplotlib.pyplot as plt
    ids = np.loadtxt(offset_file, usecols=0, dtype=str)
    d = np.loadtxt(offset_file, usecols=(1, 2, 3, 4), ndmin=2)
    if ax is None:
        _, ax = plt.subplots(2, 1, figsize=(max(8, 0.25*len(ids)), 7), sharex=True)
    k = np.arange(len(ids))
    ax[0].bar(k - 0.2, d[:, 0], 0.4, label='dx'); ax[0].bar(k + 0.2, d[:, 1], 0.4, label='dy')
    ax[0].set_ylabel('position error [um]'); ax[0].legend(); ax[0].axhline(0, c='gray', lw=.5)
    ax[1].bar(k - 0.2, d[:, 2], 0.4, label='dalpha'); ax[1].bar(k + 0.2, d[:, 3], 0.4, label='dbeta')
    ax[1].set_ylabel('rotation to apply [deg]'); ax[1].legend(); ax[1].axhline(0, c='gray', lw=.5)
    ax[1].set_xticks(k); ax[1].set_xticklabels(ids, rotation=90, fontsize=7)
    return ax


# ============================================================================ 단계별 확인
def inspect_trial(run, itrial, n=8):
    """MetrologyRun 이 분석한 trial 의 이미지를 같은 설정으로 단계별로 다시 돌려 그림으로
    보여 준다 (검출 -> fiducial 매칭/변환 fit -> fiber 매칭). json 은 건드리지 않는다.
    반환: dict(x, y, fid_flag, xfocal, yfocal, coeff, ids)"""
    import matplotlib.pyplot as plt
    from kspec_metrology.analysis.mtlcal import (load_configuration, measure_fiducials,
                                                 measure_fibers)
    xo, yo, x, y, fid_flag, nfib = load_configuration(
        run.target_file, apply_pinhole_offset=run.apply_pinhole_offset, use_measured=run.use_measured)
    ids = list(cfg.fiber_ids()) + [i for i, f in zip(_table_ids(), _table_flags()) if f == 1
                                   and i not in cfg.FIDUCIAL_EXCLUDE]
    pk = dict(threshold=run.threshold, nwindow=run.nwindow, niter_recenter=run.niter_recenter,
              finder=run.finder, finder_opts=run.finder_opts, background=run.background,
              background_opts=run.background_opts)
    head = run.image_head(itrial)
    heads = run.fid_image_heads(itrial)
    im = load_frames(run.data_dir, head, run.nexposure)
    if heads:
        xo_f, yo_f, choice, th, (ox, oy) = measure_fiducials(
            heads, x, y, fid_flag, run.data_dir, run.fid_nexposure, m_init=run.magnification,
            sat_level=run.sat_level, fid_min_peak=run.fid_min_peak, **pk)
        fid_ids = [i for i, f in zip(ids, fid_flag) if f]
        for k, h in enumerate(heads):
            imf = load_frames(run.data_dir, h, run.fid_nexposure)
            show_detections(imf, xo_f, yo_f, labels=fid_ids, flag=(choice != k), n=n,
                            title=f'fiducial set {k} ({h}): cyan = used from this set, orange = taken from another set')
        xf, yf, _, _, coeff = fitdistortion(
            x[fid_flag], y[fid_flag], np.ones(int(fid_flag.sum()), bool), xo_f, yo_f,
            np.arange(int(fid_flag.sum())), th, ox, oy, model=run.model, m_init=run.magnification)
        xfocal, yfocal, dist, (xobs, yobs) = measure_fibers(
            head, x, y, nfib, coeff, run.data_dir, run.nexposure, **pk)
        xfocal[fid_flag], yfocal[fid_flag] = xf, yf
        show_detections(im, xobs, yobs, labels=ids, flag=fid_flag, n=n,
                        title=f'fiber image ({head}): matched peaks (orange = fiducial positions, not used)')
    else:
        _, xobs, yobs, peak, nsat = findpeak(x.size, data_dir=run.data_dir, head=head,
                                             nexposure=run.nexposure, mode=run.mode, x=x, y=y,
                                             sat_level=run.sat_level, ReturnNsat=True, **pk)
        imatch, th, (ox, oy) = matchfiber(x, y, xobs, yobs, m_init=run.magnification)
        xfocal, yfocal, _, _, coeff = fitdistortion(x, y, fid_flag, xobs, yobs, imatch, th, ox, oy,
                                                    model=run.model, m_init=run.magnification)
        show_detections(im, xobs[imatch], yobs[imatch], labels=ids, flag=fid_flag, n=n,
                        title=f'{head}: matched peaks (orange = fiducial); saturated: {int((nsat[imatch] > 0).sum())}')
    show_matching(x, y, xfocal, yfocal, fid_flag, ids)
    r_fid = np.hypot(xfocal - x, yfocal - y)[fid_flag]*1e3
    r_fib = np.hypot(xfocal - x, yfocal - y)[~fid_flag]*1e3
    print(f"fiducial fit residual: rms {np.sqrt((r_fid**2).mean()):.1f} um, max {r_fid.max():.1f} um   |   "
          f"fiber - target: median {np.median(r_fib):.1f} um, max {r_fib.max():.1f} um")
    return dict(x=x, y=y, fid_flag=fid_flag, xfocal=xfocal, yfocal=yfocal, coeff=coeff, ids=ids)


def _table_ids():
    tab = ascii.read(cfg.FIBER_TABLE_PATH)
    return [str(v).strip() for v in tab[tab.colnames[0]]]


def _table_flags():
    tab = ascii.read(cfg.FIBER_TABLE_PATH)
    return [int(v) for v in tab[tab.colnames[3]]]
