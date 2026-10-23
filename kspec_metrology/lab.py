"""실험실 시험용 도우미: 목표 파일 만들기, fiber 표 줄이기, 단계별 확인 그림.

실험실에서는 목표가 관측 타일이 아니라 사람이 정한 (x, y) 나 (alpha, beta) 로
주어진다. 이 모듈이 그것을 mtlcal 이 읽는 target 파일(object.info 형식) 로
바꾸고, 꽂힌 positioner 만 남긴 fiber 표를 만들고, 각 단계의 결과를 그림으로
보여 준다. 분석 자체는 analysis/ 의 함수를 그대로 쓴다.

    measure            목표 파일 + 이미지 1-3장 -> fiber 위치 오차와 회전량 (한 번 측정)
    show_steps         measure 의 결과를 단계별 그림으로 (검출 / fiducial 매칭 / fiber 매칭)
    read_targets       object.info 형식(xp, yp 또는 각도) 목표 파일 읽기
    make_target        같은 각도 / 같은 offset / 개별 (x, y) 로 목표 파일 만들기
    check_peaks        이미지 한 장의 peak 검출 점검 (개수, 밝기, 쪼개짐)
    make_fiber_table   꽂힌 홀만 positioner 로 둔 Fiber_Configuration 표
    ab2xy, xy2ab       positioner 각도 <-> 초점면 좌표 (패키지 각도 규약 그대로)
    make_target_file   'ID x y' 또는 'ID alpha beta' 텍스트 -> target json
    inspect_trial      MetrologyRun 의 한 trial 을 단계별 그림으로
    show_image, show_detections, show_matching, show_offsets   확인 그림
"""

import json
import os
import re

import numpy as np
from astropy.io import ascii, fits

from kspec_metrology.analysis import fiber_config as cfg
from kspec_metrology.analysis.utils import find_angle_double_method2_select_elbow
from kspec_metrology.analysis.findpeak import findpeak
from kspec_metrology.analysis.matchfiber import matchfiber
from kspec_metrology.analysis.fitdistortion import fitdistortion
from kspec_metrology.naming import image_path
from kspec_metrology import naming

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
    sa, sb = cfg.ANGLE_SIGN
    theta = sa*np.asarray(alpha_deg, float)/rtod          # 부호 x 부호 = 1 이므로 곱해서 되돌린다
    phi = sb*np.asarray(beta_deg, float)/rtod
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
    sa, sb = cfg.ANGLE_SIGN
    return sa*th*rtod, sb*ph*rtod


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
    # test_260917.ipynb 와 같은 방향: origin='lower' (v 가 위로 커진다). 이미지는 findpeak 처럼
    # [::-1, :] 로 뒤집어 읽은 것이라 이렇게 그리면 배치표(초점면 x, y)와 바로 견줄 수 있다.
    ax.imshow(small, cmap='gray', vmin=lo, vmax=hi, origin='lower',
              extent=(-.5, small.shape[1]*n - .5, -.5, small.shape[0]*n - .5))
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
def _all_ids():
    """load_configuration 과 같은 순서의 홀 ID (fiber 표 순서 + fiducial)."""
    return list(cfg.fiber_ids()) + [i for i, f in zip(_table_ids(), _table_flags())
                                   if f == 1 and i not in cfg.FIDUCIAL_EXCLUDE]


def _inspect(data_dir, head, fid_heads, target_file, nexposure=1, fid_nexposure=None, mode='Raw',
             threshold=3e3, nwindow=40, niter_recenter=0, finder='sep', finder_opts=None,
             background='crop', background_opts=None, magnification=None, model='full',
             sat_level=65000., fid_min_peak=0., n_extra=0, use_measured=None,
             apply_pinhole_offset=False, cross_exclude=True, min_dist=40, fid_fallback_fiber=True, n=8):
    """mtlcal 과 같은 설정으로 단계별로 다시 분석하며 그림을 그린다. json 은 건드리지 않는다."""
    from kspec_metrology.analysis.mtlcal import (load_configuration, measure_fiducials,
                                                 measure_fibers)
    xo, yo, x, y, fid_flag, nfib = load_configuration(
        target_file, apply_pinhole_offset=apply_pinhole_offset, use_measured=use_measured)
    ids = _all_ids()
    # 끊어진 fiber 는 빼고 본다 (mtlcal 과 같다)
    dead = np.r_[cfg.dead_mask(), np.zeros(x.size-nfib, dtype=bool)]
    x, y, fid_flag = x[~dead], y[~dead], fid_flag[~dead]
    ids = [i for i, d in zip(ids, dead) if not d]
    nfib = int((~dead[:nfib]).sum())
    pk = dict(threshold=threshold, nwindow=nwindow, niter_recenter=niter_recenter, finder=finder,
              finder_opts=finder_opts, background=background, background_opts=background_opts,
              min_dist=min_dist)
    fid_nexposure = nexposure if fid_nexposure is None else fid_nexposure
    im = load_frames(data_dir, head, nexposure)
    if fid_heads:
        heads_all = list(fid_heads) + ([head] if fid_fallback_fiber else [])
        nexp_all = [fid_nexposure]*len(fid_heads) + ([nexposure] if fid_fallback_fiber else [])
        fb_all = [False]*len(fid_heads) + ([True] if fid_fallback_fiber else [])
        xo_f, yo_f, choice, th, (ox, oy) = measure_fiducials(
            heads_all, x, y, fid_flag, data_dir, fid_nexposure, set_nexposure=nexp_all, fallback=fb_all,
            m_init=magnification, sat_level=sat_level, fid_min_peak=fid_min_peak, n_extra=n_extra, **pk)
        fid_ids = [i for i, f in zip(ids, fid_flag) if f]
        for k, h in enumerate(heads_all):
            imf = load_frames(data_dir, h, nexp_all[k])
            show_detections(imf, xo_f, yo_f, labels=fid_ids, flag=(choice != k), n=n,
                            title=f'fiducial set {k} ({h}): cyan = used from this set, orange = taken from another set')
        xf, yf, _, _, coeff = fitdistortion(
            x[fid_flag], y[fid_flag], np.ones(int(fid_flag.sum()), bool), xo_f, yo_f,
            np.arange(int(fid_flag.sum())), th, ox, oy, model=model, m_init=magnification)
        xfocal, yfocal, dist, (xobs, yobs) = measure_fibers(
            head, x, y, nfib, coeff, data_dir, nexposure, n_extra=n_extra,
            exclude=(xo_f, yo_f) if cross_exclude else None, **pk)
        if cross_exclude:                       # fiber 자리를 지우고 fiducial 을 다시 (mtlcal 과 같다)
            xo_f, yo_f, choice, th, (ox, oy) = measure_fiducials(
                heads_all, x, y, fid_flag, data_dir, fid_nexposure, set_nexposure=nexp_all, fallback=fb_all,
                m_init=magnification, sat_level=sat_level, fid_min_peak=fid_min_peak, n_extra=n_extra,
                exclude=(xobs[:nfib], yobs[:nfib]), **pk)
            xf, yf, _, _, coeff = fitdistortion(
                x[fid_flag], y[fid_flag], np.ones(int(fid_flag.sum()), bool), xo_f, yo_f,
                np.arange(int(fid_flag.sum())), th, ox, oy, model=model, m_init=magnification)
            xfocal, yfocal, dist, (xobs, yobs) = measure_fibers(
                head, x, y, nfib, coeff, data_dir, nexposure, n_extra=n_extra,
                exclude=(xo_f, yo_f), **pk)
        xfocal[fid_flag], yfocal[fid_flag] = xf, yf
        show_detections(im, xobs, yobs, labels=ids, flag=fid_flag, n=n,
                        title=f'fiber image ({head}): matched peaks (orange = fiducial positions, not used)')
    else:
        npk = x.size if (mode != 'Raw' or n_extra == 0) else (x.size, x.size + n_extra)
        _, xobs, yobs, peak, nsat = findpeak(npk, data_dir=data_dir, head=head, nexposure=nexposure,
                                             mode=mode, x=x, y=y, sat_level=sat_level,
                                             ReturnNsat=True, **pk)
        imatch, th, (ox, oy) = matchfiber(x, y, xobs, yobs, m_init=magnification)
        xfocal, yfocal, _, _, coeff = fitdistortion(x, y, fid_flag, xobs, yobs, imatch, th, ox, oy,
                                                    model=model, m_init=magnification)
        show_detections(im, xobs[imatch], yobs[imatch], labels=ids, flag=fid_flag, n=n,
                        title=f'{head}: matched peaks (orange = fiducial); saturated: {int((nsat[imatch] > 0).sum())}')
    show_matching(x, y, xfocal, yfocal, fid_flag, ids)
    r_fid = np.hypot(xfocal - x, yfocal - y)[fid_flag]*1e3
    r_fib = np.hypot(xfocal - x, yfocal - y)[~fid_flag]*1e3
    print(f"fiducial fit residual: rms {np.sqrt((r_fid**2).mean()):.1f} um, max {r_fid.max():.1f} um   |   "
          f"fiber - target: median {np.median(r_fib):.1f} um, max {r_fib.max():.1f} um")
    return dict(x=x, y=y, fid_flag=fid_flag, xfocal=xfocal, yfocal=yfocal, coeff=coeff, ids=ids)


def inspect_trial(run, itrial, n=8):
    """MetrologyRun 이 분석한 trial 을 같은 설정으로 단계별로 다시 돌려 그림으로 보여 준다."""
    return _inspect(run.data_dir, run.image_head(itrial), run.fid_image_heads(itrial), run.target_file,
                    nexposure=run.nexposure, fid_nexposure=run.fid_nexposure, mode=run.mode,
                    threshold=run.threshold, nwindow=run.nwindow, niter_recenter=run.niter_recenter,
                    finder=run.finder, finder_opts=run.finder_opts, background=run.background,
                    background_opts=run.background_opts, magnification=run.magnification,
                    model=run.model, sat_level=run.sat_level, fid_min_peak=run.fid_min_peak,
                    n_extra=run.n_extra, use_measured=run.use_measured,
                    apply_pinhole_offset=run.apply_pinhole_offset, n=n)


# ============================================================================ 목표 파일 만들기
def make_target(out_path, angles=None, offset=None, xy=None, tile='lab', use_measured=None):
    """target 파일(object.info 형식) 을 만든다. 셋 중 하나를 준다.
        angles=(alpha, beta)  모든 positioner 를 같은 각도로 [deg, json 부호 규약]
        offset=(dx, dy)       모든 fiber 끝을 각자의 회전 중심에서 같은 (dx, dy) 만큼 [mm]
        xy=(xp, yp)           positioner 마다 다른 목표 [mm] (표 순서)
    각도 파일은 'alpha'/'beta' 키로 저장되어 read_targets 가 그대로 읽는다.
    반환: 저장한 dict."""
    ids = list(cfg.fiber_ids()); n = len(ids)
    if angles is not None:
        data = dict(tile_id=tile, alpha=[float(angles[0])]*n, beta=[float(angles[1])]*n)
    elif offset is not None:
        from kspec_metrology.analysis.mtlcal import load_configuration
        tmp = out_path + '.tmp'
        json.dump(dict(tile_id=tile, xp=[0.]*n, yp=[0.]*n), open(tmp, 'w'))
        try:
            xo, yo, _, _, _, nfib = load_configuration(tmp, apply_zenith_offset=False, use_measured=use_measured)
        finally:
            os.remove(tmp)
        data = dict(tile_id=tile, xp=[float(v) for v in xo[:nfib] + offset[0]],
                    yp=[float(v) for v in yo[:nfib] + offset[1]])
    elif xy is not None:
        xp, yp = np.asarray(xy[0], float), np.asarray(xy[1], float)
        if xp.size != n or yp.size != n:
            raise ValueError(f"xy must have {n} entries (table order {ids})")
        data = dict(tile_id=tile, xp=[float(v) for v in xp], yp=[float(v) for v in yp])
    else:
        raise ValueError("give angles=(alpha, beta), offset=(dx, dy) or xy=(xp, yp)")
    with open(out_path, 'w') as ff:
        json.dump(data, ff, indent=1)
    return data


# ============================================================================ 검출 점검
def check_peaks(fits_path, threshold=3e3, min_dist=40, finder='sep', finder_opts=None, background='crop',
                nwindow=40, exclude=None, exclude_r=60, n=8, crop=60, show_crops=24):
    """이미지 한 장에서 peak 검출만 해 보고 (findpeak 과 같은 길), 개수·밝기·서로 가까운 쌍을
    표와 그림으로 보여 준다. 한 spot 이 여러 peak 으로 쪼개지는지, threshold 가 맞는지,
    min_dist / finder_opts 를 어떻게 잡을지 정할 때 쓴다.

    문턱값 탐색 없이 threshold 그대로 한 번 돌린다.
    반환: pandas DataFrame(u, v, peak, nearest_px) -- nearest_px 는 가장 가까운 다른 검출까지 거리."""
    import matplotlib.pyplot as plt
    import pandas as pd
    from scipy.spatial import cKDTree
    d = os.path.dirname(os.path.abspath(fits_path))
    base = os.path.basename(fits_path)
    if not base.endswith('0.fits'):
        raise ValueError("파일 이름이 ...0.fits 로 끝나야 한다 (findpeak 의 이름 규칙)")
    head = base[:-len('0.fits')]
    im, xo, yo, pk = findpeak((1, 100000), data_dir=d, head=head, nexposure=1, threshold=threshold,
                              niter_max=1, min_dist=min_dist, finder=finder, finder_opts=finder_opts,
                              background=background, nwindow=nwindow, exclude=exclude, exclude_r=exclude_r)
    u, v = chip2pix(xo, yo, im.shape)
    dd = cKDTree(np.c_[u, v]).query(np.c_[u, v], k=2)[0][:, 1] if u.size > 1 else np.full(u.size, np.nan)
    T = pd.DataFrame(dict(u=u, v=v, peak=pk, nearest_px=dd)).sort_values('peak', ascending=False).reset_index(drop=True)
    close = T[T.nearest_px < 3*min_dist]
    print(f"{base}: {len(T)} peaks at threshold {threshold:.0f} (min_dist {min_dist} px); "
          f"peak {T.peak.min():.0f}..{T.peak.max():.0f}; {len(close)} within {3*min_dist} px of another peak "
          f"({'split spots?' if len(close) else 'ok'})")
    show_detections(im, xo, yo, flag=(dd < 3*min_dist), n=n, r=max(60, min_dist),
                    title=f'{base}: {len(T)} peaks (orange = closer than {3*min_dist} px to another)')
    k = min(show_crops, len(T))
    if k:
        ncol = 8; nrow = int(np.ceil(k/ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(2.2*ncol, 2.3*nrow), squeeze=False)
        order = np.argsort(T.nearest_px.values)[:k]            # 가까운 쌍부터
        for a, i in zip(axs.ravel(), order):
            uu, vv = int(round(T.u[i])), int(round(T.v[i]))
            sub = im[max(vv-crop, 0):vv+crop+1, max(uu-crop, 0):uu+crop+1]
            a.imshow(np.log10(np.clip(sub - np.median(sub), 3, None)), cmap='inferno', origin='lower')
            a.set_title(f'peak {T.peak[i]:.0f}  near {T.nearest_px[i]:.0f}px', fontsize=7); a.set_xticks([]); a.set_yticks([])
        for a in axs.ravel()[k:]:
            a.axis('off')
        fig.tight_layout()
    return T


# ============================================================================ 한 번 측정
def read_targets(path):
    """목표 파일 -> (xp, yp, kind). object.info 와 같은 json 이고 셋 중 하나를 담는다.
        {"xp": [...], "yp": [...]}            초점면 목표 [mm]            (kind 'xy')
        {"alpha": [...], "beta": [...]}       positioner 각도 [deg]      (kind 'ab')
        {"C4_a": .., "C4_b": .., ...}         각도 json (mtlcal 출력 형식) (kind 'ab')
    배열은 fiber 표의 positioner 순서여야 한다. 각도는 json 과 같은 부호 규약."""
    data = json.load(open(path))
    low = {str(k).lower(): v for k, v in data.items()}
    ids = list(cfg.fiber_ids())
    if 'xp' in low and 'yp' in low:
        xp, yp, kind = np.asarray(low['xp'], float), np.asarray(low['yp'], float), 'xy'
    elif 'alpha' in low and 'beta' in low:
        xp, yp, kind = np.asarray(low['alpha'], float), np.asarray(low['beta'], float), 'ab'
    elif all(f"{i}_a" in data and f"{i}_b" in data for i in ids):
        xp = np.array([data[f"{i}_a"] for i in ids], float)
        yp = np.array([data[f"{i}_b"] for i in ids], float)
        kind = 'ab'
    else:
        raise ValueError(f"{path}: 'xp'/'yp', 'alpha'/'beta', 또는 '{{ID}}_a'/'{{ID}}_b' 키가 필요하다")
    if xp.size != len(ids) or yp.size != len(ids):
        raise ValueError(f"{path}: 값이 {xp.size}/{yp.size} 개인데 표의 positioner 는 {len(ids)} 개다 "
                         f"({cfg.FIBER_TABLE_PATH})")
    return xp, yp, kind


def next_trial(out_dir, tile):
    """out_dir 에 있는 {tile}_MetrologyTrial_{k}.json 중 가장 큰 k + 1 (없으면 1)."""
    k = 0
    for f in os.listdir(out_dir) if os.path.isdir(out_dir) else []:
        m = re.fullmatch(re.escape(tile) + r'_MetrologyTrial_(\d+)\.json', f)
        if m:
            k = max(k, int(m.group(1)))
    return k + 1


def measure(target_file, fiber_image, fid_images=(), out_dir='./lab_out', tile='lab', itrial=None,
            magnification=None, model='full', nwindow=40, threshold=3e3, niter_recenter=0,
            finder='sep', finder_opts=None, background='crop', background_opts=None,
            sat_level=65000., fid_min_peak=0., n_extra=0, use_measured=None, mode='Raw',
            cross_exclude=True, min_dist=40, fid_fallback_fiber=True):
    """한 배치를 한 번 잰다: 목표 파일과 이미지 1-3장 -> fiber 위치 오차와 positioner 회전량.

    target_file : object.info 형식 (read_targets 참고; xp/yp 또는 각도)
    fiber_image : fiber 가 찍힌 FITS 경로 (fiducial 도 같이 찍혔으면 이것 한 장으로 충분)
    fid_images  : fiducial 전용 노출 FITS 경로 0-2 개. 비면 fiber_image 에서 fiducial 도 찾고,
                  하나면 그것으로, 둘이면 fiducial 마다 포화되지 않은 밝은 쪽을 쓴다.
    itrial      : 몇 번째 측정인지. None 이면 out_dir 에 있는 json 다음 번호 (처음이면 1).
                  누적 각도 = 직전 trial json + 이번 회전량 이므로, positioner 를 옮기고 다시
                  찍은 이미지로 같은 out_dir/tile 로 부르면 각도가 누적된다. trial 0 json
                  (목표 각도) 은 처음 부를 때 함께 만든다.
    나머지 인자는 mtlcal 과 같다. 이미지는 out_dir/data 에 패키지 이름 규칙으로 링크(복사)되고
    결과는 out_dir 에 {tile}_MetrologyTrial_{itrial}.json (누적 각도) 과 _offset.txt 로 남는다.

    반환 dict: itrial, ids, dx, dy [mm] (fiber + fiducial, 측정 - 목표), angle_rot (deg, fiber
    순서로 alpha, beta 번갈아), angle_cum, json, offsets, 그리고 show_steps 에 넘길 설정."""
    import shutil
    out_dir = os.path.abspath(out_dir)
    data_dir = os.path.join(out_dir, 'data')
    os.makedirs(data_dir, exist_ok=True)
    if itrial is None:
        itrial = next_trial(out_dir, tile)

    # 목표: 각도면 (x, y) 로 바꿔 표준 target json 으로 저장
    xp, yp, kind = read_targets(target_file)
    if kind == 'ab':
        from kspec_metrology.analysis.mtlcal import load_configuration
        tmp = os.path.join(out_dir, '_tmp_target.json')
        json.dump(dict(tile_id=tile, xp=[0.]*xp.size, yp=[0.]*yp.size), open(tmp, 'w'))
        try:
            xo, yo, _, _, _, nfib = load_configuration(tmp, apply_zenith_offset=False,
                                                       use_measured=use_measured)
        finally:
            os.remove(tmp)
        a1, a2 = cfg.arm_lengths(use_measured)
        xp, yp = ab2xy(xo[:nfib], yo[:nfib], xp, yp, a1, a2)
    tfile = os.path.join(out_dir, f'{tile}_target.info')
    json.dump(dict(tile_id=tile, func='lab', xp=[float(v) for v in xp], yp=[float(v) for v in yp],
                   zenith_angle=0.0), open(tfile, 'w'), indent=1)

    def place(src, dst):
        if os.path.lexists(dst):
            os.remove(dst)
        try:
            os.symlink(os.path.abspath(src), dst)
        except OSError:
            shutil.copy2(src, dst)

    if itrial == 1:
        from kspec_metrology.analysis.mtlcal import save_initial_angles, save_angles
        if kind == 'ab':
            # 각도로 준 목표는 그 값 그대로 trial 0 으로 (xy 로 갔다 오면 -330 처럼 360 도 돌아간 값이 된다)
            a = np.empty(2*xp.size); a[0::2], a[1::2] = read_targets(target_file)[:2]
            save_angles(out_dir, tile, 0, a)
        else:
            save_initial_angles(out_dir, tile, tfile, use_measured=use_measured)   # trial 0 = 목표 각도
    head = naming.image_head(tile, itrial)
    place(fiber_image, image_path(data_dir, head, 0))
    fid_heads = []
    for k, f in enumerate(fid_images or ()):
        h = naming.fid_image_head(tile, itrial, k)
        place(f, image_path(data_dir, h, 0))
        fid_heads.append(h)

    from kspec_metrology.analysis.mtlcal import mtlcal
    params = dict(mode=mode, threshold=threshold, nwindow=nwindow, niter_recenter=niter_recenter,
                  finder=finder, finder_opts=finder_opts, background=background,
                  background_opts=background_opts, magnification=magnification, model=model,
                  sat_level=sat_level, fid_min_peak=fid_min_peak, n_extra=n_extra,
                  use_measured=use_measured, cross_exclude=cross_exclude, min_dist=min_dist,
                  fid_fallback_fiber=fid_fallback_fiber)
    dx, dy, angle_rot, angle_cum = mtlcal(data_dir=data_dir, head=head, nexposure=1,
                                          fid_head=fid_heads or None, target_file=tfile,
                                          json_dir=out_dir, target_name=tile, itrial=itrial, **params)
    return dict(itrial=itrial, ids=_all_ids(), nfib=xp.size, dx=dx, dy=dy, angle_rot=angle_rot, angle_cum=angle_cum,
                json=naming.json_path(out_dir, tile, itrial), offsets=naming.offset_path(out_dir, tile, itrial),
                data_dir=data_dir, head=head, fid_heads=fid_heads, target_file=tfile, params=params)


def show_steps(res, n=8):
    """measure 의 결과를 단계별 그림으로: fiducial 검출(세트별), fiber 매칭, 초점면 잔차."""
    return _inspect(res['data_dir'], res['head'], res['fid_heads'], res['target_file'], n=n, **res['params'])


def history(out_dir, tile):
    """out_dir 의 trial 별 offset 표를 읽어 fiber 오차의 최대/중앙값 [um] 표로 (pandas DataFrame)."""
    import pandas as pd
    rows = []
    for k in range(1, next_trial(out_dir, tile)):
        f = naming.offset_path(out_dir, tile, k)
        if not os.path.exists(f):
            continue
        d = np.loadtxt(f, usecols=(1, 2), ndmin=2)
        r = np.hypot(d[:, 0], d[:, 1]); r = r[np.isfinite(r)]
        rows.append(dict(trial=k, err_max_um=r.max(), err_median_um=np.median(r),
                         json=os.path.basename(naming.json_path(out_dir, tile, k))))
    return pd.DataFrame(rows).set_index('trial')


def offsets_table(res):
    """measure 결과를 fiber 별 표(pandas 가 있으면 DataFrame, 없으면 리스트)로."""
    nfib = res['nfib']
    rows = [(i, res['dx'][k]*1e3, res['dy'][k]*1e3, res['angle_rot'][2*k], res['angle_rot'][2*k+1],
             res['angle_cum'][2*k], res['angle_cum'][2*k+1])
            for k, i in enumerate(res['ids'][:nfib])]
    cols = ('ID', 'dx_um', 'dy_um', 'dalpha_deg', 'dbeta_deg', 'alpha_deg', 'beta_deg')
    try:
        import pandas as pd
        return pd.DataFrame(rows, columns=cols).set_index('ID')
    except ImportError:
        return [dict(zip(cols, r)) for r in rows]


def _table_ids():
    tab = ascii.read(cfg.FIBER_TABLE_PATH)
    return [str(v).strip() for v in tab[tab.colnames[0]]]


def _table_flags():
    tab = ascii.read(cfg.FIBER_TABLE_PATH)
    return [int(v) for v in tab[tab.colnames[3]]]
