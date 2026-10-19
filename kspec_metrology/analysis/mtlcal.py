import json
import os

import numpy as np
from astropy.io import ascii
from scipy.spatial import cKDTree

from kspec_metrology.analysis.findpeak import findpeak
from kspec_metrology.analysis.matchfiber import matchfiber
from kspec_metrology.analysis.fitdistortion import fitdistortion
from kspec_metrology.analysis.utils import find_angle_double_method2_select_elbow, transform
from kspec_metrology.analysis.zenith_offset import zenith_offset
from kspec_metrology.analysis.target_info import read_target, zenith_angle
from kspec_metrology.analysis import fiber_config as cfg
from kspec_metrology import naming
from kspec_metrology.logging.log import get_logger

rtod = 180. / np.pi

# 출력 json의 key 형식: {hole ID}_{arm}. alpha arm은 'a', beta arm은 'b'.
# 예) hole D0 -> "D0_a" (alpha), "D0_b" (beta)
ALPHA_SUFFIX = 'a'
BETA_SUFFIX = 'b'


def load_configuration(target_file=None, apply_zenith_offset=True,
                       apply_pinhole_offset=False, use_measured=None):
    """설계 좌표와 목표 위치를 읽는다.

    target 파일의 목표 위치 (xp, yp)는 천정거리 0도 기준이라, 관측 고도가 함께
    적혀 있으면 ADC 광학계의 sky -> focal plane 매핑 변화만큼 옮겨준다
    (zenith_offset 참고). 고도 정보가 없으면 보정 없이 그대로 쓴다.
    fiducial은 초점면에 박혀 있는 물체라 보정하지 않는다.

    apply_pinhole_offset=True 면 fiducial 기준 좌표에 핀홀이 홀 중심에서 치우친
    양(dx, dy)을 더한다 (fiber_config.pinhole_offsets). 이 값이 fitdistortion의
    기준이 되므로 맞게 넣으면 그만큼 왜곡 fit이 정확해진다.
    xorigin, yorigin은 홀 설계 좌표 그대로 둔다.

    use_measured 면 설계 좌표와 arm 길이를 실측값으로 바꾼다
    (fiber_config.measured_geometry). None이면 fiber_config.USE_MEASURED를
    따르고 그 기본값은 True다 (KSPEC_USE_MEASURED=0 으로 끈다). 실측 위치가
    있는 fiducial에는 핀홀 오프셋을 더하지 않는다 (이중 계산이 되므로).
    실측 파일이 없으면, 인자로 명시해 켠 경우에만 예외를 내고 아니면 경고만
    남기고 설계값을 쓴다. 실측 좌표가 설계 좌표에서 MEASURED_MAX_DEV 보다 멀면
    다른 plate 의 파일로 보고 멈춘다.

    기본값이 False인 이유: 핀홀 측정값은 plate마다 다시 재야 하는데, 표를 바꿔도
    fiducial ID가 상당수 겹치는 탓에 옛 plate의 값이 조용히 적용될 수 있다.
    지금 들어 있는 파일은 d = 16.8 mm plate 기준이다. 새로 측정해 파일을 갈아
    끼운 뒤 True로 켜면 된다.

    반환:
        xorigin, yorigin : positioner 회전 중심 (fiber -> fiducial 순)
        x, y             : 목표 위치. fiber는 target 파일 값, fiducial은 설계 값
        fid_flag         : fiducial이면 True
        nfib             : fiber 개수 (배열 앞쪽 nfib개가 fiber)
    """
    log = get_logger()

    tab = ascii.read(cfg.FIBER_TABLE_PATH)
    tab.rename_columns(tab.colnames[:4], ["ID", "X", "Y", "FiducialFlag"])

    ids = cfg.fiber_ids()
    nfib = len(ids)

    fid_ids = [i for i in tab["ID"][tab["FiducialFlag"] == 1]
               if i not in cfg.FIDUCIAL_EXCLUDE]
    pick_ids = list(dict.fromkeys(list(ids) + list(fid_ids)))

    log.info("Fiber table %s: %d positioners + %d fiducials",
             cfg.FIBER_TABLE_PATH, nfib, len(fid_ids))

    xy_map = {r["ID"]: (r["X"], r["Y"]) for r in tab}
    xy = np.array([xy_map[i] for i in pick_ids], dtype=float)
    xorigin, yorigin = xy[:, 0], xy[:, 1]

    #---설계 좌표를 실측값으로 교체------------------------------------------
    # positioner는 회전 중심, fiducial은 핀홀 위치다. 안 잰 홀은 설계값을 쓴다.
    measured = {}
    if cfg._use_measured(use_measured):
        measured = cfg.measured_geometry()
        if not measured:
            msg = ("use_measured가 켜져 있는데 실측 기하 파일이 없습니다: "
                   f"{cfg.MEASURED_GEOMETRY_PATH}")
            if use_measured:                     # 인자로 명시해 켠 경우만 멈춘다
                raise FileNotFoundError(msg)
            log.warning("%s -- 설계값을 씁니다", msg)

        nhit = 0
        for k, hole in enumerate(pick_ids):
            m = measured.get(hole)
            if m is None or m['x'] is None or m['y'] is None:
                continue
            dev = float(np.hypot(m['x']-xorigin[k], m['y']-yorigin[k]))
            if dev > cfg.MEASURED_MAX_DEV:
                raise ValueError(
                    f"실측 기하 {cfg.MEASURED_GEOMETRY_PATH}: {hole} 의 실측 좌표가 "
                    f"설계 좌표에서 {dev:.2f} mm 떨어져 있습니다 (한계 "
                    f"{cfg.MEASURED_MAX_DEV} mm). 다른 plate 의 실측 파일이거나 "
                    f"표({cfg.FIBER_TABLE_PATH})가 다릅니다. 맞는 파일로 바꾸거나 "
                    "KSPEC_USE_MEASURED=0 (use_measured=False) 로 끄십시오.")
            xorigin[k], yorigin[k] = m['x'], m['y']
            nhit += 1

        unknown = [h for h, m in measured.items() if h not in set(pick_ids)
                   and (m['x'] is not None or m['arm1'] is not None)]
        log.info("Measured geometry %s: %d/%d holes replaced",
                 cfg.MEASURED_GEOMETRY_PATH, nhit, len(pick_ids))
        if unknown:
            log.warning("Measured file lists %d hole(s) not in the fiber "
                        "table; ignored: %s", len(unknown),
                        ', '.join(sorted(unknown)[:10]))

    x, y = np.copy(xorigin), np.copy(yorigin)
    fid_flag = np.zeros(x.size, dtype=bool)
    fid_flag[nfib:] = True

    #---Fiducial 핀홀 오프셋----------------------------------------------------
    if apply_pinhole_offset:
        off = cfg.pinhole_offsets()
        if not off:
            log.warning("No fiducial pinhole offsets found (%s); "
                        "using the hole centres as the fiducial reference",
                        os.path.basename(cfg.FIDUCIAL_PINHOLE_MAP_PATH))
        else:
            # fiducial에만 더한다. fiber 자리는 곧 target 파일 값으로 덮인다.
            def _done(h):
                m = measured.get(h)
                return m is not None and m['x'] is not None and m['y'] is not None

            # 실측 위치로 이미 바꾼 fiducial에는 오프셋을 또 더하지 않는다
            dxy = np.array([off[i] if (f and i in off and not _done(i))
                            else (0., 0.)
                            for i, f in zip(pick_ids, fid_flag)])
            x += dxy[:, 0]
            y += dxy[:, 1]

            nhit = int(np.any(dxy != 0., axis=1).sum())
            log.info("Fiducial pinhole offsets applied to %d/%d fiducials "
                     "(max %.0f um)", nhit, int(fid_flag.sum()),
                     np.hypot(*dxy.T).max()*1e3)

            unused = [i for i in off
                      if i not in {h for h, f in zip(pick_ids, fid_flag) if f}]
            if unused:
                log.warning("Pinhole map lists %d hole(s) that are not used as "
                            "fiducials; ignored: %s", len(unused),
                            ', '.join(sorted(unused)[:10]))

    target = read_target(target_file)
    xp = np.asarray(target['xp'], dtype=float)
    yp = np.asarray(target['yp'], dtype=float)

    # 개수가 다르면 예전에는 앞에서 잘라 쓰며 조용히 지나갔다. fiber 표를 바꾸고
    # target 파일을 안 바꾼 경우가 대표적인데, 엉뚱한 fiber에 엉뚱한 목표가
    # 붙으므로 바로 멈춘다.
    if xp.size != nfib or yp.size != nfib:
        raise ValueError(
            f"target 파일의 xp/yp 개수({xp.size}/{yp.size})가 fiber 표의 "
            f"positioner 개수({nfib})와 다릅니다.\n"
            f"  표    : {cfg.FIBER_TABLE_PATH}\n"
            f"  target: {target_file or cfg.TARGET_INFO_PATH}\n"
            "표를 바꿨다면 target 파일도 같은 순서로 다시 만들어야 합니다 "
            "(KSPEC_FIBER_TABLE로 표를 지정할 수 있습니다).")

    x[:nfib] = xp
    y[:nfib] = yp

    #---관측 고도에 따른 보정 (fiber만)-----------------------------------------
    if apply_zenith_offset:
        za = zenith_angle(target)
        if za is None:
            log.warning("No altitude in %s; skipping the zenith angle correction",
                        target_file or cfg.TARGET_INFO_PATH)
        else:
            dx, dy = zenith_offset(x[:nfib], y[:nfib], za)
            x[:nfib] += dx
            y[:nfib] += dy

            log.info("Zenith angle %.2f deg: target shifted by max %.1f um, "
                     "median %.1f um", za,
                     np.sqrt(dx**2 + dy**2).max()*1e3,
                     np.median(np.sqrt(dx**2 + dy**2))*1e3)

    return xorigin, yorigin, x, y, fid_flag, nfib


def _arm_angles(xorigin, yorigin, xtarget, ytarget, use_measured=None):
    """각 positioner가 (xtarget, ytarget)을 가리키는 (theta, phi)를 구한다."""
    arm1s, arm2s = cfg.arm_lengths(use_measured)

    theta = np.full(arm1s.size, np.nan)
    phi = np.full(arm1s.size, np.nan)
    for i in range(arm1s.size):
        if not (np.isfinite(xtarget[i]) and np.isfinite(ytarget[i])):
            continue                                            # 끊어진 fiber: 측정값 없음
        theta[i], phi[i] = find_angle_double_method2_select_elbow(
            xorigin[i], yorigin[i], xtarget[i], ytarget[i], arm1s[i], arm2s[i],
            elbow=cfg.ELBOW,
            clockwise=cfg.CLOCKWISE,
        )
    return theta, phi


def _wrap(theta, phi):
    """theta는 (-pi, pi], phi는 (-pi/2, pi/2] 범위로 접는다."""
    theta = theta % (2.*np.pi)
    theta = np.where(theta > np.pi, theta - 2.*np.pi, theta)

    phi = phi % np.pi
    phi = np.where(phi > np.pi/2., phi - np.pi, phi)

    return theta, phi


def _flatten_angles(theta, phi):
    """fiber 순서의 (theta, phi)를 (alpha, beta)가 번갈아 놓인 degree 배열로 편다.

    출력 부호는 positioner 제어 규약에 맞춰 뒤집는다 (출력 각도 = -계산각).
    """
    angles = np.empty(2*theta.size)
    angles[0::2] = -theta*rtod
    angles[1::2] = -phi*rtod

    return angles


def angle_dict(angles):
    """각도 배열을 {"B5_a": ..., "B5_b": ...} 형태로 만든다.

    key는 Fiber_Configuration 표의 hole ID에 arm 접미사를 붙인 것이고,
    순서는 FIBERS 표 순서(= target 파일의 fiber 순서)를 그대로 따른다.
    """
    data = {}
    for k, hole in enumerate(cfg.fiber_ids()):
        data[f"{hole}_{ALPHA_SUFFIX}"] = float(angles[2*k])
        data[f"{hole}_{BETA_SUFFIX}"] = float(angles[2*k+1])
    return data


def angle_array(data):
    """angle_dict의 역변환. json dict을 각도 배열로 되돌린다."""
    return np.array([data[f"{hole}_{suffix}"]
                     for hole in cfg.fiber_ids()
                     for suffix in (ALPHA_SUFFIX, BETA_SUFFIX)], dtype=float)


def json_path(json_dir, target_name, itrial):
    """각도 json 경로. 이미지 파일 이름과 같은 규칙을 쓴다 (naming 참고)."""
    return naming.json_path(json_dir, target_name, itrial)


def save_angles(json_dir, target_name, itrial, angles):
    """각도 배열을 trial별 json으로 저장한다."""
    os.makedirs(json_dir, exist_ok=True)

    path = json_path(json_dir, target_name, itrial)
    with open(path, 'w') as ff:
        json.dump(angle_dict(angles), ff, indent=2)

    return path


def save_initial_angles(json_dir, target_name, target_file=None,
                        use_measured=None):
    """관측 전에 목표 각도를 Trial 0 json으로 저장한다.

    이후 trial의 누적 각도는 이 파일에서 출발한다.
    """
    log = get_logger()

    xorigin, yorigin, x, y, _, _ = load_configuration(
        target_file, use_measured=use_measured)
    theta_true, phi_true = _arm_angles(xorigin, yorigin, x, y, use_measured)

    angles = _flatten_angles(theta_true, phi_true)
    path = save_angles(json_dir, target_name, 0, angles)
    log.info("Saved target angles to %s", path)

    return angles


def save_offsets(json_dir, target_name, itrial, dx, dy, angle_rot, nfib):
    """fiber 마다 남은 위치 오차와 이번 trial 회전량을 표로 저장한다.

    positioner 를 옮기는 쪽이 json(누적 각도) 대신 "얼마나 더 가야 하는지" 를
    바로 보고 싶을 때 쓴다. 열: ID, dx, dy [um] (측정 - 목표), dalpha, dbeta [deg]
    (이번 trial 회전량, json 과 같은 부호 규약).
    """
    os.makedirs(json_dir, exist_ok=True)
    path = naming.offset_path(json_dir, target_name, itrial)
    with open(path, 'w') as ff:
        ff.write("# fiber position error (measured - target) and this trial's rotation; "
                 "nan = excluded (broken) fiber\n")
        ff.write(f"# {'ID':<6} {'dx[um]':>9} {'dy[um]':>9} {'dalpha[deg]':>12} {'dbeta[deg]':>12}\n")
        for k, hole in enumerate(cfg.fiber_ids()[:nfib]):
            ff.write(f"  {hole:<6} {dx[k]*1e3:9.1f} {dy[k]*1e3:9.1f} "
                     f"{angle_rot[2*k]:12.4f} {angle_rot[2*k+1]:12.4f}\n")
    return path


def measure_fiducials(fid_heads, x, y, fid_flag, data_dir, nexposure,
                      m_init=None, sat_level=65000., fid_min_peak=0., n_extra=0, **pk):
    """fiducial 전용 노출 세트들에서 fiducial 의 카메라 좌표를 잰다.

    x, y 는 fiber 목표 + fiducial 설계 위치 전체, fid_flag 는 그중 fiducial.
    세트마다 peak 을 찾고 (fiber 가 보일 수도 있으므로 nfid..nfid+nfib 개 허용)
    설계 패턴에 매칭한다. fiber 가 보이면 (검출이 fiducial 수보다 많으면) fiber
    목표까지 넣어 전체 패턴으로 매칭한다 -- fiducial 30개만으로는 그 사이에
    fiber 가 빽빽이 있을 때 회전각이 정해지지 않는다. 세트가 여럿이면 fiducial
    마다 포화되지 않은 것 중 가장 밝은 세트를 고르고, 전부 포화면 포화 화소가
    가장 적은 세트를 쓴다. fid_min_peak 보다 어두운 검출은 (다른 세트가 있으면)
    피한다.

    반환: xobs, yobs (fiducial 순서), choice (fiducial 마다 고른 세트 번호),
          theta_guess, (xoff, yoff) (마지막 매칭의 초기값; fitdistortion 용)
    """
    log = get_logger()
    nfid, nfib = int(fid_flag.sum()), int((~fid_flag).sum())
    sets = []
    for k, head in enumerate(fid_heads):
        _, xobs, yobs, peak, nsat = findpeak((nfid, nfid+nfib+n_extra), data_dir=data_dir, head=head,
                                             nexposure=nexposure, sat_level=sat_level,
                                             ReturnNsat=True, **pk)
        if xobs.size > nfid:
            # fiber 도 보인다: 전체 패턴으로 매칭하고 fiducial 자리만 쓴다
            imatch, theta_guess, (xoff, yoff) = matchfiber(x, y, xobs, yobs, m_init=m_init)
            imatch = imatch[fid_flag]
        else:
            imatch, theta_guess, (xoff, yoff) = matchfiber(x[fid_flag], y[fid_flag], xobs, yobs,
                                                           m_init=m_init)
        if np.unique(imatch).size != nfid:
            log.warning("Fiducial set %d (%s): %d of %d fiducials matched to the same peak",
                        k, head, nfid-np.unique(imatch).size, nfid)
        sets.append(dict(x=xobs[imatch], y=yobs[imatch], peak=peak[imatch], nsat=nsat[imatch],
                         theta=theta_guess, off=(xoff, yoff)))
        log.info("Fiducial set %d (%s): %d peaks, %d fiducials saturated, peak %.0f..%.0f",
                 k, head, xobs.size, int((nsat[imatch] > 0).sum()),
                 peak[imatch].min(), peak[imatch].max())

    choice = np.zeros(nfid, dtype=int)
    for i in range(nfid):
        ok = [k for k, st in enumerate(sets) if st['nsat'][i] == 0 and st['peak'][i] >= fid_min_peak]
        if ok:
            choice[i] = max(ok, key=lambda k: sets[k]['peak'][i])
        else:
            choice[i] = min(range(len(sets)), key=lambda k: (sets[k]['nsat'][i], -sets[k]['peak'][i]))
    xobs = np.array([sets[choice[i]]['x'][i] for i in range(nfid)])
    yobs = np.array([sets[choice[i]]['y'][i] for i in range(nfid)])
    bad = [i for i in range(nfid) if sets[choice[i]]['nsat'][i] > 0]
    if len(sets) > 1:
        log.info("Fiducial exposure chosen per fiducial: %s",
                 ", ".join(f"set {k}: {int((choice == k).sum())}" for k in range(len(sets))))
    if bad:
        log.warning("%d fiducial(s) saturated in every set: %s", len(bad),
                    ", ".join(str(i) for i in bad))
    last = sets[-1]
    return xobs, yobs, choice, last['theta'], last['off']


def measure_fibers(head, x, y, nfib, coeff, data_dir, nexposure, max_dist=None, n_extra=0, **pk):
    """fiber 노출에서 peak 을 찾아 (fiducial 로 맞춘) 변환으로 focal 좌표로 옮기고
    목표 위치에 가장 가까운 것을 그 fiber 로 본다. x, y 는 fiber + fiducial 전체
    목표 (fiducial 은 그 자리의 검출을 fiber 로 잘못 잡지 않게 하려고 넣는다).

    반환: xfocal, yfocal (x, y 순서), dist (매칭 거리 [mm]), (xobs, yobs) (같은 순서의 칩 좌표)
    """
    log = get_logger()
    nfid = x.size - nfib
    _, xobs, yobs, _ = findpeak((nfib, nfib+nfid+n_extra), data_dir=data_dir, head=head,
                                nexposure=nexposure, **pk)
    xf, yf = transform(xobs, yobs, coeff)
    dist, j = cKDTree(np.column_stack((xf, yf))).query(np.column_stack((x, y)))
    if np.unique(j[:nfib]).size != nfib:
        log.warning("Fiber matching is not unique: %d fibers share a peak",
                    nfib-np.unique(j[:nfib]).size)
    if max_dist is not None and (dist[:nfib] > max_dist).any():
        log.warning("%d fiber(s) matched farther than %.2f mm from the target (max %.2f mm)",
                    int((dist[:nfib] > max_dist).sum()), max_dist, dist[:nfib].max())
    log.info("Fiber image (%s): %d peaks; fiber-target distance median %.1f um, max %.1f um",
             head, xobs.size, np.median(dist[:nfib])*1e3, dist[:nfib].max()*1e3)
    return xf[j], yf[j], dist, (xobs[j], yobs[j])


def mtlcal(data_dir='./MTL/data/'
           , head='test'
           , mode='Raw'
           , threshold=3e3
           , nwindow=40
           , niter_recenter=0
           , apply_pinhole_offset=False
           , use_measured=None
           , finder='sep'
           , finder_opts=None
           , background='crop'
           , background_opts=None
           , nexposure=1
           , fid_head=None
           , fid_nexposure=None
           , model='full'
           , magnification=None
           , sat_level=65000.
           , fid_min_peak=0.
           , n_extra=0
           , target_file=None
           , json_dir=None
           , target_name=None
           , itrial=1):
    """이미지에서 fiber 위치를 재고 positioner 보정각을 구한다.

    fid_head 가 None 이면 이미지 한 벌(head)에 fiber 와 fiducial 이 같이 찍혀
    있다고 보고 한 번에 처리한다 (기본). fiber 와 fiducial 의 밝기 차이가 커서
    노출을 나눠 찍었으면 fiducial 노출의 head 를 fid_head 에 준다 (하나 또는
    여러 개의 리스트: fiducial 밝기 차이가 심해 두 노출로 나눠 찍은 경우). 그러면

        1. fiducial 노출(들)에서 fiducial 위치를 재고  (measure_fiducials)
        2. fiducial 로 camera -> focal 변환을 맞춘 뒤 (fitdistortion)
        3. fiber 노출(head)의 peak 을 그 변환으로 옮겨 목표에 매칭한다 (measure_fibers)

    두 노출 사이에 카메라와 plate 가 움직이지 않았다고 가정한다.

    model 은 왜곡 다항식에서 풀 항이다 (utils.DISTORTION_MODELS: 'full', 'optics',
    'radial', 'similarity'). 렌즈가 적은 실험실 광학계에서는 'radial' 처럼 줄인다.
    magnification 은 camera -> focal 배율 초기값으로, 본 장비와 광학계가 다르면
    (실험실) 꼭 주어야 fiber 매칭이 된다 (utils.CAMERA2FOCAL_M 참고; 부호 포함).
    sat_level 은 포화 판정 수준 [ADU], fid_min_peak 은 fiducial 노출을 고를 때 이보다
    어두운 검출을 피하는 기준이다. n_extra 는 표에 없는데 이미지에 보이는 밝은 점의 수
    (예: 고정된 sky fiber 9개). 그만큼 검출을 더 허용하고 매칭에서 자연히 빠진다.

    mode는 findpeak에 그대로 넘어간다. "Raw"는 threshold를 조절해 가며 peak을
    찾고, "Predict"는 목표 위치에서 출발해 가장 밝은 픽셀로 정렬한다.

    nwindow는 findpeak이 center of mass를 잴 때 쓰는 crop 반폭 [pixel]이다.
    이 창 안에 이웃 spot이 들어오면 무게중심이 끌려가 위치를 잘못 잰다. 창의
    반폭을 focal plane 거리로 환산하면 nwindow * 3.76e-3 mm/px * 배율(~8.47)
    이므로, fiber 사이 최소 이격 거리의 절반보다 작아야 한다.
    최소 이격 3mm 기준이면 nwindow <= 40 이어야 한다.

    use_measured 는 설계 좌표와 arm 길이를 실측값으로 바꿀지다
    (load_configuration 참고). None이면 fiber_config.USE_MEASURED를 따른다.

    apply_pinhole_offset 은 fiducial 기준 좌표에 핀홀 오프셋을 더할지다
    (load_configuration 참고). 기본값은 False이고, 쓰는 plate에서 측정한
    핀홀 파일을 넣은 뒤에 True로 켠다.

    niter_recenter 는 창을 무게중심으로 옮겨 다시 재는 횟수다 (findpeak 참고).
    PSF가 창에 비해 크면 켜는 것이 좋다.

    finder / finder_opts 는 peak을 찾는 방법과 그 인자다 ('find_peaks',
    'daofind', 'sep', 'segmentation'). background / background_opts 는 배경
    제거 방식이다. 둘 다 findpeak으로 그대로 넘어간다 (peakfind 참고).

    json_dir와 target_name을 주면 누적 각도를 trial별 json으로 저장한다.
    누적 각도는 직전 trial의 json에 이번 회전량을 더한 값이라, itrial-1 파일이
    있어야 한다 (없으면 목표 각도에서 출발한다).

    반환:
        dx, dy    : fiber와 fiducial 전체의 목표 대비 위치 오차 (mm)
        angle_rot : fiber 순서의 이번 trial 회전량 (deg)
        angle_cum : fiber 순서의 누적 각도 (deg). json에 저장하는 값과 같다.
    """
    log = get_logger()

    xorigin, yorigin, x, y, fid_flag, nfib = load_configuration(
        target_file, apply_pinhole_offset=apply_pinhole_offset,
        use_measured=use_measured)

    pk = dict(threshold=threshold, nwindow=nwindow, niter_recenter=niter_recenter,
              finder=finder, finder_opts=finder_opts, background=background,
              background_opts=background_opts)

    # 끊어진 fiber (fiber_config.FIBER_EXCLUDE): 찾지 않고 매칭에서 뺀다. 배열 순서는 그대로.
    dead = np.r_[cfg.dead_mask(), np.zeros(x.size-nfib, dtype=bool)]
    if dead.any():
        log.info("Excluded fiber(s) (not expected in the image): %s",
                 ", ".join(h for h, d in zip(cfg.fiber_ids(), dead[:nfib]) if d))
    xa, ya, fa = x[~dead], y[~dead], fid_flag[~dead]
    nfib_a = int((~dead[:nfib]).sum())

    if fid_head is None:
        #---이미지 한 벌: fiber 와 fiducial 을 같이 찾는다-------------------------
        npeaks = xa.size if (mode != "Raw" or n_extra == 0) else (xa.size, xa.size + n_extra)
        _, xobs, yobs, _, nsat = findpeak(npeaks, data_dir=data_dir, head=head,
                                          nexposure=nexposure, mode=mode, x=xa, y=ya,
                                          sat_level=sat_level, ReturnNsat=True, **pk)

        imatch, theta_guess, (xoff_guess, yoff_guess) = matchfiber(
            xa, ya, xobs, yobs, m_init=magnification)

        nsat_fid = int((nsat[imatch][fa] > 0).sum())
        if nsat_fid:
            log.warning("%d fiducial(s) saturated; consider a separate fiducial "
                        "exposure (fid_head / fid_exptime)", nsat_fid)

        xf_a, yf_a, _, _, _ = fitdistortion(xa, ya, fa
                                            , xobs, yobs
                                            , imatch, theta_guess
                                            , xoff_guess, yoff_guess
                                            , model=model, m_init=magnification)
    else:
        #---노출을 나눠 찍은 경우: fiducial 노출로 변환을 맞추고 fiber 노출에 적용----
        fid_heads = [fid_head] if isinstance(fid_head, str) else list(fid_head)
        fid_nexposure = nexposure if fid_nexposure is None else fid_nexposure
        x_fid, y_fid = x[fid_flag], y[fid_flag]

        xo_fid, yo_fid, _, theta_guess, (xoff_guess, yoff_guess) = measure_fiducials(
            fid_heads, xa, ya, fa, data_dir, fid_nexposure, n_extra=n_extra,
            m_init=magnification, sat_level=sat_level, fid_min_peak=fid_min_peak, **pk)

        xf_fid, yf_fid, dxf, dyf, coeff = fitdistortion(
            x_fid, y_fid, np.ones(x_fid.size, dtype=bool), xo_fid, yo_fid,
            np.arange(x_fid.size), theta_guess, xoff_guess, yoff_guess,
            model=model, m_init=magnification)

        xf_a, yf_a, _, _ = measure_fibers(head, xa, ya, nfib_a, coeff, data_dir, nexposure,
                                          n_extra=n_extra, **pk)
        xf_a[fa], yf_a[fa] = xf_fid, yf_fid                      # fiducial 은 전용 노출 값

    xfocal = np.full(x.size, np.nan); yfocal = np.full(y.size, np.nan)
    xfocal[~dead], yfocal[~dead] = xf_a, yf_a                   # 끊어진 fiber 는 NaN
    dx, dy = xfocal-x, yfocal-y

    #---보정각 계산 (fiber만)---------------------------------------------------
    theta_true, phi_true = _arm_angles(xorigin, yorigin, x, y, use_measured)
    theta_obs, phi_obs = _arm_angles(xorigin, yorigin, xfocal, yfocal,
                                     use_measured)

    dtheta, dphi = _wrap(theta_true - theta_obs, phi_true - phi_obs)
    dtheta[dead[:nfib]] = 0.; dphi[dead[:nfib]] = 0.             # 끊어진 fiber 는 안 움직인다

    angle_rot = _flatten_angles(dtheta, dphi)

    #---누적 각도---------------------------------------------------------------
    prev = None
    if json_dir is not None and target_name is not None:
        prev_path = json_path(json_dir, target_name, itrial-1)
        if os.path.exists(prev_path):
            with open(prev_path, 'r') as ff:
                prev_data = json.load(ff)
            prev = angle_array(prev_data)
        else:
            log.warning("Previous trial file not found: %s", prev_path)

    if prev is None:
        prev = _flatten_angles(theta_true, phi_true)

    angle_cum = prev + angle_rot

    if json_dir is not None and target_name is not None:
        log.info("Saved angles to %s",
                 save_angles(json_dir, target_name, itrial, angle_cum))
        log.info("Saved fiber offsets to %s",
                 save_offsets(json_dir, target_name, itrial, dx, dy, angle_rot, nfib))

    dist = np.sqrt(dx[~fid_flag & ~dead]**2 + dy[~fid_flag & ~dead]**2)*1e3
    log.info("Fiber position error: max %.1f um, median %.1f um",
             dist.max(), np.median(dist))

    return dx, dy, angle_rot, angle_cum
