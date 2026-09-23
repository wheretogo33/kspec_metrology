"""사용하는 fiber / fiducial 목록과 positioner arm 길이 정의.

fiber 구성이 바뀌면 이 파일만 고치면 된다. fiber 개수와 arm 길이가 모두 여기서
나오므로 다른 코드는 손댈 필요가 없다.

지금은 Fiber_Configuration 표의 FiducialFlag == 0 인 홀을 "표에 적힌 순서대로"
모두 positioner로 쓴다 (load_fibers 참고). 어느 홀에 positioner를 꽂았는지가
바뀌면 코드가 아니라 표의 FiducialFlag를 고치면 된다 (0 = positioner,
1 = fiducial, -9 = 미장착). 설정별로 표를 여러 개 두고 쓰려면
KSPEC_FIBER_TABLE 환경변수를 보라.

표의 순서가 곧

    - object.info의 xp, yp 순서
    - mtlcal.angle_dict가 만드는 출력 json의 key 순서

라서, 표를 건드리면 target 파일도 같은 순서로 다시 만들어야 한다.

FIBERS의 각 행:
    ID      Fiber_Configuration 표의 fiber ID (= Hole ID)
    arm1    positioner 안쪽 arm 길이 (mm)
    arm2    positioner 바깥쪽 arm 길이 (mm)
"""

import os
from pathlib import Path

import numpy as np
from astropy.io import ascii

# fiber/fiducial의 설계 좌표가 들어 있는 표.
#
# 기본값은 패키지 안에 같이 들어 있는 파일이다. 설정(어느 홀에 positioner를
# 꽂았는지)이 여럿이면 표를 패키지 밖에 이름 붙여 따로 두고 KSPEC_FIBER_TABLE
# 환경변수로 고르는 편이 낫다. 패키지 안의 파일을 덮어쓰면 pip install 때
# 날아가고, 지금 어느 설정으로 돌고 있는지도 안 보인다.
#
#   KSPEC_FIBER_TABLE=~/configs/Fiber_Configuration_phase1.txt python ...
#
# 표를 바꾸면 positioner 개수와 순서가 바뀌므로 target 파일(object.info)의
# xp, yp도 같은 순서로 다시 만들어야 한다. xp가 더 많으면 앞에서 잘라 쓰면서
# 조용히 지나가므로(mtlcal.load_configuration) 특히 주의할 것.
FIBER_TABLE_ENV = 'KSPEC_FIBER_TABLE'
DEFAULT_FIBER_TABLE_PATH = str(Path(__file__).with_name('Fiber_Configuration_250415.txt'))

FIBER_TABLE_PATH = os.environ.get(FIBER_TABLE_ENV) or DEFAULT_FIBER_TABLE_PATH

if not os.path.exists(FIBER_TABLE_PATH):
    raise FileNotFoundError(
        f"Fiber configuration table not found: {FIBER_TABLE_PATH}"
        + (f" (from ${FIBER_TABLE_ENV})" if os.environ.get(FIBER_TABLE_ENV)
           else ""))

# 관측할 타일/타겟 정보 (fiber별 목표 위치 xp, yp)
TARGET_INFO_PATH = '/home/kspecmtl/work/KSPEC_ICS/MTL/target/object.info'

# Fiducial 핀홀의 실제 위치.
#
# fiducial은 홀 한가운데가 아니라 조금 치우쳐 박혀 있고, 그 치우친 양(dx, dy)을
# 따로 재서 npz로 들고 있다. 그 npz의 index가 어느 홀인지는 map 파일에 적는다
# (현장에서 핀홀을 옮기면 map 파일의 ID만 고치면 된다).
#
# 두 파일이 다 있으면 load_configuration이 fiducial 기준 좌표를
# '표의 홀 좌표 + dx, dy' 로 쓴다. 없으면 표의 홀 좌표를 그대로 쓴다.
FIDUCIAL_PINHOLE_LOC_PATH = os.environ.get('KSPEC_FIDUCIAL_PINHOLE_LOC') or str(
    Path(__file__).with_name('Fiducial_pinhole_loc.npz'))
FIDUCIAL_PINHOLE_MAP_PATH = os.environ.get('KSPEC_FIDUCIAL_PINHOLE_MAP') or str(
    Path(__file__).with_name('Fiducial_pinhole_map.txt'))


def pinhole_offsets(map_path=None, loc_path=None):
    """Hole ID -> (dx, dy) [mm]. 파일이 없으면 빈 dict을 돌려준다.

    map 파일은 "index  ID" 두 칸이고 '#' 뒤는 주석이다. ID가 '-' 인 줄은
    쓰지 않는 index로 보고 건너뛴다.
    """
    map_path = map_path or FIDUCIAL_PINHOLE_MAP_PATH
    loc_path = loc_path or FIDUCIAL_PINHOLE_LOC_PATH

    missing = [q for q in (map_path, loc_path) if not os.path.exists(q)]
    if missing:
        return {}

    with np.load(loc_path) as d:
        off = {int(i): (float(a), float(b))
               for i, a, b in zip(d['index'], d['dx'], d['dy'])}

    out = {}
    with open(map_path) as ff:
        for nline, line in enumerate(ff, 1):
            line = line.split('#')[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(
                    f"{map_path}:{nline}: 'index ID' 두 칸이어야 한다: {line!r}")
            idx, hole = parts
            if hole == '-':
                continue
            try:
                key = int(idx)
            except ValueError:
                raise ValueError(
                    f"{map_path}:{nline}: index가 정수가 아니다: {idx!r}") from None
            if key not in off:
                raise ValueError(
                    f"{map_path}:{nline}: index {key} 가 "
                    f"{os.path.basename(loc_path)} 에 없다")
            if hole in out:
                raise ValueError(f"{map_path}:{nline}: {hole} 가 중복이다")
            out[hole] = off[key]

    return out


# arm 길이 측정값이 없는 positioner에 쓰는 설계값 (arm1, arm2) [mm]
DEFAULT_ARM = (5.2, 11.6)

# 측정된 arm 길이가 있는 positioner만 여기에 적는다. 나머지는 DEFAULT_ARM을 쓴다.
#
# 2025-08 시점에 14개 fiber만 측정값이 있었고, 그 값은 아래와 같았다. 지금은 전체
# positioner를 쓰는데 나머지 측정값이 없어서 전부 설계값으로 통일한다. 측정이
# 끝나면 아래 주석을 되살려 넣으면 된다.
#
#   B5  5.13 11.35 | C2  5.29 11.36 | D0  5.21 11.31 | E5  5.23 11.39
#   E8  5.23 11.31 | G4  5.28 11.34 | G11 5.21 11.49 | H7  5.23 11.43
#   I2  5.29 11.34 | I9  5.20 11.32 | K3  5.28 11.44 | K6  5.04 11.18
#   L10 4.98 11.42 | M7  4.98 11.49
ARM_OVERRIDES = {}


def load_fibers(path=None, arm_overrides=None):
    """Fiber_Configuration 표에서 positioner 행을 읽어 FIBERS 형식으로 돌려준다.

    FiducialFlag == 0 인 행이 positioner다. 표에 적힌 순서를 그대로 지킨다.
    """
    tab = ascii.read(path or FIBER_TABLE_PATH)
    tab.rename_columns(tab.colnames[:4], ["ID", "X", "Y", "FiducialFlag"])

    overrides = ARM_OVERRIDES if arm_overrides is None else arm_overrides

    rows = []
    for row in tab:
        if int(row["FiducialFlag"]) != 0:
            continue
        hole = str(row["ID"])
        arm1, arm2 = overrides.get(hole, DEFAULT_ARM)
        rows.append((hole, float(arm1), float(arm2)))

    return tuple(rows)


FIBERS = load_fibers()

# Fiducial은 Fiber_Configuration 표의 FiducialFlag == 1 에서 가져오되,
# 아래 ID는 사용하지 않는다. 지금은 전부 쓴다.
FIDUCIAL_EXCLUDE = ()

# find_angle_double_method2_select_elbow에 넘기는 자세 정의
ELBOW = "down"
CLOCKWISE = False


def fiber_ids():
    """FIBERS 순서 그대로의 fiber ID(= Hole ID) 리스트.

    출력 json의 key가 이 값에서 나온다 (mtlcal.angle_dict 참고).
    """
    return [row[0] for row in FIBERS]


def arm_lengths():
    """FIBERS 순서 그대로의 (arm1, arm2) 배열."""
    return (np.array([row[1] for row in FIBERS], dtype=float),
            np.array([row[2] for row in FIBERS], dtype=float))
