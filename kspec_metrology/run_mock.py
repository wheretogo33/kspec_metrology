"""mock 데이터로 metrology 분석 한 바퀴를 끝까지 돌려 보는 예제 스크립트.

카메라 없이, 미리 만들어 둔 mock 이미지 한 장과 target 파일만 가지고
"이미지 -> peak 검출 -> fiber 매칭 -> 왜곡 보정 -> positioner 보정각 json"
전체 과정을 실행한다. 패키지가 제대로 설치됐는지 확인하는 용도로도 쓸 수 있다.

필요한 파일은 두 개뿐이다.

    <image-dir>/test_MetrologyTrial_1_0.fits    mock 이미지 (약 400MB)
    <target-dir>/object.info                    fiber 목표 위치

둘 다 기본값은 이 파일 옆의 tmp/ 이지만, 이미지는 용량이 커서 저장소에 들어
있지 않고 따로 받게 된다. 그래서 폴더를 따로 지정할 수 있다. 파일 이름은
바꾸지 말 것 (이미지 이름은 naming.py의 규칙에서 나온다).

실행:

    # 이미지만 다른 곳에 받아 둔 경우 (가장 흔한 경우)
    python -m kspec_metrology.run_mock --image-dir /받은/폴더

    # 두 파일을 한 폴더에 모아 둔 경우
    python -m kspec_metrology.run_mock --data-dir /어디든/mock

    # 둘 다 tmp/ 에 있는 경우
    python -m kspec_metrology.run_mock

    # object.info 경로를 파일째로 지정하거나, 결과만 따로 받기
    python -m kspec_metrology.run_mock --target /경로/object.info --out-dir ./result

결과:

    <out-dir>/test_MetrologyTrial_0.json   관측 전 목표 각도
    <out-dir>/test_MetrologyTrial_1.json   이 trial의 누적 보정각  <- 최종 산출물

json의 key는 "{hole ID}_a" (alpha arm), "{hole ID}_b" (beta arm)이고 값은
degree다. 자세한 것은 README의 "출력 형식" 참고.

mock을 만든 쪽에서 정답 파일(test_MetrologyTrial_1_truth.npz)을 같이 줬다면
자동으로 찾아서 결과를 검증한다 (target-dir -> image-dir 순으로 찾는다).
없어도 그냥 넘어간다.
"""

import argparse
import json
import os
import sys

import numpy as np

from kspec_metrology import naming
from kspec_metrology.analysis import fiber_config as cfg
from kspec_metrology.analysis.mtlcal import (angle_array, mtlcal,
                                             save_initial_angles)

# 이 파일 옆의 tmp/
DEFAULT_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tmp')

TILE = 'test'
ITRIAL = 1


def check_inputs(image_dir, target, tile, itrial, nexposure):
    """필요한 파일이 다 있는지 먼저 확인하고, 없으면 어디를 봤는지 알려준다."""
    images = naming.image_paths(image_dir, tile, itrial, nexposure)

    missing = [p for p in [target] + images if not os.path.exists(p)]
    if missing:
        print("[!] 다음 파일을 찾을 수 없습니다:", file=sys.stderr)
        for p in missing:
            print(f"      {p}", file=sys.stderr)
        print("\n    이미지 폴더는 --image-dir, object.info 폴더는 --target-dir 로",
              file=sys.stderr)
        print("    지정합니다 (한 폴더에 다 있으면 --data-dir 하나로 충분).",
              file=sys.stderr)
        print("    파일 이름은 바꾸지 말아야 합니다.", file=sys.stderr)
        raise SystemExit(1)

    return images


def find_truth(dirs, tile, itrial):
    """정답 파일을 주어진 폴더들에서 순서대로 찾는다. 없으면 None."""
    name = f'{tile}_MetrologyTrial_{itrial}_truth.npz'
    for d in dirs:
        path = os.path.join(d, name)
        if os.path.exists(path):
            return path
    return None


def check_against_truth(truth_path, tile, itrial, out_dir, target_file):
    """정답 파일이 있으면 결과를 검증한다. 없으면 조용히 넘어간다."""
    if truth_path is None:
        return

    from kspec_metrology.analysis.mtlcal import load_configuration

    xorigin, yorigin, _, _, fid_flag, nfib = load_configuration(target_file)
    arm1, arm2 = cfg.arm_lengths()

    tgt = json.load(open(target_file))
    xp, yp = np.asarray(tgt['xp'], float), np.asarray(tgt['yp'], float)

    def forward(angles):
        """json 각도 -> fiber 끝 위치. 출력각 = -계산각 규약을 되돌린다."""
        th = -np.deg2rad(angles[0::2])
        ph = -np.deg2rad(angles[1::2])
        return (xorigin[:nfib] + arm1*np.cos(th) + arm2*np.cos(th + ph),
                yorigin[:nfib] + arm1*np.sin(th) + arm2*np.sin(th + ph))

    a0 = angle_array(json.load(open(naming.json_path(out_dir, tile, 0))))
    a1 = angle_array(json.load(open(naming.json_path(out_dir, tile, itrial))))

    # Trial 0 각도를 순기구학으로 풀면 target의 xp, yp가 그대로 나와야 한다.
    x0, y0 = forward(a0)
    e0 = np.hypot(x0 - xp, y0 - yp)*1e3

    # 구동계에 고정 offset이 있다고 보는 모델에서, 보정각을 명령하면 목표에 닿는다.
    delta = a0 - a1                      # = -(이번 회전량)
    xa, ya = forward(a1 + delta)
    e1 = np.hypot(xa - xp, ya - yp)*1e3
    xb, yb = forward(a0 + delta)         # 보정 전에 가 있던 자리
    eb = np.hypot(xb - xp, yb - yp)*1e3

    print("--- 정답 파일로 검증 -------------------------------------------")
    print(f"  Trial 0 각도 -> 목표 위치 재현      : max {e0.max():.2e} um")
    print(f"  보정 전 위치 오차                   : median {np.median(eb):7.2f} um")
    print(f"  보정각 적용 후 목표까지             : max {e1.max():.2e} um")
    print()


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='mock 데이터로 metrology 분석 전체를 돌린다.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--data-dir', default=None,
                    help='이미지와 object.info가 한 폴더에 다 있을 때 그 폴더. '
                         '--image-dir와 --target-dir의 기본값이 된다')
    ap.add_argument('--image-dir', default=None,
                    help='mock 이미지(.fits)가 있는 폴더 '
                         f'(기본: --data-dir, 없으면 {DEFAULT_DATA_DIR})')
    ap.add_argument('--target-dir', default=None,
                    help='object.info가 있는 폴더 '
                         f'(기본: --data-dir, 없으면 {DEFAULT_DATA_DIR})')
    ap.add_argument('--target', default=None,
                    help='object.info 경로를 파일째로 지정 (--target-dir보다 우선)')
    ap.add_argument('--out-dir', default=None,
                    help='결과 json을 쓸 폴더 (기본: object.info가 있는 폴더)')
    ap.add_argument('--tile', default=TILE, help='타일 이름 (파일 이름에 쓰인다)')
    ap.add_argument('--itrial', type=int, default=ITRIAL, help='trial 번호')
    ap.add_argument('--nexposure', type=int, default=1, help='이미지 장수')
    ap.add_argument('--mode', default='Raw', choices=['Raw', 'Predict'],
                    help='peak 검출 방식')
    ap.add_argument('--threshold', type=float, default=3e3,
                    help='Raw mode의 검출 문턱값')
    ap.add_argument('--nwindow', type=int, default=40,
                    help='center of mass crop 반폭 [pixel]. fiber 최소 이격 '
                         '3mm 기준이 40이다')
    args = ap.parse_args(argv)

    # --data-dir 가 두 폴더의 기본값이고, 각각 따로 덮어쓸 수 있다.
    base = args.data_dir or DEFAULT_DATA_DIR
    image_dir = os.path.abspath(args.image_dir or base)
    target_dir = os.path.abspath(args.target_dir or base)

    target_file = (os.path.abspath(args.target) if args.target
                   else os.path.join(target_dir, 'object.info'))

    out_dir = (os.path.abspath(args.out_dir) if args.out_dir
               else os.path.dirname(target_file))

    images = check_inputs(image_dir, target_file, args.tile, args.itrial,
                          args.nexposure)
    truth_path = find_truth([os.path.dirname(target_file), image_dir],
                            args.tile, args.itrial)

    print("=== 입력 ========================================================")
    print(f"  target : {target_file}")
    for p in images:
        print(f"  image  : {p}")
    if truth_path:
        print(f"  truth  : {truth_path}")
    print(f"  fiber  : {len(cfg.fiber_ids())} positioners "
          f"(arm {cfg.DEFAULT_ARM[0]}, {cfg.DEFAULT_ARM[1]} mm)")
    print(f"  출력   : {out_dir}")
    print()

    os.makedirs(out_dir, exist_ok=True)

    # 1) 관측 전 목표 각도 (Trial 0).
    #    object.info의 xp,yp + 홀 설계 좌표 + arm 길이로 역기구학을 풀어 얻는다.
    #    itrial >= 2 에서는 직전 trial json이 반드시 있어야 하므로,
    #    trial을 이어서 돌 때는 이 파일을 지우지 말 것.
    save_initial_angles(out_dir, args.tile, target_file)

    # 2) 이미지에서 위치를 재고 보정각을 계산해 Trial N json으로 저장
    dx, dy, angle_rot, angle_cum = mtlcal(
        data_dir=image_dir,
        head=naming.image_head(args.tile, args.itrial),
        mode=args.mode,
        threshold=args.threshold,
        nwindow=args.nwindow,
        nexposure=args.nexposure,
        target_file=target_file,
        json_dir=out_dir,
        target_name=args.tile,
        itrial=args.itrial,
    )

    nfib = len(cfg.fiber_ids())
    err = np.hypot(dx[:nfib], dy[:nfib])*1e3

    print()
    print("=== 결과 ========================================================")
    print(f"  fiber 위치 오차 : median {np.median(err):.2f} um, max {err.max():.2f} um")
    print(f"  회전량          : median {np.median(np.abs(angle_rot)):.4f} deg, "
          f"max {np.abs(angle_rot).max():.4f} deg")
    print()

    out_json = naming.json_path(out_dir, args.tile, args.itrial)
    data = json.load(open(out_json))
    print(f"  최종 json : {out_json}")
    print(f"              key {len(data)}개 (= fiber {nfib} x arm 2)")
    for k in list(data)[:4]:
        print(f"                {k:>8s} : {data[k]:12.4f} deg")
    print("                     ...")
    print()

    check_against_truth(truth_path, args.tile, args.itrial, out_dir, target_file)

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
