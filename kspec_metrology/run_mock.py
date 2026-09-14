"""mock 데이터로 metrology 분석을 한 번 돌려 보는 예제.

하는 일은 mtlcal 한 번 호출하는 것뿐이다. 카메라 없이 미리 만들어 둔 이미지와
target 파일만 가지고 "이미지 -> peak 검출 -> fiber 매칭 -> 왜곡 보정 ->
positioner 보정각 json"이 끝까지 도는지 확인하는 용도다. 나머지는 파일 경로를
받아 넘기고 결과를 보기 좋게 찍는 것뿐이다.

필요한 파일은 두 개다.

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

Trial 0 (관측 전 목표 각도) json은 만들지 않는다. 실제 운용에서는
MetrologyRun.start() 가 관측 전에 그것부터 저장하지만, 여기서는 이미 찍혀 있는
이미지를 분석만 하므로 필요 없다. 직전 trial json이 없다고 경고가 한 줄 나오는
것은 정상이고, mtlcal이 목표 각도를 그 자리에서 계산해 쓴다.
"""

import argparse
import json
import os
import sys

import numpy as np

from kspec_metrology import naming
from kspec_metrology.analysis import fiber_config as cfg
from kspec_metrology.analysis.mtlcal import mtlcal

# 이 파일 옆의 tmp/
DEFAULT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tmp')


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='mock 데이터로 metrology 분석을 돌린다 (mtlcal 한 번).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--data-dir', default=None,
                    help='이미지와 object.info가 한 폴더에 다 있을 때 그 폴더. '
                         '--image-dir와 --target-dir의 기본값이 된다')
    ap.add_argument('--image-dir', default=None,
                    help=f'이미지(.fits) 폴더 (기본: --data-dir, 없으면 {DEFAULT_DIR})')
    ap.add_argument('--target-dir', default=None,
                    help=f'object.info 폴더 (기본: --data-dir, 없으면 {DEFAULT_DIR})')
    ap.add_argument('--target', default=None,
                    help='object.info 경로를 파일째로 지정 (--target-dir보다 우선)')
    ap.add_argument('--out-dir', default=None,
                    help='결과 json을 쓸 폴더 (기본: object.info가 있는 폴더)')
    ap.add_argument('--tile', default='test', help='타일 이름 (파일 이름에 쓰인다)')
    ap.add_argument('--itrial', type=int, default=1, help='trial 번호')
    ap.add_argument('--nexposure', type=int, default=1, help='이미지 장수')
    ap.add_argument('--mode', default='Raw', choices=['Raw', 'Predict'],
                    help='peak 검출 방식')
    ap.add_argument('--threshold', type=float, default=3e3,
                    help='Raw mode의 검출 문턱값')
    ap.add_argument('--nwindow', type=int, default=40,
                    help='center of mass crop 반폭 [pixel]. fiber 최소 이격 '
                         '3mm 기준이 40이다')
    args = ap.parse_args(argv)

    #---경로: --data-dir가 기본값이고 각각 따로 덮어쓸 수 있다------------------
    base = args.data_dir or DEFAULT_DIR
    image_dir = os.path.abspath(args.image_dir or base)
    target_file = (os.path.abspath(args.target) if args.target
                   else os.path.join(os.path.abspath(args.target_dir or base),
                                     'object.info'))
    out_dir = (os.path.abspath(args.out_dir) if args.out_dir
               else os.path.dirname(target_file))

    images = naming.image_paths(image_dir, args.tile, args.itrial, args.nexposure)
    missing = [p for p in [target_file] + images if not os.path.exists(p)]
    if missing:
        print("[!] 다음 파일을 찾을 수 없습니다:", file=sys.stderr)
        for p in missing:
            print(f"      {p}", file=sys.stderr)
        print("\n    이미지 폴더는 --image-dir, object.info 폴더는 --target-dir 로\n"
              "    지정합니다 (한 폴더에 다 있으면 --data-dir 하나로 충분).\n"
              "    파일 이름은 바꾸지 말아야 합니다.", file=sys.stderr)
        return 1

    print(f"target : {target_file}")
    for p in images:
        print(f"image  : {p}")
    print(f"fiber  : {len(cfg.fiber_ids())} positioners "
          f"(arm {cfg.DEFAULT_ARM[0]}, {cfg.DEFAULT_ARM[1]} mm)")
    print()

    os.makedirs(out_dir, exist_ok=True)

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

    #---결과-------------------------------------------------------------------
    # mtlcal이 이미 fiber 위치 오차를 로그로 찍는다. 여기서는 json만 확인한다.
    out_json = naming.json_path(out_dir, args.tile, args.itrial)
    data = json.load(open(out_json))

    print()
    print(f"json  : {out_json}")
    print(f"        key {len(data)}개 (= fiber {len(cfg.fiber_ids())} x arm 2)")
    for k in list(data)[:4]:
        print(f"          {k:>8s} : {data[k]:12.4f} deg")
    print("               ...")
    print(f"회전량 : median {np.median(np.abs(angle_rot)):.4f} deg, "
          f"max {np.abs(angle_rot).max():.4f} deg")

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
