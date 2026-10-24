#!/usr/bin/env python
"""시험 촬영: 원하는 노출시간으로 찍어 원하는 이름으로 저장한다 (fiducial 밝기 조정 등).

    python take_exposure.py --exptime 0.010 --out /data/fid_10ms.fits
    python take_exposure.py --exptime 0.002 0.010 --out /data/{date}_fid_{ms}ms.fits
    python take_exposure.py --exptime 0.4 --n 3 --out /data/fib_400ms.fits        # 3장: ..._0 _1 _2

- exptime 은 MetrologyRun 의 exptime 과 같은 값·단위로 카메라(mtlexp)에 그대로 넘긴다.
  여러 개를 주면 차례로 찍는다.
- 저장 이름: --out 에 {ms} (노출시간 ms), {date} (YYMMDD), {time} (HHMMSS), {i} (장 번호) 를 쓸 수 있다.
  분석 함수(findpeak, lab.check_peaks, lab.fiducial_brightness)가 '...{i}.fits' 꼴을 기대하므로
  {i} 가 없으면 '.fits' 앞에 '_{i}' 를 붙인다 (fid_10ms.fits -> fid_10ms_0.fits).
- 찍은 뒤 최대 화소값과 포화 화소 수를 바로 찍어 준다 (--sat-level, 기본 65000).
- 카메라 설정은 MetrologyRun 기본값과 같다 (--gain 10 --offset 30 --readmode 1 --usb-traffic 40).
- --dry-run 은 카메라 없이 저장될 이름만 보여 준다.
"""
import argparse
import datetime
import os
import sys


def output_paths(pattern, exptime, nexposure, now=None):
    """--out 패턴에서 실제 저장 경로(head, 각 장의 경로)를 만든다."""
    now = now or datetime.datetime.now()
    name = pattern.format(ms=f'{exptime*1e3:g}', date=now.strftime('%y%m%d'), time=now.strftime('%H%M%S'), i='{i}')
    if '{i}' not in name:
        root, ext = os.path.splitext(name)
        name = f'{root}_{{i}}{ext or ".fits"}'
    head = name.replace('{i}', '')
    paths = [name.replace('{i}', str(i)) for i in range(nexposure)]
    return os.path.dirname(os.path.abspath(head)), os.path.basename(head), paths


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__.split('\n', 1)[1])
    ap.add_argument('--exptime', type=float, nargs='+', required=True, help='노출시간 (MetrologyRun exptime 과 같은 단위). 여러 개 가능')
    ap.add_argument('--out', required=True, help="저장 이름. {ms} {date} {time} {i} 사용 가능. 예: /data/{date}_fid_{ms}ms.fits")
    ap.add_argument('--n', type=int, default=1, help='노출시간마다 찍을 장 수 (기본 1)')
    ap.add_argument('--gain', type=float, default=10); ap.add_argument('--offset', type=float, default=30)
    ap.add_argument('--readmode', type=int, default=1); ap.add_argument('--usb-traffic', type=int, default=40)
    ap.add_argument('--sat-level', type=float, default=65000.)
    ap.add_argument('--note', default='', help='FITS header NOTE 에 남길 메모 (예: "G1 up")')
    ap.add_argument('--dry-run', action='store_true', help='찍지 않고 저장될 이름만')
    a = ap.parse_args(argv)

    plan = [(t, *output_paths(a.out, t, a.n)) for t in a.exptime]
    for t, d, head, paths in plan:
        print(f'exptime {t:g}: ' + ', '.join(paths))
    if a.dry_run:
        return 0

    import numpy as np
    from astropy.io import fits
    from kspec_metrology.exposure.mtlexp import mtlexp        # 카메라 SDK 는 여기서야 로드된다
    for t, d, head, paths in plan:
        hdr = {'NOTE': a.note} if a.note else {}
        saved = mtlexp(exptime=t, readmode=a.readmode, usb_traffic=a.usb_traffic, gain=a.gain, offset=a.offset,
                       nexposure=a.n, data_dir=d, head=head, extra_header=hdr)
        for p in saved:
            im = fits.getdata(p)
            nsat = int((im >= a.sat_level).sum())
            print(f'{p}: max {im.max():.0f}, saturated pixels {nsat}' + ('  <-- SATURATED' if nsat else ''))
    return 0


if __name__ == '__main__':
    sys.exit(main())
