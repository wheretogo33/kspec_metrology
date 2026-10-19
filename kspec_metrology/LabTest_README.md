# kspec_metrology_LabTest — 실험실 시험용 묶음

인터넷이 없는 실험실 컴퓨터용. 설치하지 않고 폴더째 쓴다. 이미 깔린 옛
`kspec_metrology` 와 섞이지 않도록 패키지 이름이 `kspec_metrology_LabTest` 다
(import 문도 그 이름).

```
kspec_metrology_LabTest/     패키지 (analysis, exposure, lab, mtlrun ...)
LabTest_Metrology.ipynb      순서대로 실행하는 노트북 (한 배치를 한 번 측정)
LabTest_README.md            이 파일
```

## 쓰는 법

1. zip 을 아무 폴더에나 푼다 (예: `/home/lab/metrology/`).
2. `LabTest_Metrology.ipynb` 를 열고 첫 셀의 `PKG_DIR`(푼 폴더) 과 `OUT_DIR`(결과 폴더),
   2절의 목표 파일·이미지 경로를 바꾼다.
3. 셀을 순서대로 실행한다. 설정 셀(표)을 바꾸면 커널을 재시작한다.

필요한 파이썬 패키지: numpy, scipy, astropy, photutils, sep, matplotlib, pandas, jupyter
(옛 kspec_metrology 가 돌던 환경이면 있다).

## 입력

**fiber 표**: 패키지 기본 `analysis/Fiber_Configuration_config73.txt`
(positioner 73 + fiducial 30). 바깥 sky fiber 9개 홀은 −9 로 표에서 빠져 있고,
켜져 있으면 `N_EXTRA = 9` 로 "표에 없는 밝은 점" 으로 처리한다.

**목표 파일**: `object.info` 와 같은 json. 배열은 표의 positioner 순서(73개).

```json
{"xp": [...], "yp": [...]}              초점면 목표 [mm]
{"alpha": [...], "beta": [...]}         positioner 각도 [deg]
{"C4_a": ..., "C4_b": ..., ...}         각도 json (이 패키지의 출력 형식)
```

**이미지**: FITS 경로 1–3 개.

| 경우 | `FIBER_IMAGE` | `FID_IMAGES` |
|---|---|---|
| 한 노출에 fiber + fiducial | 그 이미지 | `[]` |
| fiducial 따로 (case 1) | fiber 노출 | `[fiducial 노출]` |
| fiducial 두 노출 (case 2) | fiber 노출 | `[t1_1, t1_2]` (fiducial 마다 포화 안 된 밝은 쪽) |

## 출력 (`OUT_DIR`)

- `{TILE}_MetrologyTrial_1_offset.txt` — fiber 별 남은 위치 오차 dx, dy [µm] 와 회전량 dalpha, dbeta [deg]
- `{TILE}_MetrologyTrial_1.json` — 목표 각도 + 회전량 (positioner 로 보낼 누적 각도)
- `data/` — 분석에 쓴 이미지 링크, `{TILE}_target.info` — (x, y) 로 바꾼 목표

## 파이썬에서 직접

```python
from kspec_metrology_LabTest import lab
res = lab.measure(target_file, fiber_image, [fid1, fid2], out_dir='lab_out',
                  magnification=8.40, n_extra=9)
lab.offsets_table(res)      # fiber 별 표
lab.show_steps(res)         # 단계별 그림
```

자세한 인자 설명은 패키지 안 `README.md`.
