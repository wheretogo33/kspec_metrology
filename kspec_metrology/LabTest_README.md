# kspec_metrology_LabTest — 실험실 시험용 묶음

인터넷이 없는 실험실 컴퓨터용. 설치하지 않고 폴더째 쓴다. 이미 깔린 옛
`kspec_metrology` 와 섞이지 않도록 패키지 이름이 `kspec_metrology_LabTest` 다
(import 문도 그 이름).

```
kspec_metrology_LabTest/     패키지 (analysis, exposure, lab, mtlrun ...)
LabTest_Metrology.ipynb      순서대로 실행하는 노트북
LabTest_README.md            이 파일
```

## 쓰는 법

1. zip 을 아무 폴더에나 푼다 (예: `/home/lab/metrology/`).
2. `LabTest_Metrology.ipynb` 를 열고 첫 셀의 `PKG_DIR` 을 그 폴더로, `WORK_DIR` 을
   이미지·결과를 쓸 폴더로, `POSITIONERS` 를 꽂힌 홀 ID 로 바꾼다.
3. 셀을 순서대로 실행한다. 설정 셀을 바꾸면 커널을 재시작하고 처음부터 다시.

필요한 파이썬 패키지: numpy, scipy, astropy, photutils, sep, matplotlib, jupyter
(옛 kspec_metrology 가 돌던 환경이면 있다). 카메라는 `exposure/lib/libqhyccd.so`
를 쓰며 패키지에 들어 있다.

## 목표 위치 파일

한 줄에 `ID  v1  v2` (공백/쉼표 구분, `#` 뒤 주석). 표의 positioner 전부가 있어야 한다.

```
# ID   x[mm]    y[mm]
C4    -1.4812   82.3105
C5    17.5000   82.2724
```

- `KIND = 'xy'`: 초점면 목표 (x, y) [mm]
- `KIND = 'ab'`: positioner 각도 (alpha, beta) [deg]. json 과 같은 부호 규약이며,
  회전 중심(실측값)과 arm 길이로 (x, y) 로 바꾼다.

## 이미지 이름 (이미 찍어 둔 이미지로 계산할 때)

`WORK_DIR/data/` 에 다음 이름으로 둔다 (`i` = trial 번호, 1 부터).

```
{TILE}_MetrologyTrial_{i}_0.fits          fiber 노출
{TILE}_MetrologyTrial_{i}_fid0_0.fits     fiducial 전용 노출 세트 0
{TILE}_MetrologyTrial_{i}_fid1_0.fits     fiducial 전용 노출 세트 1 (FID_EXPTIME 이 리스트일 때)
```

## 결과

`WORK_DIR/json/` 에 trial 마다

- `{TILE}_MetrologyTrial_{i}.json` — 누적 각도 (positioner 로 보낼 값)
- `{TILE}_MetrologyTrial_{i}_offset.txt` — fiber 별 남은 위치 오차 [µm] 와 이번 회전량 [deg]

자세한 설명은 패키지 안 `README.md` (본 장비용 설명이지만 인자는 같다).
