# tmp — mock metrology 데이터

카메라 없이 분석 전체를 돌려 볼 수 있는 mock 한 세트다.
K-SPEC 초점면을 그대로 본뜬 것으로, positioner 151개와 fiducial 30개가 들어 있다.

| 파일 | 크기 | 저장소 | 설명 |
|---|---|---|---|
| `object.info` | 6 KB | 포함 | fiber 목표 위치 (`xp`, `yp`) 151개. `tile_id="test"`, `zenith_angle=0` |
| `test_MetrologyTrial_1_truth.npz` | 20 KB | 포함 | 정답. 검증용이라 분석에는 없어도 된다 |
| `test_MetrologyTrial_1_0.fits` | 416 MB | **제외** | mock metrology 이미지 (8842 x 11760, int32). 용량 때문에 따로 공유한다 |

fits는 `.gitignore`로 빠져 있으니 따로 받아서 아무 폴더에나 두면 된다.
파일 이름은 바꾸지 말 것 (`naming.py`의 규칙에서 나온다).

## 돌려보기

```bash
# 이미지를 따로 받아 둔 경우 (가장 흔한 경우)
python -m kspec_metrology.run_mock --image-dir /이미지를/받은/폴더

# 두 파일을 한 폴더에 모아 둔 경우
python -m kspec_metrology.run_mock --data-dir /어디든/mock

# 이미지까지 이 tmp/ 에 넣어 둔 경우
python -m kspec_metrology.run_mock
```

`--target` 으로 `object.info` 경로를 파일째 지정할 수도 있고, `--out-dir` 로
결과를 다른 곳에 쓸 수도 있다. 전체 옵션은 `--help` 참고.

## 나오는 것

각도 json 두 개가 생긴다 (기본 위치는 `object.info`가 있는 폴더).

```
test_MetrologyTrial_0.json   관측 전 목표 각도
test_MetrologyTrial_1.json   보정각 (최종 산출물)
```

key는 `{hole ID}_a` (alpha arm), `{hole ID}_b` (beta arm), 값은 degree다.

## 기대값

이 mock은 positioner 위치 오차를 sigma = 40 um로 심어 두었다. 제대로 돌면

```
peak 검출        181 / 181
fiber 위치 오차  median 약 53 um, max 약 122 um
```

가 나온다. 정답 파일이 함께 있으면 `run_mock`이 보정각을 순기구학으로 되풀어
목표 위치를 재현하는지까지 확인해 준다 (오차 1e-11 um 수준이면 정상).

## mock을 다시 만들려면

생성기는 이 저장소가 아니라 image_simulation 쪽에 있다.

```
image_simulation/src/Total/Mock_Image_Generator_kspec_metrology.py
```

`fiber_config`와 `Fiber_Configuration_250415.txt`를 그대로 읽어 만들기 때문에,
fiber 구성을 바꾸면 mock도 같이 다시 만들어야 한다. 용량이 부담되면 생성기의
`COMPRESS = True`로 두면 RICE 무손실 압축이 걸려 92MB가 되고, 확장자가 `.fits`
그대로라 읽는 쪽 코드는 바꿀 필요가 없다.
