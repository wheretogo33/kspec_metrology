def hot_pixel_removal_median_ratio(
    img: np.ndarray,
    factor: float = 5.0,            # 주변 median의 몇 배 이상이면 제거할지
    n_iter: int = 1,
    mode: str = "mirror",
    saturated_value: int | float | None = None,
    eps: float = 1e-6,              # median이 0 근처일 때 0나눔/과잉검출 방지
    abs_threshold: float | None = None,  # (선택) |P - median|이 이것보다 커야 제거
    keep_dtype: bool = True,
):
    """
    주변 8픽셀 median을 기준으로:
      P > median_neighbors * factor  이면 핫픽셀로 보고,
    해당 픽셀을 median_neighbors로 치환.

    abs_threshold를 같이 쓰면 (P - median) 절대 차이도 커야 제거되므로
    어두운 배경에서 과잉 검출되는 것을 줄일 수 있습니다.
    """
    img = np.asarray(img)
    work = img.astype(np.float32, copy=True)

    # 중앙 제외한 8-neighbor footprint
    footprint = np.array([[1, 1, 1],
                          [1, 0, 1],
                          [1, 1, 1]], dtype=np.uint8)

    for _ in range(max(1, int(n_iter))):
        med_nb = median_filter(work, footprint=footprint, mode=mode)

        # ratio 조건: P > median * factor
        denom = np.maximum(np.abs(med_nb), eps)   # median이 0일 때 폭주 방지
        mask = work > denom * float(factor)

        # (선택) 절대 차이 조건도 추가
        if abs_threshold is not None:
            mask &= (work - med_nb) > float(abs_threshold)

        # saturation 값은 제외하고 싶다면
        if saturated_value is not None:
            mask &= (work != float(saturated_value))

        # 치환
        work[mask] = med_nb[mask]

    # dtype 복구
    if keep_dtype:
        if np.issubdtype(img.dtype, np.integer):
            info = np.iinfo(img.dtype)
            work = np.clip(np.rint(work), info.min, info.max).astype(img.dtype)
        else:
            work = work.astype(img.dtype)

    return work

cleaned = hot_pixel_removal_median_ratio(im, factor=1.5, n_iter=2)

from astropy.io import fits

hdr = fits.Header()
hdr['Gain'] = 0
hdr['offset'] = 0
hdr['texp'] = 0

empty_primary = fits.PrimaryHDU(header=hdr, data=cleaned)
empty_primary.writeto(f'hotpix_test_1.5median.fits', overwrite=True)
