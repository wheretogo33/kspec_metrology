import numpy as np
from photutils.detection import find_peaks
from astropy.io import fits
from kspec_metrology.analysis.utils import com, dedupe_peaks_kdtree
from kspec_metrology.logging.log import get_logger
from astropy.table import Table
from scipy.spatial import cKDTree

def findpeak(npeaks
            , data_dir='./MTL/data/'
            , head='test'
            , nexposure=1
            , threshold=5e3
            , boxsize=40
            , nwindow=40
            , SaveFiberImage=False
            , mode="Raw"
            , x=None, y=None):
    log = get_logger()

    xchip = np.linspace(-5879.5, 5879.5, 11760)*3.76e-3
    ychip = np.linspace(-4420.5, 4420.5, 8842)*3.76e-3
    
    #---Stack "nframe" Image----------------------------------------------------------------------------------------------------------
    im = np.zeros((8842, 11760))
    for iframe in range(nexposure):
        im += fits.getdata(data_dir + f'{head}{iframe}.fits').astype(np.float64)[::-1,:] / nexposure


    #---Find Peaks Using any method---------------------------------------------------------------------------------------------------    
    log.info("Start finding peaks with %s", mode)
    if mode == "Raw":
        peak_table_raw = find_peaks(im, threshold=threshold, box_size=boxsize)
        peak_table = dedupe_peaks_kdtree(peak_table_raw, min_dist=40)
        xf, yf = peak_table['x_peak'].data, peak_table['y_peak'].data   
        log.info(f"Found {xf.size} peaks")
        if xf.size > npeaks:
            log.warning("Too many peaks detected")
        elif xf.size < npeaks:
            log.warning(f"{npeaks-xf.size} peaks are not detected")
    elif mode=="Predict":
        from kspec_metrology.analysis.utils import transform, focal2camera_coeff
        coeff_temp = np.copy(focal2camera_coeff)
        coeff_temp[1] = 5.21
        xpredict, ypredict = transform(x, y, coeff_temp)

    #---Calculate Center--------------------------------------------------------------------------------------------------------------
    xobs = np.zeros(npeaks)
    yobs = np.copy(xobs)

    if SaveFiberImage:
        im_crop_full = np.zeros( (npeaks, nwindow*2, nwindow*2))

    for ifiber in range(xf.size):
        if mode != "Predict":
            ipredict = np.argmin( np.abs(xchip[xf[ifiber]]-xchip) )
            jpredict = np.argmin( np.abs(ychip[yf[ifiber]]-ychip) )
        else:
            ipredict = np.argmin( np.abs(xpredict[ifiber]-xchip) )
            jpredict = np.argmin( np.abs(ypredict[ifiber]-ychip) )

        im_crop = im[jpredict-nwindow:jpredict+nwindow, ipredict-nwindow:ipredict+nwindow]

        x_crop = xchip[ipredict-nwindow:ipredict+nwindow]
        y_crop = ychip[jpredict-nwindow:jpredict+nwindow]

        xobs[ifiber], yobs[ifiber] = com(im_crop, x_crop, y_crop)

        if SaveFiberImage:
            im_crop_full[ifiber] = im_crop

    if SaveFiberImage:
        np.save(data_dir+f'fiberimage.npy', im_crop_full)

    return im, xobs, yobs




def merge_peaks_kdtree(t1, t2, tol=1.0,
                       xcol='x_peak', ycol='y_peak', vcol='peak_value',
                       sat_value=65534):

    xy1 = np.vstack([t1[xcol], t1[ycol]]).T.astype(float)
    xy2 = np.vstack([t2[xcol], t2[ycol]]).T.astype(float)

    tree2 = cKDTree(xy2)
    dist12, j = tree2.query(xy1, k=1)

    tree1 = cKDTree(xy1)
    dist21, i_back = tree1.query(xy2, k=1)

    # mutual nearest neighbor + 거리 제한
    mutual = (dist12 <= tol) & (i_back[j] == np.arange(len(xy1)))

    i1 = np.where(mutual)[0]
    i2 = j[mutual]

    used2 = np.zeros(len(t2), dtype=bool)
    used2[i2] = True

    # --- 여기서부터 리스트에 쌓고 마지막에 Table 만들면 dtype 문제 100% 해결 ---
    xs, ys, peaks = [], [], []
    chosens = []
    peak_t1s, peak_t2s = [], []
    sat_t1s, sat_t2s = [], []

    # matched pairs
    for a, b in zip(i1, i2):
        p1 = float(t1[vcol][a])
        p2 = float(t2[vcol][b])
        sat1 = (p1 == sat_value)
        sat2 = (p2 == sat_value)

        if sat1 and (not sat2):
            chosen = 't2'
        elif sat2 and (not sat1):
            chosen = 't1'
        elif (not sat1) and (not sat2):
            chosen = 't1' if (p1 >= p2) else 't2'
        else:
            # 둘 다 sat이면 더 어두운 쪽(더 작은 값) 선택
            chosen = 't1' if (p1 <= p2) else 't2'

        if chosen == 't1':
            x, y, peak = float(t1[xcol][a]), float(t1[ycol][a]), p1
        else:
            x, y, peak = float(t2[xcol][b]), float(t2[ycol][b]), p2

        xs.append(x); ys.append(y); peaks.append(peak)
        chosens.append(chosen)
        peak_t1s.append(p1); peak_t2s.append(p2)
        sat_t1s.append(sat1); sat_t2s.append(sat2)

    # unmatched from t1
    matched1 = np.zeros(len(t1), dtype=bool)
    matched1[i1] = True
    for a in np.where(~matched1)[0]:
        p1 = float(t1[vcol][a])
        xs.append(float(t1[xcol][a]))
        ys.append(float(t1[ycol][a]))
        peaks.append(p1)
        chosens.append('t1')
        peak_t1s.append(p1); peak_t2s.append(np.nan)
        sat_t1s.append(p1 == sat_value); sat_t2s.append(False)

    # unmatched from t2
    for b in np.where(~used2)[0]:
        p2 = float(t2[vcol][b])
        xs.append(float(t2[xcol][b]))
        ys.append(float(t2[ycol][b]))
        peaks.append(p2)
        chosens.append('t2')
        peak_t1s.append(np.nan); peak_t2s.append(p2)
        sat_t1s.append(False); sat_t2s.append(p2 == sat_value)

    out = Table()
    out['x'] = np.array(xs, dtype=float)
    out['y'] = np.array(ys, dtype=float)
    out['peak'] = np.array(peaks, dtype=float)

    out['chosen'] = np.array(chosens, dtype='U2')

    out['peak_t1'] = np.array(peak_t1s, dtype=float)


    out['peak_t2'] = np.array(peak_t2s, dtype=float)
    out['sat_t1']  = np.array(sat_t1s, dtype=bool)
    out['sat_t2']  = np.array(sat_t2s, dtype=bool)

    return out

def Fiducial_findpeak(npeaks
            , data_dir='./MTL/data/'
            , head1='Fiducial_short_Tile_', head2='Fiducial_long_Tile_'
            , nexposure=1
            , threshold=5e3
            , boxsize=40
            , nwindow=40
            , SaveFiberImage=False
            , mode="Raw"
            , x=None, y=None):
    log = get_logger()

    xchip = np.linspace(-5879.5, 5879.5, 11760)*3.76e-3
    ychip = np.linspace(-4420.5, 4420.5, 8842)*3.76e-3
    
    #---Stack "nframe" Image----------------------------------------------------------------------------------------------------------
    im_short = np.zeros((8842, 11760))
    for iframe in range(nexposure):
        im_short += fits.getdata(data_dir + f'{head1}{iframe}.fits').astype(np.float64)[::-1,:] / nexposure

    im_long = np.zeros((8842, 11760))
    for iframe in range(nexposure):
        im_long += fits.getdata(data_dir + f'{head2}{iframe}.fits').astype(np.float64)[::-1,:] / nexposure

    #---Find Peaks Using any method---------------------------------------------------------------------------------------------------    
    peak_table_raw = find_peaks(im_short, threshold=threshold, box_size=boxsize)
    peak_table_short = dedupe_peaks_kdtree(peak_table_raw, min_dist=40)

    peak_table_raw = find_peaks(im_long, threshold=threshold, box_size=boxsize)
    peak_table_long = dedupe_peaks_kdtree(peak_table_raw, min_dist=40)

    merged_peak_table = merge_peaks_kdtree(
    peak_table_short, peak_table_long,
    tol=3.0,sat_value=65534
    )

    xf, yf = merged_peak_table['x'].data, merged_peak_table['y'].data
    use_flag = merged_peak_table['chosen']
    log.info(f"Found {xf.size} peaks")
    if xf.size > npeaks:
        log.warning("Too many peaks detected")
    elif xf.size < npeaks:
        log.warning(f"{npeaks-xf.size} peaks are not detected")

    #---Calculate Center--------------------------------------------------------------------------------------------------------------
    xobs = np.zeros(npeaks)
    yobs = np.copy(xobs)

    if SaveFiberImage:
        im_crop_full = np.zeros( (npeaks, nwindow*2, nwindow*2))

    for ifiber in range(xf.size):
        ipredict = np.argmin( np.abs(xchip[xf[ifiber]]-xchip) )
        jpredict = np.argmin( np.abs(ychip[yf[ifiber]]-ychip) )

        if use_flag[ifiber]=='t1':
            im_use = im_short
        elif use_flag[ifiber]=='t2':
            im_use = im_long
        im_crop = im_use[jpredict-nwindow:jpredict+nwindow, ipredict-nwindow:ipredict+nwindow]

        x_crop = xchip[ipredict-nwindow:ipredict+nwindow]
        y_crop = ychip[jpredict-nwindow:jpredict+nwindow]

        xobs[ifiber], yobs[ifiber] = com(im_crop, x_crop, y_crop)

        if SaveFiberImage:
            im_crop_full[ifiber] = im_crop

    if SaveFiberImage:
        np.save(data_dir+f'fiberimage.npy', im_crop_full)

    return im, xobs, yobs
