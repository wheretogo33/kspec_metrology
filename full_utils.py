import numpy as np
from astropy.table import Table
xccd_global = (np.arange(-int(11760/2)+1, int(11760/2)+1)-0.5)*3.76e-3
yccd_global = (np.arange(-int(8842/2)+1, int(8842/2)+1)-0.5)*3.76e-3

def com(im_crop, x_crop, y_crop):
    xsum = np.sum(im_crop, axis=0)
    ysum = np.sum(im_crop, axis=1)

    xcom = np.sum( x_crop*xsum ) / np.sum(xsum)
    ycom = np.sum( y_crop*ysum ) / np.sum(ysum)

    return xcom, ycom

focal2camera_coeff = np.array([-1.18e-1, 0., 0., 0.
                            , -1e-4 , -1e-8 , -1e-10, -1e-13
                            ,  -1e-7 ,  1e-10, 1e-3, 1e-4
                            , -1e-10, -1e-10,  1e-8 , -1e-8
                            ])

camera2focal_coeff = np.array([-8.5, 0. , 0., 0. 
                            , 1e-6 , 1e-11,  1e-16, 1e-20
                            , 1e-8 , 1e-11, -1e-3, 1e-7
                            , 1e-14, 1e-14, -1e-10, 1e-10
                            ])

def z9(x,y):
    return y*(3.*(x**2) - y**2)

def z10(x,y):
    return  x*(x**2 - 3.*(y**2))

def z20(x,y):
    return x*(16.*(x**4) - 20.*(x**2)*(x**2+y**2) + 5.*((x**2+y**2)**2))

def z21(x,y):
    return y*(16.*(y**4) - 20.*(y**2)*(x**2+y**2) + 5.*((x**2+y**2)**2))


def distrad_SingleTelescope(xin
               , m
               , theta
               , x0, y0
               , r1
               , r2
               , r3
               ):

    xtemp, ytemp = xin
    xtemp = xtemp[:int(xtemp.size/2)]; ytemp = ytemp[:int(ytemp.size/2)] 

    x = m*xtemp
    y = m*ytemp
    
    r = (x**2 + y**2)**0.5
    xrad_term = x*(r1*(r**2) + r2*(r**4) + r3*(r**6) )
    xnew = (x + xrad_term)
    
    yrad_term = y*(r1*(r**2) + r2*(r**4) + r3*(r**6) )
    ynew = (y + yrad_term) 
        
    new = np.concatenate(((xnew*np.cos(theta)-ynew*np.sin(theta))-x0, (xnew*np.sin(theta)+ynew*np.cos(theta))-y0))
    
    return new

def transform_polynomial(xin
               , m
               , theta
               , x0, y0
               , r1x, r2x, r3x, r4x
               , t1x, t2x, t3x, t4x
               , a1x, a2y, a4x, a3y
               ):
    
    xtemp, ytemp = xin
    xtemp = xtemp[:int(xtemp.size/2)]
    ytemp = ytemp[:int(ytemp.size/2)]

    x = m*(xtemp*np.cos(theta)-ytemp*np.sin(theta)) 
    y = m*(xtemp*np.sin(theta)+ytemp*np.cos(theta)) 
    
    #x = m*xtemp
    #y = m*ytemp

    r = (x**2 + y**2)**0.5
    xrad_term = x*(r1x*(r**2) + r2x*(r**4) + r3x*(r**6) + r4x*(r**8))
    xtan_term = (t1x*(r**2 + 2*(x**2)) + t2x*2.*x*y)*(1.+ t3x*(r**2)+ t4x*(r**4))#+ t5x*(r**6))

    xnew = (x + xrad_term + xtan_term
            + a1x*z20(x,y)
            + a4x*z10(x,y)         
            )

    yrad_term = y*(r1x*(r**2) + r2x*(r**4) + r3x*(r**6) + r4x*(r**8))
    ytan_term = (t2x*(r**2 + 2*(y**2)) + t1x*2.*x*y)*(1.+ t3x*(r**2)+ t4x*(r**4))#+ t5x*(r**6))

    ynew = (y + yrad_term + ytan_term
            + a2y*z21(x,y)
            + a3y*z9(x,y)
            )

    #new = np.concatenate((xnew*np.cos(theta)-ynew*np.sin(theta) - x0, xnew*np.sin(theta)+ynew*np.cos(theta) - y0))
    new = np.concatenate((xnew - x0, ynew - y0))

    return new

def transform(x, y, coeff):
    xin = np.concatenate((x, x))
    yin = np.concatenate((y, y))
    temp = transform_polynomial((xin, yin)
                        , *coeff
                        )

    xout, yout = temp[:int(temp.size/2)], temp[int(temp.size/2):]

    return xout, yout


def transform_single(x, y, coeff):
    xin = np.concatenate((x, x))
    yin = np.concatenate((y, y))
    temp = distrad_SingleTelescope((xin, yin)
                        , *coeff
                        )

    xout, yout = temp[:int(temp.size/2)], temp[int(temp.size/2):]

    return xout, yout


from scipy.spatial import cKDTree

def Match_Fiber(x, y
                , xobs_raw, yobs_raw
                , nbuffer):

    coeff_temp = camera2focal_coeff
    xobs, yobs = transform(xobs_raw, yobs_raw, coeff_temp)

    nhunt = 720
    theta_grid = np.linspace(0., 2.*np.pi, nhunt)
    dd_sum = np.zeros(nhunt)
    for ihunt, theta_temp in enumerate(theta_grid):
        xobs_rot = np.cos(theta_temp)*xobs - np.sin(theta_temp)*yobs
        yobs_rot = np.sin(theta_temp)*xobs + np.cos(theta_temp)*yobs
        obs_flag = np.concatenate( (np.full(xobs.size, 0), np.full(xobs.size, 1)) )
        pos_tot = np.concatenate( (np.vstack((x, y)).T
                                 , np.vstack((xobs_rot, yobs_rot)).T) )

        tree = cKDTree(pos_tot)
        dd, ii = tree.query(pos_tot, k=nbuffer)

        for ipeak in range(xobs.size):
            dd_sum[ihunt] += dd[ipeak][obs_flag[ii[ipeak]] == 1].min()

    theta_guess = theta_grid[dd_sum.argmin()]

    xobs_rot = np.cos(theta_guess)*xobs - np.sin(theta_guess)*yobs
    yobs_rot = np.sin(theta_guess)*xobs + np.cos(theta_guess)*yobs

    ngrid = 41
    offset_grid = np.linspace(-10., 10., ngrid)
    dsum_temp = np.zeros((ngrid, ngrid))
    for i in range(ngrid):
        for j in range(ngrid):
            tree = cKDTree(np.c_[xobs_rot+offset_grid[i], yobs_rot+offset_grid[j]])
            d, _ = tree.query(np.c_[x, y], k=1)

            dsum_temp[i,j] = d.sum()

    imin, jmin = np.unravel_index(dsum_temp.argmin(), dsum_temp.shape)

    xobs_rot += offset_grid[imin]
    yobs_rot += offset_grid[jmin]
    #==============================================================================


    pos_tot = np.concatenate( (np.vstack((x, y)).T
                                 , np.vstack((xobs_rot, yobs_rot)).T) )
    tree = cKDTree(pos_tot)
    dd, ii = tree.query(pos_tot, k=nbuffer)   

    imatch = np.zeros(xobs.size, dtype=np.int32)
    for ipeak in range(xobs.size):
        imatch[ipeak] = ii[ipeak][obs_flag[ii[ipeak]] == 1][dd[ipeak][obs_flag[ii[ipeak]] == 1].argmin()] - xobs.size

    if np.unique(imatch).size != xobs.size:
        print("WARNING : Matching is not Unique")

    return imatch, theta_guess, (offset_grid[imin], offset_grid[jmin])



def dedupe_peaks_kdtree(peak_table: Table, min_dist: float) -> Table:
    """
    peak_table: photutils.find_peaks 결과 (x_peak, y_peak, peak_value 컬럼 포함)
    min_dist  : 같은 소스군으로 간주할 최소 거리(픽셀)
    반환값    : 중복 제거된 Table (원본 메타 유지)
    """

    # 배열 추출
    x = np.asarray(peak_table['x_peak'], dtype=float)
    y = np.asarray(peak_table['y_peak'], dtype=float)
    v = np.asarray(peak_table['peak_value'], dtype=float)

    n = x.size
    if n == 0:
        return peak_table.copy()

    # 더 밝은 순서(내림차순)로 처리 → 이 순서로 선택하면 주변 약한 피크를 지움
    order = np.argsort(-v)

    pts = np.column_stack((x, y))
    tree = cKDTree(pts)

    keep = np.ones(n, dtype=bool)
    selected_idx = []

    for idx in order:
        if not keep[idx]:
            continue
        # 이 피크는 채택
        selected_idx.append(idx)
        # min_dist 이내 이웃(자기 자신 포함)
        neighbors = tree.query_ball_point(pts[idx], r=min_dist)
        # 자신보다 어두운 이웃들을 제거(우리는 내림차순으로 돌고 있으므로 이웃은 항상 같거나 더 어둠)
        for nb in neighbors:
            if nb == idx:
                continue
            keep[nb] = False

    selected_idx = np.array(selected_idx, dtype=int)
    selected_idx.sort()  # 원하면 원래 순서로 복원하려면 주석 해제/유지 선택
    return peak_table[selected_idx]


def scatter_contour(x, y,
                    levels=10,
                    threshold=100,
                    log_counts=False,
                    histogram2d_args=None,
                    plot_args=None,
                    contour_args=None,
                    filled_contour=True,
                    ax=None,
                    normalize=None  # <<< NEW: None | 'fraction' | 'max'
                    ):
    """
    Scatter plot with contour over dense regions
    (…중략…)

    Extra
    -----
    normalize : {'fraction', 'max', None}, optional
        'fraction' -> H를 전체 점수로 나눠 각 빈이 차지하는 비율(0~1).
                      threshold와 levels는 0~1 비율로 해석.
        'max'      -> H를 H.max()로 나눠 0~1로 스케일만 통일.
        None       -> 원래대로(카운트 단위).
    """
    x = np.asarray(x)
    y = np.asarray(y)

    default_contour_args = dict(zorder=2)
    default_plot_args = dict(marker='.', linestyle='none', zorder=1)

    if plot_args is not None:
        default_plot_args.update(plot_args)
    plot_args = default_plot_args

    if contour_args is not None:
        default_contour_args.update(contour_args)
    contour_args = default_contour_args

    if histogram2d_args is None:
        histogram2d_args = {}

    if ax is None:
        from matplotlib import pyplot as plt
        ax = plt.gca()

    # 히스토그램 (density=False가 기본이 되도록)
    histogram2d_args = {'density': False, **histogram2d_args}
    H, xbins, ybins = np.histogram2d(x, y, **histogram2d_args)

    # --- NEW: 정규화 옵션 ---
    if normalize == 'fraction':
        total = H.sum()
        if total > 0:
            H = H / total
        # fraction 모드에서는 log 변환을 비활성화(의미가 약함)
        log_counts = False
    elif normalize == 'max':
        hmax = H.max()
        if hmax > 0:
            H = H / hmax

    # (옵션) 로그 변환 (원래 동작 유지)
    if log_counts:
        H = np.log10(1 + H)
        # threshold가 카운트 기준일 때를 가정한 원래 코드.
        # 필요하면 사용자가 log 공간 임계값을 직접 전달하면 됩니다.
        threshold = np.log10(1 + threshold)

    levels = np.asarray(levels)
    if levels.size == 1:
        # threshold ~ H.max() 사이를 균등 분할
        levels = np.linspace(threshold, H.max(), int(levels))

    extent = [xbins[0], xbins[-1], ybins[0], ybins[-1]]
    i_min = np.argmin(levels)

    # 외곽 폴리곤 잡기
    outline = ax.contour(H.T, levels[i_min:i_min + 1],
                         linewidths=0, extent=extent,
                         alpha=0)

    # 등고선/채움
    if filled_contour:
        contours = ax.contourf(H.T, levels, extent=extent, **contour_args)
    else:
        contours = ax.contour(H.T, levels, extent=extent, **contour_args)

    # 외곽선 안쪽은 산점도 생략
    X = np.hstack([x[:, None], y[:, None]])
    if len(outline.allsegs[0]) > 0:
        outer_poly = outline.allsegs[0][0]
        try:
            from matplotlib.path import Path
            points_inside = Path(outer_poly).contains_points(X)
        except ImportError:
            import matplotlib.nxutils as nx
            points_inside = nx.points_inside_poly(X, outer_poly)
        Xplot = X[~points_inside]
    else:
        Xplot = X

    points = ax.plot(Xplot[:, 0], Xplot[:, 1], **plot_args)
    return points, contours 