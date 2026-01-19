import numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree


def com(im_crop, x_crop, y_crop):
    xsum = np.sum(im_crop, axis=0)
    ysum = np.sum(im_crop, axis=1)

    xcom = np.sum( x_crop*xsum ) / np.sum(xsum)
    ycom = np.sum( y_crop*ysum ) / np.sum(ysum)

    return xcom, ycom


focal2camera_coeff = np.array([-1.18e-1, 0., 0., 0.
                            , -1e-4 , -1e-8 , -1e-10, -1e-13
                            ,  1e-7 ,  1e-11, -1e-1
                            , -1e-10, -1e-10,  1e-8 , -1e-8])

camera2focal_coeff = np.array([-8.5, 0. , 0., 0.
                            , 1e-6 , 1e-11,  1e-16, 1e-20
                            , 1e-8 , 1e-11, -1e-3 
                            , 1e-14, 1e-14, -1e-10, 1e-10])

def transform_polynomial(xin
               , m
               , theta
               , x0, y0
               , r1x, r2x, r3x, r4x
               , t1x, t2x, t3x
               , a1x, a2y, a4x, a3y
               ):
    
    xtemp, ytemp = xin
    xtemp = xtemp[:int(xtemp.size/2)]; ytemp = ytemp[:int(ytemp.size/2)]

    x = m*(xtemp*np.cos(theta)-ytemp*np.sin(theta))
    y = m*(xtemp*np.sin(theta)+ytemp*np.cos(theta))
    
    r = (x**2 + y**2)**0.5
    xrad_term = x*(r1x*(r**2) + r2x*(r**4) + r3x*(r**6) + r4x*(r**8) )
    xtan_term = (t1x*(r**2 + 2*(x**2)) + t2x*2.*x*y)*(1.+t3x*(r**2))

    xnew = (x + xrad_term + xtan_term
            + a1x*z20(x,y)
            + a4x*z10(x,y)            
            )

    yrad_term = y*(r1x*(r**2) + r2x*(r**4) + r3x*(r**6) + r4x*(r**8))
    ytan_term = (t2x*(r**2 + 2*(y**2)) + t1x*2.*x*y)*(1.+t3x*(r**2))

    ynew = (y + yrad_term + ytan_term
            + a2y*z21(x,y)
            + a3y*z9(x,y)
            )

    new = np.concatenate((xnew-x0, ynew-y0))

    return new

def transform(x, y, coeff):
    xin = np.concatenate((x, x))
    yin = np.concatenate((y, y))
    temp = transform_polynomial((xin, yin)
                        , *coeff
                        )

    xout, yout = temp[:int(temp.size/2)], temp[int(temp.size/2):]

    return xout, yout

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


def find_angle_double_method2_elbow_down(
    ori_x, ori_y,
    target_x_ori, target_y_ori
):
    """
    elbow-down 해만 반환:
      phi(=beta) >= 0 (0~pi)
      theta는 [0, 2pi)로 정규화

    반환: theta, phi (rad)
    """
    arm1 = 5.2; arm2 = 11.6
    Positioner_p = arm1 + arm2

    # base 기준 상대좌표
    target_x = target_x_ori - ori_x
    target_y = target_y_ori - ori_y

    d = float(np.hypot(target_x, target_y))
    gamma = float(np.arctan2(target_y, target_x))  # [-pi, pi]

    # ---- 도달 불가능: 너무 멀면 완전 펼침(phi=0) ----
    if d >= Positioner_p:
        theta = gamma
        phi = 0.0

    # ---- 도달 불가능: 너무 가까우면 완전 접힘(phi=pi) ----
    elif d <= abs(arm2 - arm1):
        theta = gamma + np.pi
        phi = np.pi

    # ---- 도달 가능: elbow-down (phi in [0,pi]) ----
    else:
        # 원래 MATLAB 그대로: phi = beta in [0,pi]
        arg = (d**2 - arm1**2 - arm2**2) / (-2.0 * arm1 * arm2)
        arg = float(np.clip(arg, -1.0, 1.0))

        theta_beta = float(np.arccos(arg))
        phi = float(np.pi - theta_beta)  # <= 이게 beta(0~pi)

        theta_alpha = float(np.arctan2(
            arm2 * np.sin(phi),
            arm1 + arm2 * np.cos(phi)
        ))
        theta = float(gamma - theta_alpha)

    # ---- theta를 항상 [0, 2pi)로 정규화 ----
    theta = float(theta % (2.0 * np.pi))

    # ---- phi는 elbow-down이라서 항상 [0,pi] 유지 ----
    # (수치 오차로 아주 미세하게 -가 나올 가능성만 방지)
    phi = float(np.clip(phi, 0.0, np.pi))

    return theta, phi

def z2(x,y):
    return x

def z3(x,y):
    return y

def z4(x, y):
    return 2.*(x**2 + y**2) - 1.

def z5(x, y):
    return x*y

def z6(x, y):
    return (x**2 - y**2)

def z7(x,y):
    return y*(3.*(x**2 + y**2) - 2.)

def z8(x,y):
    return x*(3.*(x**2 + y**2) - 2.)

def z9(x,y):
    return y*(3.*(x**2) - y**2)

def z10(x,y):
    return  x*(x**2 - 3.*(y**2))

def z11(x,y):
    return 6.*( (x**2+y**2)**2 ) - 6.*(x**2+y**2) +1

def z12(x,y):
    return (x**2 - y**2)*(4.*(x**2+y**2) - 3.)

def z13(x,y):
    return x*y*(4.*(x**2+y**2) - 3.)

def z14(x,y):
    return (x**2+y**2)**2 - 8.*(x**2)*(y**2)

def z15(x,y):
    return x*y*(x**2 - y**2) 

def z16(x,y):
    return x*(10.*((x**2+y**2)**2) - 12.*(x**2+y**2) +3.)

def z17(x,y):
    return y*(10.*((x**2+y**2)**2) - 12.*(x**2+y**2) +3.)

def z18(x,y):
    return x*(x**2 - 3.*(y**2))*(5.*(x**2+y**2)-4.)

def z19(x,y):
    return y*(3.*(x**2) - y**2)*(5.*(x**2+y**2)-4.)

def z20(x,y):
    return x*(16.*(x**4) - 20.*(x**2)*(x**2+y**2) + 5.*((x**2+y**2)**2))

def z21(x,y):
    return y*(16.*(y**4) - 20.*(y**2)*(x**2+y**2) + 5.*((x**2+y**2)**2))
