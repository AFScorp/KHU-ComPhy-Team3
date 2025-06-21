import numpy as np
import matplotlib.pyplot as plt

# 상수 설정
G = 6.67430e-11  # 중력 상수 [m^3/kg/s^2]
s_M = 1.9885e30  # 태양 질량 [kg]
e_M = 5.972e24   # 지구 질량 [kg]
m_M = 7.348e22   # 달 질량 [kg]
AU = 1.496e11    # 천문단위 [m]
oe = 0.0167      # 지구 궤도 이심률

# 미분 방정식 정의
def f(t, p):
	e_x, e_y, e_vx, e_vy = p[0:4]
	m_x, m_y, m_vx, m_vy = p[4:8]
	a_x, a_y, a_vx, a_vy = p[8:12]

	r_es = np.sqrt(e_x**2 + e_y**2)
	r_em = np.sqrt((e_x - m_x)**2 + (e_y - m_y)**2)
	r_ms = np.sqrt(m_x**2 + m_y**2)
	r_as = np.sqrt(a_x**2 + a_y**2)
	r_ae = np.sqrt((a_x - e_x)**2 + (a_y - e_y)**2)
	r_am = np.sqrt((a_x - m_x)**2 + (a_y - m_y)**2)

	e_ax = -G * s_M * e_x / r_es**3 - G * m_M * (e_x - m_x) / r_em**3
	e_ay = -G * s_M * e_y / r_es**3 - G * m_M * (e_y - m_y) / r_em**3
	m_ax = -G * s_M * m_x / r_ms**3 - G * e_M * (m_x - e_x) / r_em**3
	m_ay = -G * s_M * m_y / r_ms**3 - G * e_M * (m_y - e_y) / r_em**3
	a_ax = -G * s_M * a_x / r_as**3 - G * e_M * (a_x - e_x) / r_ae**3 - G * m_M * (a_x - m_x) / r_am**3
	a_ay = -G * s_M * a_y / r_as**3 - G * e_M * (a_y - e_y) / r_ae**3 - G * m_M * (a_y - m_y) / r_am**3

	return np.array([e_vx, e_vy, e_ax, e_ay,
				m_vx, m_vy, m_ax, m_ay,
				a_vx, a_vy, a_ax, a_ay], float)

# 시뮬레이션 시간 설정
T = 60 * 60 * 24 * 366  # 1년
dt = 60 * 60            # 1시간
dt_c = 2				# 2초
t_list = np.arange(0.0, T, dt)
N = len(t_list) + 1

# 초기 조건
# 지구
e_x_i = AU * (1 - oe)
e_y_i = 0
e_vx_i = 0
e_vy_i = np.sqrt(G * s_M * (1 + oe) / (AU * (1 - oe)))

# 달
moon_distance = 3.844e8
m_x_i = e_x_i + moon_distance
m_y_i = 0
m_vx_i = e_vx_i
m_vy_i = e_vy_i + np.sqrt(G * e_M / moon_distance)

# 소행성 입력
a_x_i=float(input("Asteroid initial x(AU)"))*1.496e11
a_y_i=float(input("Asteroid initial y(AU)"))*1.496e11
a_speed_i=float(input("Asteroid initial speed(km/s), ex)earth speed = 29.78 km/s"))*1000
a_theta_mean= np.radians(float(input("Asteroid initial theta (deg): ")))
a_theta_kappa= float(input("Asteroid theta error (kappa), ex) 90% cut +- 5 degree kappa = 67.33: "))

#-180도 ~ 179도 1도 단위로 분할
thetas_d = np.linspace(-180, 179, 360)
thetas_r = np.radians(thetas_d)

# 상태 벡터 초기화
p_base = np.array([e_x_i, e_y_i, e_vx_i, e_vy_i,
		m_x_i, m_y_i, m_vx_i, m_vy_i,
		a_x_i, a_y_i, 0, 0], float)		#a_vx_i, a_vy_i를 아직 안정했기에 0, 0 대입

#360개의 각도에 대한 최소 r_ae 구하기기
R_earth = 6.371e7  # 지구 반지름 [m] 실제 지구 크기의 10배로 설정
min_R_p = []     # 각도별 최소 거리 리스트

for i, theta in enumerate(thetas_r):
	a_vx_i = a_speed_i * np.sin(theta)
	a_vy_i = -a_speed_i * np.cos(theta)

	p = p_base.copy()
	p[10], p[11] = a_vx_i, a_vy_i

	min_R = np.inf  # 최소 거리 무한 설정

	for j, t in enumerate(t_list):
		k1 = dt * f(t, p)
		k2 = dt * f(t + 0.5 * dt, p + 0.5 * k1)
		k3 = dt * f(t + 0.5 * dt, p + 0.5 * k2)
		k4 = dt * f(t + dt, p + k3)
		p += (k1 + 2 * k2 + 2 * k3 + k4) / 6

		# 최소 거리 계산
		r_ae = np.sqrt((p[8] - p[0])**2 + (p[9] - p[1])**2)
		if r_ae < min_R:
			min_R = r_ae

	# 각도별 최소 거리 저장 및 출력
	min_R_p.append(min_R)
	theta_deg = np.degrees(theta)
	print(f"Theta {theta_deg:.2f} deg: min R = {min_R:.2e} m")

# 시각화 마무리
plt.figure()
plt.plot(thetas_d, min_R_p, marker='o')
plt.axhline(R_earth, color='red', linestyle='--', label='Earth Radius')
plt.xlabel("Initial Theta (deg)")
plt.ylabel("Minimum Distance (m)")
plt.title("Minimum Distance vs Initial Launch Angle")
plt.legend(loc='upper right')
plt.grid(True)
plt.show()

#/////////////////////////////////////////////////
def lisim(theta):								#선형 근사가 포함된 최소거리 찾기 함수
	a_vx_i = a_speed_i * np.sin(theta)
	a_vy_i = -a_speed_i * np.cos(theta)

	# 상태 벡터
	p = p_base.copy()
	p[10], p[11] = a_vx_i, a_vy_i

	#상태 저장 리스트 (초기조건)
	e_px, e_py = np.zeros(N), np.zeros(N)
	a_px, a_py = np.zeros(N), np.zeros(N)
	e_px[0], e_py[0] = e_x_i, e_y_i
	a_px[0], a_py[0] = a_x_i, a_y_i

	min_R = np.inf  # 최소 거리 무한 설정
	min_TF = True

	#4th RKM
	for j, t in enumerate(t_list):
		k1 = dt * f(t, p)
		k2 = dt * f(t + 0.5 * dt, p + 0.5 * k1)
		k3 = dt * f(t + 0.5 * dt, p + 0.5 * k2)
		k4 = dt * f(t + dt, p + k3)
		p += (k1 + 2 * k2 + 2 * k3 + k4) / 6

		a_px[j + 1], a_py[j + 1] = p[8], p[9]
		e_px[j + 1], e_py[j + 1] = p[0], p[1]

		r_ae = np.sqrt((p[8] - p[0])**2 + (p[9] - p[1])**2)
		if r_ae < min_R:
			min_R = r_ae
			min_TF = True

		if min_TF == True and r_ae > min_R:
			collision_N = j
			e_vx_c = p[2]
			e_vy_c = p[3]
			a_vx_c = p[10]
			a_vy_c = p[11]
			min_TF = False

	#충돌 예상 부근 정밀 검사
	e_v_c = np.sqrt((e_vx_c)**2 + (e_vy_c)**2)
	a_v_c = np.sqrt((a_vx_c)**2 + (a_vy_c)**2)

	e_x_c_i = e_px[collision_N-2]
	e_y_c_i = e_py[collision_N-2]
	e_x_c_f = e_px[collision_N]
	e_y_c_f = e_py[collision_N]

	a_x_c_i = a_px[collision_N-2]
	a_y_c_i = a_py[collision_N-2]
	a_x_c_f = a_px[collision_N]
	a_y_c_f = a_py[collision_N]

	e_c_p = np.array([e_x_c_i,e_y_c_i])
	a_c_p = np.array([a_x_c_i,a_y_c_i])

	r_e_c = np.sqrt((e_x_c_f - e_x_c_i)**2 + (e_y_c_f - e_y_c_i)**2)
	r_a_c = np.sqrt((a_x_c_f - a_x_c_i)**2 + (a_y_c_f - a_y_c_i)**2)

	e_direc_x = (e_x_c_f - e_x_c_i)/r_e_c
	e_direc_y = (e_y_c_f - e_y_c_i)/r_e_c
	a_direc_x = (a_x_c_f - a_x_c_i)/r_a_c
	a_direc_y = (a_y_c_f - a_y_c_i)/r_a_c

	e_direc = np.array([e_direc_x,e_direc_y])
	a_direc = np.array([a_direc_x,a_direc_y])

	T_e_c = r_e_c/e_v_c
	T_a_c = r_a_c/a_v_c

	if T_e_c > T_a_c:
		T_c = T_e_c
	else:
		T_c = T_a_c

	t_c_list = np.arange(0.0,T_c,dt_c)

	min_c_R = np.inf

	for t in t_c_list:
		e_c_p = e_c_p + e_direc * e_v_c * dt_c
		a_c_p = a_c_p + a_direc * a_v_c * dt_c

		r_ae_c = np.sqrt((e_c_p[0] - a_c_p[0])**2 + (e_c_p[1] - a_c_p[1])**2)

		if r_ae_c < min_c_R:
				min_c_R = r_ae_c
		else:
			break
	print(f"{min_c_R/1000} km")
	return min_c_R
#//////////////////////////////////////////////////////////////////////////////////////////
error = 1e-6
delta = 0.001
solution_theta_p=[]
for i in np.arange(1,359,1):			# 각도에 따른 최소거리의 극솟값 찾기

	if min_R_p[i] < 1e10 and min_R_p[i-1] > min_R_p[i] < min_R_p[i+1]:
		min_R_the = lisim(thetas_r[i])		#극솟값
			
		theta_left = thetas_r[i-1]
		theta_right = thetas_r[i+1]
		theta_middle = thetas_r[i]
		j=0
		while theta_right - theta_left > error:
			min_R_middle = lisim(theta_middle)
			if min_R_middle < R_earth:										#근이 양쪽에 위치
		
				theta_B_left = theta_middle									##Bisection Method로 오른쪽 근 탐색
				theta_B_right = theta_right

				while theta_B_right - theta_B_left > error:
					theta_B_middle = (theta_B_left + theta_B_right)/2
					min_R_middle = lisim(theta_B_middle)

					if min_R_middle < R_earth:
						theta_B_left = theta_B_middle
					else:
						theta_B_right = theta_B_middle
					
				solution_theta_right = np.degrees(theta_B_middle)
					
				theta_B_left = theta_left									#Bisection Method로 왼쪽 근 탐색
				theta_B_right = theta_middle
		
				while theta_B_right - theta_B_left > error:
					theta_B_middle = (theta_B_left + theta_B_right)/2
					min_R_middle = lisim(theta_B_middle)

					if min_R_middle < R_earth:
						theta_B_left = theta_B_middle
					else:
						theta_B_right = theta_B_middle
				
				solution_theta_left = np.degrees(theta_B_middle)

				solution_theta_p.append(solution_theta_left)
				solution_theta_p.append(solution_theta_right)
				break
			else:
				theta_delta = theta_middle + delta/(2**j)
				min_R_delta = lisim(theta_delta)
		
				if min_R_middle > min_R_delta:							#근 2개가 theta_middle 오른쪽에 위치
					theta_left = theta_middle
					theta_middle = (theta_left + theta_right)/2
				else:
					theta_right = theta_middle
					theta_middle = (theta_left + theta_right)/2
				j+=1

solution_theta_p.append(0)
i=0
while solution_theta_p[i] != 0:
	solution_deg_L = solution_theta_p[i]
	solution_deg_R = solution_theta_p[i+1]
	print(f"({solution_deg_L} deg,{solution_deg_R} deg)")			#근 범위 출력
	i += 2

def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result *= i
    return result

def I0(kappa, terms=50):
    sum_val = 0.0
    for n in range(terms):
        num = (kappa / 2) ** (2 * n)
        denom = factorial(n) ** 2
        sum_val += num / denom
    return sum_val

def von_mises_pdf(theta, mu, kappa):
    denom = 2 * np.pi * I0(kappa)
    return np.exp(kappa * np.cos(theta - mu)) / denom

i=0
possibility = 0
while solution_theta_p[i] != 0:
	solution_deg_L = solution_theta_p[i]
	solution_deg_R = solution_theta_p[i+1]
	solution_rad_L = np.radians(solution_deg_L)
	solution_rad_R = np.radians(solution_deg_R)
	h = 1e-5
	for j in np.arange(solution_rad_L,solution_rad_R,h):
		alpha = (von_mises_pdf(j,a_theta_mean,a_theta_kappa) + von_mises_pdf(j+h,a_theta_mean,a_theta_kappa))*h/2
		possibility += alpha
	i += 2

print(f"소행성과 충돌할 확률: {possibility*100} % ")

		



