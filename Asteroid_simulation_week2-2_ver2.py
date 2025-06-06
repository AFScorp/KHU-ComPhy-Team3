import numpy as np
from pylab import plot, xlabel, ylabel, show, scatter, legend, title

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

# 소행성
a_x_i=float(input("Asteroid initial x(AU)"))*1.496e11
a_y_i=float(input("Asteroid initial y(AU)"))*1.496e11
a_speed_i=float(input("Asteroid initial speed(km/s), ex)earth speed = 29.78 km/s"))*1000
a_theta_i= np.radians(float(input("Asteroid initial theta (deg): ")))

a_vx_i = a_speed_i * np.sin(a_theta_i)
a_vy_i = -a_speed_i * np.cos(a_theta_i)

# 상태 벡터 초기화
p = np.array([e_x_i, e_y_i, e_vx_i, e_vy_i,
		m_x_i, m_y_i, m_vx_i, m_vy_i,
		a_x_i, a_y_i, a_vx_i, a_vy_i], float)

# 상태 저장 리스트
e_px, e_py = np.zeros(N), np.zeros(N)
m_px, m_py = np.zeros(N), np.zeros(N)
a_px, a_py = np.zeros(N), np.zeros(N)
e_px[0], e_py[0] = e_x_i, e_y_i
m_px[0], m_py[0] = m_x_i, m_y_i
a_px[0], a_py[0] = a_x_i, a_y_i

# 4th Runge-Kutta Method
i=0
for t in t_list:
	k1 = dt * f(t, p)
	k2 = dt * f(t + 0.5 * dt, p + 0.5 * k1)
	k3 = dt * f(t + 0.5 * dt, p + 0.5 * k2)
	k4 = dt * f(t + dt, p + k3)
	p += (k1 + 2 * k2 + 2 * k3 + k4) / 6

	e_px[i + 1], e_py[i + 1] = p[0], p[1]
	m_px[i + 1], m_py[i + 1] = p[4], p[5]
	a_px[i + 1], a_py[i + 1] = p[8], p[9]
	i += 1

# 시각화
scatter(0, 0, color='yellow', s=200, label='Sun')
plot(e_px, e_py, 'b', label='Earth')
plot(m_px, m_py, 'r', label='Moon')
plot(a_px, a_py, 'k', label='Asteroid')
title("Sun-Fixed Asteroid Simulation", fontsize=16)
xlabel(r"$x\ (m)$", fontsize=16)
ylabel(r"$y\ (m)$", fontsize=16)
legend(loc='upper right')
show()

