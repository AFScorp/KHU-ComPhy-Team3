# 시뮬레이션 자체를 담당하는 모듈입니다.

import numpy as np
from numpy import radians, sqrt
from pylab import plot,xlabel,ylabel,show,scatter,legend,title
from generalfns import *
from constants import *
#import time

#time_start = time.time()

# 미분방정식: r''=-G*m_i*(r-r_i)/dist_i (아인슈타인 표기 사용됨)
# suggestion: 추후 일반화를 하게 된다면 array를 3층으로 만들어서 위치 텐서, 속도 텐서, 가속도 텐서 등으로 분리해도 좋을 것 같습니다.
def f(t,p):
	origin = np.array([0,0]) #태양
	e_pos = p[0]
	m_pos = p[1]
	a_pos = p[2]
	e_vel = p[3]
	m_vel = p[4]
	a_vel = p[5]

	e_acc = gravitation(e_pos, origin, s_M) + gravitation(e_pos, m_pos, m_M)
	m_acc = gravitation(m_pos, origin, s_M) + gravitation(m_pos, e_pos, e_M)
	a_acc = gravitation(a_pos, origin, s_M) + gravitation(a_pos, e_pos, e_M) + gravitation(a_pos, m_pos, m_M)

	return np.array([e_vel, m_vel, a_vel, m_acc, e_acc, a_acc], float)

# 외부에서 실행할 수 있도록 시뮬레이션 자체를 메서드화
def simulate(p0, T, dt):
	'''
	p0: array-like, 초기 조건 벡터/텐서
	T: float, 시뮬레이션의 총 주기
	dt: float, 시뮬레이션의 각 프레임별 시간 간격
	'''
	t_array = np.arange(0, T, dt)
	return RKM4(t_array, p0, f)

'''T = year #1 year = 365.25636 day
dt1 = hour
dt2 = minute	###
t_i,t_f = 0.0,T
t_list1=np.arange(t_i,t_f,dt1)
t_list2=np.arange(t_i,t_f,dt2)

#초기값 설정
#지구
e_pos = np.array([AU*(1-oe), 0])
e_vel = np.array([0,sqrt(G*s_M*(1+oe)/(AU*(1-oe)))])

#달
m_pos = np.array([moon_distance, 0]) + e_pos
m_vel = np.array(0, sqrt(G * e_M / moon_distance)) + e_vel

#소행성 변수
a_x0=float(input("Asteroid initial x(AU): "))*AU
a_y0=float(input("Asteroid initial y(AU): "))*AU

a_speed=float(input("Asteroid initial speed(km/s), ex)earth speed = 29.78 km/s"))*1000
a_theta=radians(float(input("Asteroid initial theta(dgree)")))

a_pos = np.array([a_x0, a_y0])
a_vel = a_speed * np.array([np.sin(a_theta), np.sin(a_theta)])

# 초기값 벡터
p = np.array([e_pos, m_pos, a_pos, e_vel, m_vel, a_vel], float)

# 4th order RKM
result_hour = RKM4(t_list1, p, f)[:,2]
result_minute = RKM4(t_list2, p, f)[:,2]

# 오차 계산
error = np.zeros_like(result_hour)

len_min = len(result_minute)
for i in range(len(result_hour)-1):
	idx_min = 60*i
	if idx_min > len_min:
		break
	error[i] = result_hour[i] - result_minute[60*i]

#시각화
plot(t_list1[:len(error)], error[:,0], 'r', label='x error')
plot(t_list1[:len(error)], error[:,1], 'b', label='y error')
xlabel("Time (s)")
ylabel("Position Error (m)")
legend(loc='upper right')
title("RK4 Position Error of Asteroid (dt=1hr vs dt=1min)")
#print("소요시간:", time.time()-time_start, "초")
show()
'''