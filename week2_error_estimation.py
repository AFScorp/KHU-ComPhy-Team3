import numpy as np
from simulation import simulate
from constants import *
from pylab import plot,xlabel,ylabel,show,scatter,legend,title

T = year #1 year = 365.25636 day
dt1 = hour
dt2 = minute	###
t_i,t_f = 0.0,T
t_list1=np.arange(t_i,t_f,dt1)
t_list2=np.arange(t_i,t_f,dt2)

#초기값 설정
#지구
e_pos = np.array([AU*(1-e_ecc), 0])
e_vel = np.array([0,np.sqrt(G*s_M*(1+e_ecc)/(AU*(1-e_ecc)))])

#달
m_pos = np.array([m_radius, 0]) + e_pos
m_vel = np.array(0, np.sqrt(G * e_M / m_radius)) + e_vel

#소행성 변수
a_x0=float(input("Asteroid initial x(AU): "))*AU
a_y0=float(input("Asteroid initial y(AU): "))*AU

a_speed=float(input("Asteroid initial speed(km/s)\nex)earth speed = 29.78 km/s\n: "))*1000
a_theta=np.radians(float(input("Asteroid initial direction(deg): ")))

a_pos = np.array([a_x0, a_y0])
a_vel = a_speed * np.array([np.sin(a_theta), np.sin(a_theta)])

# 초기값 벡터
p = np.array([e_pos, m_pos, a_pos, e_vel, m_vel, a_vel], float)

# 4th order RKM
result_hour = simulate(p, 0, T, dt1)[:,2]
result_minute = simulate(p, 0, T, dt2)[:,2]

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
