import numpy as np
from numpy.linalg import norm
from simulation import simulate
from constants import *
import matplotlib.pyplot as plt

T = year #1 year = 365.25636 day
dt = hour
dt_c = 2 #2초

#초기값 설정
#지구
e_pos = np.array([AU*(1-e_ecc), 0])
e_vel = np.array([0,np.sqrt(G*s_M*(1+e_ecc)/(AU*(1-e_ecc)))])

#달
m_pos = np.array([m_radius, 0]) + e_pos
m_vel = np.array([0, np.sqrt(G * e_M / m_radius)]) + e_vel

#소행성 변수
a_x0=float(input("Asteroid initial x(AU): "))*AU
a_y0=float(input("Asteroid initial y(AU): "))*AU

a_speed=float(input("Asteroid initial speed(km/s)\nex)earth speed = 29.78 km/s\n: "))*1000
a_theta=np.radians(float(input("Asteroid initial direction(deg): ")))

a_pos = np.array([a_x0, a_y0])
a_vel = a_speed * np.array([np.cos(a_theta), np.sin(a_theta)])

#초기값 벡터
p = np.array([e_pos, m_pos, a_pos, e_vel, m_vel, a_vel], float)

#거리 구하기

#시뮬레이션 결과로부터 지구와 천체 간의 거리 계산
result = simulate(p,0,T,dt)
e_pos = result[:,0]
m_pos = result[:,1]
a_pos = result[:,2]
a_pos_rel = a_pos - e_pos #지구로부터의 상대 위치
#a_vel_rel = result[:,5]-result[:,3] #지구로부터의 상대 속도
min_idx = np.argmin(norm(a_pos_rel, axis=1))

#최근점 부근에서 선형 등속 운동을 하는 것으로 interpolation
#속력은 최근점에서 그 근방 점까지 사이의 평균 속력
#1 hr/2s=1800
before_minimum = norm(np.linspace(a_pos_rel[min_idx-1], a_pos_rel[min_idx],1801), axis=0)
after_minimum =  norm(np.linspace(a_pos_rel[min_idx],a_pos_rel[min_idx+1],1801), axis=0)
min_dist = min(np.min(before_minimum),np.min(after_minimum))

print("최소거리: ", min_dist/1000, "km")

# 시각화
plt.plot(a_pos[:,0], a_pos[:,1], 'k', linewidth=2.5, label="Asteroid")
plt.plot(e_pos[:,0], e_pos[:,1], 'b', label='Earth')
plt.plot(m_pos[:,0], m_pos[:,1], 'r', label='Moon')
plt.scatter(0, 0, color='yellow', s=200, label='Sun')
plt.scatter(e_pos[min_idx,0], e_pos[min_idx,1], color='blue', s=200)
plt.scatter(a_pos[min_idx,0], a_pos[min_idx,1], color='Red', s=200)
plt.title("Sun-Fixed Asteroid Simulation", fontsize=16)
plt.xlabel(r"$x\ (m)$", fontsize=16)
plt.ylabel(r"$y\ (m)$", fontsize=16)
plt.legend(loc='upper right')
plt.show()
