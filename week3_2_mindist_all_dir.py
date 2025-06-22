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
#a_theta=np.radians(float(input("Asteroid initial direction(deg): ")))

a_pos = np.array([a_x0, a_y0])

thetas = np.linspace(-180,180,181)
min_R = np.zeros(181)
#-180 deg ~ +180 deg, 1 deg diff.
for i, a_theta in enumerate(np.radians(thetas)):
    a_vel = a_speed * np.array([np.cos(a_theta), np.sin(a_theta)])

    #초기값 벡터
    p = np.array([e_pos, m_pos, a_pos, e_vel, m_vel, a_vel], float)
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
    res[i] = min_dist

# 시각화 마무리
plt.plot(np.degrees(thetas), min_R, marker='o')
plt.axhline(e_radius, color='red', linestyle='--', label='Earth Radius')
plt.xlabel("Initial Theta (deg)")
plt.ylabel("Minimum Distance (m)")
plt.title("Minimum Distance vs Initial Launch Angle")
plt.legend(loc='upper right')
plt.grid(True)
plt.show()
