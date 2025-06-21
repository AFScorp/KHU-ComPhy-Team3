import numpy as np
from simulation import simulate
from constants import *
from pylab import plot,xlabel,ylabel,show,scatter,legend,title

# 시뮬레이션 시간 설정
T = 60 * 60 * 24 * 366  # 1년
dt = 60 * 60            # 1시간

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

a_pos = np.array([a_x0, a_y0])

a_theta_mean= float(input("Asteroid initial direction (deg): "))
a_theta_error= float(input("Asteroid direction maximum error (deg): "))

directions = np.radians(np.linspace(a_theta_mean - a_theta_error, a_theta_mean + a_theta_error,21))

for i, a_theta in enumerate(directions):
    a_vel = a_speed * np.array([np.sin(a_theta), np.sin(a_theta)])

    # 초기값 벡터
    p = np.array([e_pos, m_pos, a_pos, e_vel, m_vel, a_vel], float)

    res = simulate(p, 0, T, dt)
    
    if i==10: #중심선
        plot(res[0,0], res[0,1], 'b', label='Earth') # 지구
        plot(res[1,0], res[1,1], 'r', label='Moon') # 달
        plot(res[3,0], res[3,1], 'k', label='Asteroid', linewidth=2.5)
    else:
        plot(res[3,0], res[3,1], color='gray', alpha=1 - 0.07 * abs(10 - i))

scatter(0, 0, color='yellow', s=200, label='Sun')
title("Sun-Fixed Asteroid Simulation", fontsize=16)
xlabel(r"$x\ (m)$", fontsize=16)
ylabel(r"$y\ (m)$", fontsize=16)
legend(loc='upper right')
show()
    