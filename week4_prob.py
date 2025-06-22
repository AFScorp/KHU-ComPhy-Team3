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
	
    result= simulate(p,0,T,dt)
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
    min_R[i] = min_dist

error = 1e-6
delta = 0.001
solution_theta_p=[]
for i in np.arange(1,359,1):			# 각도에 따른 최소거리의 극솟값 찾기

	if min_R[i] < 1e10 and min_R[i-1] > min_R[i] < min_R[i+1]:
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