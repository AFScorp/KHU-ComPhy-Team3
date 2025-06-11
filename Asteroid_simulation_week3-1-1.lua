from numpy import arange,array,sqrt,sin,cos,radians
from pylab import plot,xlabel,ylabel,show,scatter,legend,title

#상수 설정
G=6.67430e-11 #G[m^3/kg/s^2]
s_M=1.9885e30   #[Kg]
e_M = 5.972e24
m_M = 7.348e22
AU = 1.496e11   #[m]
oe = 0.0167  # 지구 궤도의 이심률 orbital_eccentricity

#변수 벡터화, 미방 함수
def f(t,p):
	e_x=p[0]
	e_y=p[1]
	e_vx=p[2]
	e_vy=p[3]
	m_x=p[4]
	m_y=p[5]
	m_vx=p[6]
	m_vy=p[7]
	a_x=p[8]
	a_y=p[9]
	a_vx=p[10]
	a_vy=p[11]
	r_es=sqrt(e_x**2+e_y**2)
	r_em=sqrt((e_x-m_x)**2+(e_y-m_y)**2)
	r_ms=sqrt(m_x**2+m_y**2)
	r_as=sqrt(a_x**2+a_y**2)
	r_ae=sqrt((a_x-e_x)**2+(a_y-e_y)**2)
	r_am=sqrt((a_x-m_x)**2+(a_y-m_y)**2)
	e_ax=-G*s_M*e_x/(r_es**3)-G*m_M*(e_x-m_x)/(r_em**3)
	e_ay=-G*s_M*e_y/(r_es**3)-G*m_M*(e_y-m_y)/(r_em**3)
	m_ax=-G*s_M*m_x/(r_ms**3)-G*e_M*(m_x-e_x)/(r_em**3)
	m_ay=-G*s_M*m_y/(r_ms**3)-G*e_M*(m_y-e_y)/(r_em**3)
	a_ax=-G*s_M*a_x/(r_as**3)-G*e_M*(a_x-e_x)/(r_ae**3)-G*m_M*(a_x-m_x)/(r_am**3)
	a_ay=-G*s_M*a_y/(r_as**3)-G*e_M*(a_y-e_y)/(r_ae**3)-G*m_M*(a_y-m_y)/(r_am**3)
	return array([e_vx,e_vy,e_ax,e_ay,m_vx,m_vy,m_ax,m_ay,a_vx,a_vy,a_ax,a_ay],float)

#시간 변수
T=60*60*24*366  #1 year = 365.25636 day
dt1=60*60
dt2=60	###
t_i,t_f = 0.0,T
t_list1=arange(t_i,t_f,dt1)
t_list2=arange(t_i,t_f,dt2)	###

#지구 변수
e_x_i,e_y_i = AU*(1-oe),0
e_vx_i,e_vy_i=0,sqrt(G*s_M*(1+oe)/(AU*(1-oe)))

#달 변수
moon_distance = 3.844e8
m_x_i,m_y_i = AU*(1-oe)+moon_distance,0
m_vx_i,m_vy_i=e_vx_i,e_vy_i+sqrt(G * e_M / moon_distance)

#소행성 변수
a_x_i=float(input("Asteroid initial x(AU)"))*1.496e11
a_y_i=float(input("Asteroid initial y(AU)"))*1.496e11
a_v_i=float(input("Asteroid initial speed(km/s), ex)earth speed = 29.78 km/s"))*1000
a_theta_i=radians(float(input("Asteroid initial theta(dgree)")))
a_vx_i=a_v_i*sin(a_theta_i)
a_vy_i=a_v_i*-cos(a_theta_i)
a_px1=[a_x_i]
a_py1=[a_y_i]
a_px2=[a_x_i]
a_py2=[a_y_i]

p=array([e_x_i,e_y_i,e_vx_i,e_vy_i,m_x_i,m_y_i,m_vx_i,m_vy_i,a_x_i,a_y_i,a_vx_i,a_vy_i],float)

#4th RKM
for t in t_list1:
	k1=dt1*f(t,p)
	k2=dt1*f(t+0.5*dt1,p+0.5*k1)
	k3=dt1*f(t+0.5*dt1,p+0.5*k2)
	k4=dt1*f(t+dt1,p+k3)
	p+=(k1+2.0*k2+2.0*k3+k4)/6.0
	a_px1.append(p[8])
	a_py1.append(p[9])

p=array([e_x_i,e_y_i,e_vx_i,e_vy_i,m_x_i,m_y_i,m_vx_i,m_vy_i,a_x_i,a_y_i,a_vx_i,a_vy_i],float)		###

for t in t_list2:		###아래는 전부 새 코드
	k1=dt2*f(t,p)
	k2=dt2*f(t+0.5*dt2,p+0.5*k1)
	k3=dt2*f(t+0.5*dt2,p+0.5*k2)
	k4=dt2*f(t+dt2,p+k3)
	p+=(k1+2.0*k2+2.0*k3+k4)/6.0
	a_px2.append(p[8])
	a_py2.append(p[9])

error_px=[0]
error_py=[0]

for n in arange(1,366*24):
	dx=a_px2[60*n]-a_px1[n]
	error_px.append(dx)
	dy=a_py2[60*n]-a_py1[n]
	error_py.append(dy)

#시각화
plot(t_list1, error_px, 'r', label='x error')
plot(t_list1, error_py, 'b', label='y error')
xlabel("Time (seconds)")
ylabel("Position Error (m)")
legend(loc='upper right')
title("RK4 Position Error of Asteroid (dt=1hr vs dt=1min)")
show()


