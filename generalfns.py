import numpy as np
from constants import G

# 일반 함수 정의
# Runge-Kutta method, 4th
def RKM4(x, y0, function):
	'''
	Input:
	x: array-like, 시간 구간
	y0: array-like, 초기조건
	function: function, 미분방정식의 음함수꼴

	초기조건에 대하여 구간 위에서 4차 Runge-Kutta 방법으로 방정식을 해결합니다.

	Output:
	y: ndarray, 해의 함수값
	'''
	def eval(t,p):
		return function(t,p)
	
	y=np.zeros((len(x),*y0.shape))
	y[0]=y0
	h=x[1]-x[0] #시뮬레이션의 시간 구간이 일정하다고 가정

	for i in range(len(x)-1):
		t,p=x[i],y[i]
		k1=h*eval(t,p)
		k2=h*eval(t+0.5*h,p+0.5*k1)
		k3=h*eval(t+0.5*h,p+0.5*k2)
		k4=h*eval(t+h,p+k3)
		p+=(k1+2.0*k2+2.0*k3+k4)/6.0
		y[i+1]=p
	
	return y

def norm(x):
	'''
	x: ndarray

	x의 제곱근 노름을 구합니다.
	'''
	res=0.0
	for x_i in x:
		res+=x_i**2
	return np.sqrt(res)

def gravitation(p1,p2,m2):
	pos_rel=p1-p2
	dist=norm(pos_rel)
	if dist <= 0.0001:
		return 0
	else:
		return -G*m2*pos_rel/dist**3