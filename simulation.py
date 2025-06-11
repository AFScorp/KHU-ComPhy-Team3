# 시뮬레이션 자체를 담당하는 모듈입니다.
import numpy as np
from generalfns import *
from constants import *
#import time

#time_start = time.time()

# 미분방정식: r''=-G*m_i*(r-r_i)/dist_i (아인슈타인 표기 사용됨)
# suggestion: 추후 일반화를 하게 된다면 array를 3층으로 만들어서 위치 텐서, 속도 텐서, 가속도 텐서 등으로 분리해도 좋을 것 같습니다.

class SimulationManager:
	def init(self):
		return


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
def simulate(p0:np.ndarray, T_start:float, T_end:float, dt:float):
	'''
	p0: array-like, 초기 조건 벡터/텐서

	T: float, 시뮬레이션의 총 주기(s)

	dt: float, 시뮬레이션의 각 프레임별 시간 간격(s)
	'''
	res=RKM4(np.arange(T_start,T_end,dt), p0, f)
	return res