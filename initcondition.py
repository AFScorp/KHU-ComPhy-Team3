import numpy as np
from generalfns import *
from constants import G
'''
import scipy as sp
import datetime

T0 = datetime.datetime(year=)
'''
def GetVelocity(apsis, ecc, pos, origin_mass, origin=None, axis=None):
    '''
    apsis: 장반경
    ecc: 이심률
    pos: 위치 좌표
    origin_mass: 모항성/모행성 질량
    origin: 모항성/모행성 좌표
    '''
    if origin:
        pos-=origin

    x = pos[0]
    cos = x/norm(pos)

    r = apsis*(1-ecc**2)/(1+ecc*cos)

    speed = np.sqrt(G*origin_mass*(2/r-1/apsis))
    
    return velocity, speed

    
