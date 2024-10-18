import os
import shutil
import subprocess
import math as m
import numpy as np

def cell2prism(cell):
    a=np.linalg.norm(cell[0][:])
    b=np.linalg.norm(cell[1][:])
    c=np.linalg.norm(cell[2][:])
    al=m.degrees(m.acos((np.dot(cell[1][:],cell[2][:]))/(b*c)))
    bt=m.degrees(m.acos((np.dot(cell[0][:],cell[2][:]))/(a*c)))
    gm=m.degrees(m.acos((np.dot(cell[0][:],cell[1][:]))/(a*b)))
    return np.array([a, b, c, al, bt, gm])

def frac2cart(pos, cell):
    return np.matmul(pos, cell)

def cart2frac(pos, cell):
    return np.matmul(pos, np.linalg.inv(cell))

def prism2cell(pri):
    a  = pri[0]
    b  = pri[1]
    c  = pri[2]
    al = m.radians(pri[3])
    bt = m.radians(pri[4])
    gm = m.radians(pri[5])
    sal = m.sin(al)
    cal = m.cos(al)
    sbt = m.sin(bt)
    cbt = m.cos(bt)
    sgm = m.sin(gm)
    cgm = m.cos(gm)

    aa = a
    ab = 0.0
    ac = 0.0
    ba = b*cgm
    bb = b*sgm
    bc = 0.0
    ca = c*cbt
    cb = c*((cal-cbt*cgm)/(sgm))
    cc = (c/sgm)*(m.sqrt(1+2*cal*cbt*cgm-cal**2-cbt**2-cgm**2))
    return np.array([[aa, ab, ac],[ba, bb, bc],[ca, cb, cc]])

def Rx(ang):
    sss=m.sin(m.radians(ang))
    ccc=m.cos(m.radians(ang))
    return np.array([[1.0, 0.0, 0.0],[0.0, ccc, -1*sss],[0.0, sss, ccc]])

def Ry(ang):
    sss=m.sin(m.radians(ang))
    ccc=m.cos(m.radians(ang))
    return np.array([[ccc, 0.0, sss],[0.0, 1.0, 0.0],[-1*sss, 0.0, ccc]])

def Rz(ang):
    sss=m.sin(m.radians(ang))
    ccc=m.cos(m.radians(ang))
    return np.array([[ccc, -1*sss, 0.0],[sss, ccc, 0.0],[0.0, 0.0, 1.0]])

def Rg(vec,ang):
    vl=np.linalg.norm(vec)
    I=np.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]])
    #give a identity matrix if vl is 0, its a bit of a singularity here
    if vl < 1.0e-10:
        return I
    w=np.array([[0, -1.0*(vec[2]/vl), (vec[1]/vl)],[(vec[2]/vl), 0, -1.0*(vec[0]/vl)],[-1.0*(vec[1]/vl), (vec[0]/vl), 0]])
    s1=m.sin(m.radians(ang))
    s2=2.0*(m.sin(m.radians(ang)/2.0))**2.0
    return I+s1*w+s2*np.matmul(w,w)
    


