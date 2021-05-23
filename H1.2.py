# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 23:12:05 2021

@author: karam
"""

# -*- coding: utf-8 -*-

flexcode = '''TITLE 'G2 H1.2 gillk34'     { the problem identification }
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
r(threshold=1e-3) = vector(rx,ry)
v(threshold=1e-3) = vector(vx,vy)
 SELECT         { method controls }
ngrid = 1
DEFINITIONS    { parameter definitions }
ag = vector(0, -9.81)

!mi = 600
thetai=45*pi/180
vi = 40

mfueli = 600! 0.8*mi
mdry = 200!mi - mfueli

q = %s
tfuel = mfueli/q

mfuel = if(t<tfuel) then mfueli-q*t else 0
vfuel = 1400
Ftmag = if(t<tfuel) then q*vfuel else 0
Ft = Ftmag*v/magnitude(v)

m = mdry + mfuel

Fg = ag*m


rho = 1.2
area = 2
Cd = 0.8
vwind = vector(-20,0)
vrel = v - vwind
Fdmag = 0.5*rho*Cd*area*(magnitude(vrel))^2
Fd = -Fdmag*vrel/magnitude(vrel)

vrelx = dot(vrel, vector(1,0))
vrely = dot(vrel, vector(0,1))

vhatperp = vector(-vrely, vrelx)/magnitude(vrel)
Areawing = 15
CL = 0.1
FLmag = 0.5*rho*CL*areaWing*(magnitude(vrel))^2
FL = FLmag*vhatperp



Fnet =  Fg + Ft + Fd !+ FL
a = Fnet/m

INITIAL VALUES
v = vi*vector(cos(thetai),sin(thetai))
r = vector(4e3,400)
EQUATIONS        { PDE's, one for each variable }
r: dt(r)=v
v: dt(v) = a

! CONSTRAINTS    { Integral constraints }
BOUNDARIES       { The domain definition }
  REGION 1       { For each material region }
    START(0,0)   { Walk the domain boundary }
    LINE TO (1,0) TO (1,1) TO (0,1) TO CLOSE
TIME 0 TO 200 halt(ry<0)    { if time dependent }
MONITORS         { show progress }
PLOTS            { save result displays }
for t = 0 by endtime/30 to endtime
history(rx,ry) at (0,0) PrintOnly Export Format '#t#b#1#b#2' file = 'test.txt'


summary
report eval(rx, 0,0)
report eval(ry, 0,0)
report tfuel
END'''



FlexFileName = "G2 H1.2 gillk34.pde"
import subprocess

import numpy as np

import matplotlib.pyplot as plt

import time

starttime = time.time()
Bestq = -1
BestRange = -1
FFRs = np.arange(15,66,10)
for q in FFRs:
    with open(FlexFileName, "w") as f:
        print(flexcode%q, file=f)
    #subprocess.run(["FlexPDE6s", "-S",FlexFileName],timeout=5)
   
    completed =subprocess.run(["C:\FlexPDE6student\FlexPDE6s.exe","-S",FlexFileName],timeout=30) #,shell=True)
    #completed = subprocess.run(["FlexPDE6s", "-S",FlexFileName],timeout=30)
    print('returned: ', completed.returncode)

    with open("test.txt") as f:

        flexoutputrawdata=np.loadtxt(f, skiprows=7)
    t = flexoutputrawdata[:,0]
    xd = flexoutputrawdata[:,1]
    yd = flexoutputrawdata[:,2]
    print('{q} kg/s landed after {t} s at x = {x} m.'.format(q=q,t=t[-1],x=xd[-1]))
    plt.plot(xd,yd,label=q)
    if(xd[-1]>BestRange):
        BestRange=xd[-1]
        Bestq=q


plt.xlabel('x-displacement [m]')
plt.ylabel('y-displacement [m]')
plt.legend()
plt.show()
print('Best flow rate was {q} kg/s, which travelled {x} m.'.format(q=Bestq,x=BestRange+3e3))
print('Took {t} seconds.'.format(t=(time.time()-starttime)))

