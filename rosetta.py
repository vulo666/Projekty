######################################
# Simulates movement of one particle #
### in 3D Plummer sphere. Run with ###
######## output to file like: ########
### python rosetta.py > output.dat ###

import numpy as np

n=9000 #nr of timesteps

#constants in potential
GM=10.0
b=5.0
b2=b**2.0

#initial conditions
x=100.0
y=50.0
z=0.0
dini=np.sqrt(x**2.0+y**2.0+z**2.0) #initial distance
vx=-0.5*np.sqrt(2*GM/(np.sqrt(dini**2.0+b2))) #escape velocity at init. distance 
vy=0.0
vz=0.0

def acc(pos, r2): #acceleration
	accel = (- GM * pos) / (r2 + b2)**1.5
	return accel
	
#inicialization at t=0
r2 = x**2.0 + y**2.0
ax = acc(x,r2)
ay = acc(y,r2)
az = acc(z,r2)
vx = vx - 0.5*ax
vy = vy - 0.5*ay
vz = vz - 0.5*az

#for timesteps
for t in range(0,n):
	r2 = x**2.0 + y**2.0 + z**2.0
	ax = acc(x,r2)
	ay = acc(y,r2)
	az = acc(z,r2)
	vx = vx + ax
	vy = vy + ay
	vz = vz + az
	x = x + vx
	y = y + vy
	z = z + vz
	print(x, y, vx, vy)
