import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from scipy import special
from scipy.ndimage.interpolation import rotate
import pyfits as pf

A = 0.0001 #parametre
B = 43000.0
M = 7000

zeta2 = B/M
ro1 = A*M*((2.0*m.pi*zeta2)**(1.5))
G = 6.67e-11

def f(v,r): #Fce Kingovho modela
    u = v[0]
    y = v[1]
    pre = -4.0*m.pi*G*(ro1/zeta2)
    ex = m.exp(y)
    sqr = m.sqrt(y)
    err = m.erf(sqr)
    f0 = pre*(ex*err-2.0*sqr*m.sqrt(1.0/m.pi)*(1.0+(2.0/3.0)*y))-(2.0/r)*u
    f1 = u
    return [f0,f1]

u0 = 0 #okrajove podmienky
y0 = 9.50
v0 = [u0,y0]

r_t = 6920 #koncovy bod integracie
r = np.linspace (1e-10, r_t, 4610)

yy = odeint(f, v0, r, mxstep = 100000)
y = yy[:,1]

expy = np.exp(y)
sqrty = np.sqrt(y)

ro = ro1*(expy*special.erf(sqrty)-2.0*sqrty*(1.0/m.sqrt(m.pi))*(1.0+(2.0/3.0)*y))
np.clip(ro, 0, 1e20, out=ro)

delta_x = r[1]-r[0]

R = np.linspace (0,r_t-(0.01*r_t), 4610)
delta_R = R[1]-R[0]

tmp_ro = ro
tmp_r = r

epsylon = np.array([])
                      
for i in range(0,len(R)):
	multiplier = tmp_r/(((tmp_r*tmp_r)-(R[i]*R[i]))**0.5)
	integrand = tmp_ro*multiplier
	g = 2.0*trapz(integrand,dx=delta_x)
	epsylon = np.append(epsylon, g)
	tmp_ro = np.delete(tmp_ro, 0)
	tmp_r = np.delete(tmp_r, 0)
epsylon = epsylon*2.5e-8

table1 = np.zeros ((3500,3000))

for a in range(0,3500):
    for b in range(0,3000):
        R_i = round(m.sqrt(a*a+b*b))
        table1[a,b] = epsylon[R_i]

table2 = np.fliplr(table1)
table3 = np.hstack((table2,table1))
table4 = np.flipud(table3)
table = np.vstack((table4,table3))

rottab = np.zeros((4200,6000))

matrix = np.array([[0.6,0],[0,1]]) 
u = np.array([[0],[0]])

for i in range(0,7000):
    for j in range(0,6000):
       v = np.array([[i],[j]])
       u = np.dot(matrix,v)
       rottab[round(u[0,0]),round(u[1,0])] = table[i,j]
       
newtab = np.zeros((4239,4212))
for a in range(0,4239):
    for b in range(0,4212):
        newtab[a,b] = otocena[a+1013,b+1574]

hdu = pf.PrimaryHDU(newtab)
hdu.writeto('namodelovane6920.fits')

hdulist1 = pf.open('j8fk01010_drz.fits')
scidata1 = hdulist1[1].data
sci = rotate(scidata1, 33)

odcitana = scidata1 - newtab
hdu = pf.PrimaryHDU(odcitana)
hdu.writeto('odcitana6920.fits')
