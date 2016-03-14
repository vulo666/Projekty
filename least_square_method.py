############################################################
# Program to fit observed data with least-square method#####
############################################################


import math as m
import scipy, scipy.linalg
import string
import numpy as np

cech = open('cech.txt','r') #otvori import subot s datami

data = []
#this whole part convets data from string to float
for line in cech:
    data.append(string.split(line))

cech.close()

n = len(data)

for i in range(0,n):
    for j in range(0,2):
        data[i][j]=string.atof(data[i][j])
#till this

x_list=[]
y_list=[]

for i in range(0,n): #makes matrices
    x = data[i][0]
    x_list.append([x**2,x,1])
    y_list.append(data[i][1])

x = scipy.matrix(x_list)
y = scipy.matrix(y_list)

y = scipy.transpose(y)
y_cech = y


v = scipy.transpose(x)*x
u = scipy.transpose(x)*y
h = np.linalg.inv(v)
b = h*u
r = scipy.transpose(y)*y-scipy.transpose(b)*u
sigma = np.sqrt(r/(n-3))
delta_b = sigma*np.sqrt(np.diagonal(h))
delta_yp = sigma*np.sqrt(np.diagonal(x*h*scipy.transpose(x)))

#output from model
print 'cech bez vah'
print 'a1: {} +- {}' .format(b[0,0],delta_b[0,0])
print 'a2: {} +- {}' .format(b[1,0],delta_b[0,1])
print 'a4: {} +- {}' .format(b[2,0],delta_b[0,2])
print 'sigma: {}' .format(float(sigma))

a1 = b[0,0]
a2 = b[1,0]
t_min = -a2/(2*a1)
dtda1 = a2/(2*a1**2) #derivative tmin / a1
dtda2 = -1/(2*a1) # derivative tmin / a2
dt = np.matrix([dtda1,dtda2,0]) #gradient tmin
delta_t_min = sigma*np.sqrt(dt*h*scipy.transpose(dt))

print 't_min: {} +- {}' .format(t_min,float(delta_t_min))

a6 = b[2,0]-b[0,0]*t_min**2
da6 = np.matrix([t_min**2,t_min,1]) #gradient a6
delta_a6 = sigma*np.sqrt(da6*h*scipy.transpose(da6))

print 'a6: {} +- {}' .format(float(a6),float(delta_a6))
print '\n'

cinan = open('cinan.txt','r')

data = []
#this whole part converts string to float
for line in cinan:
    data.append(string.split(line))

cinan.close()

n = len(data)

for i in range(0,n):
    for j in range(0,2):
        data[i][j]=string.atof(data[i][j])
#till this
x_list=[]
y_list=[]

for i in range(0,n):
    x = data[i][0]
    x_list.append([x**2,x,1])
    y_list.append(data[i][1])

x = scipy.matrix(x_list)
y = scipy.matrix(y_list)

y = scipy.transpose(y)
y_cinan = y


v = scipy.transpose(x)*x
u = scipy.transpose(x)*y
h = np.linalg.inv(v)
b = h*u
yp = x*b
r = scipy.transpose(y)*y-scipy.transpose(b)*u
sigma = np.sqrt(r/(n-3))
delta_b = sigma*np.sqrt(np.diagonal(h))

#output fomr model
print 'cinan bez vah'
print 'a1: {} +- {}' .format(b[0,0],delta_b[0,0])
print 'a2: {} +- {}' .format(b[1,0],delta_b[0,1])
print 'a3: {} +- {}' .format(b[2,0],delta_b[0,2])
print 'sigma: {}' .format(float(sigma))

a1 = b[0,0]
a2 = b[1,0]
t_min = -a2/(2*a1)
dtda1 = a2/(2*a1**2) #derivative tmin / a1
dtda2 = -1/(2*a1) # derivative tmin / a2
dt = np.matrix([dtda1,dtda2,0]) #gradient tmin
delta_t_min = sigma*np.sqrt(dt*h*scipy.transpose(dt))

print 't_min: {} +- {}' .format(t_min,float(delta_t_min))

a5 = b[2,0]-b[0,0]*t_min**2
da5 = np.matrix([t_min**2,t_min,1]) #gradient a5
delta_a5 = sigma*np.sqrt(da5*h*scipy.transpose(da5))

print 'a5: {} +- {}' .format(float(a5),float(delta_a5))
print '\n'


############################################################
############################################################
##############Joint Observations###########################
############################################################
############################################################
cech = open('cech.txt','r')
cinan = open('cinan.txt','r')

data = []
#Converts string to float
for line in cech:
    data.append(string.split(line))

cech.close()

n = len(data)

for i in range(0,n):
    for j in range(0,2):
        data[i][j]=string.atof(data[i][j])
#till this

x_list=[]
y_list=[]

for i in range(0,n):
    x = data[i][0]
    x_list.append([x**2,x,0,1])
    y_list.append(data[i][1])

data = []
#converts sting to float
for line in cinan:
    data.append(string.split(line))

cinan.close()

n = len(data)

for i in range(0,n):
    for j in range(0,2):
        data[i][j]=string.atof(data[i][j])
#till this

for i in range(0,n):
    x = data[i][0]
    x_list.append([x**2,x,1,0])
    y_list.append(data[i][1])

n = len(x_list)

x = scipy.matrix(x_list)
y = scipy.matrix(y_list)

y = scipy.transpose(y)


v = scipy.transpose(x)*x
u = scipy.transpose(x)*y
h = np.linalg.inv(v)
b = h*u
yp = x*b
r = scipy.transpose(y)*y-scipy.transpose(b)*u
sigma = np.sqrt(r/(n-4))
delta_b = sigma*np.sqrt(np.diagonal(h))

yp_cech = np.zeros((30,1))
yp_cinan = np.zeros((15,1))

for i in range(0,30):
    yp_cech[i][0] = yp[i][0]

for i in range(0,15):
    yp_cinan[i][0] = yp[i+30][0]

yp_cinan = scipy.matrix(yp_cinan)
yp_cech = scipy.matrix(yp_cech)

r_cech = scipy.transpose(y_cech)*y_cech-scipy.transpose(yp_cech)*yp_cech
r_cinan = scipy.transpose(y_cinan)*y_cinan-scipy.transpose(yp_cinan)*yp_cinan

sigma_cech = np.sqrt(r_cech/(30-3))
sigma_cinan = np.sqrt(r_cinan/(15-3))


print 'spojene bez vah'

print 'a1: {} +- {}' .format(b[0,0],delta_b[0,0])
print 'a2: {} +- {}' .format(b[1,0],delta_b[0,1])
print 'a3: {} +- {}' .format(b[2,0],delta_b[0,2])
print 'a4: {} +- {}' .format(b[3,0],delta_b[0,3])
print 'sigma: {}' .format(float(sigma))

a1 = b[0,0]
a2 = b[1,0]
t_min = -a2/(2*a1)
dtda1 = a2/(2*a1**2) #derivative tmin / a1
dtda2 = -1/(2*a1) # derivative tmin / a2
dt = np.matrix([dtda1,dtda2,0,0]) #gradient tmin
delta_t_min = sigma*np.sqrt(dt*h*scipy.transpose(dt))

print 't_min: {} +- {}' .format(t_min,float(delta_t_min))

a6 = b[3,0]-b[0,0]*t_min**2
da6 = np.matrix([t_min**2,t_min,0,1]) #gradient a6
delta_a6 = sigma*np.sqrt(da6*h*scipy.transpose(da6))
a5 = b[2,0]-b[0,0]*t_min**2
da5 = np.matrix([t_min**2,t_min,1,0]) #gradient a5
delta_a5 = sigma*np.sqrt(da5*h*scipy.transpose(da5))

print 'a5: {} +- {}' .format(float(a5),float(delta_a5))
print 'a6: {} +- {}' .format(float(a6),float(delta_a6))
print '\n'

w_vektor_cech = np.zeros((30,1))
w_vektor_cinan = np.zeros((15,1))


for i in range(0,30):
    w_vektor_cech[i][0] = (45*m.pow(sigma_cech,-2))/(30*m.pow(sigma_cech,-2)+15*m.pow(sigma_cinan,-2))

for i in range(0,15):
    w_vektor_cinan[i][0] = (45*m.pow(sigma_cinan,-2))/(30*m.pow(sigma_cech,-2)+15*m.pow(sigma_cinan,-2))

w_vektor = np.vstack((w_vektor_cech,w_vektor_cinan))

n = len(w_vektor)
w = np.zeros((n,n))
for i in range(0,n):
    w[i][i] = w_vektor[i][0]

w = scipy.matrix(w)

cech = open('cech.txt','r')
cinan = open('cinan.txt','r')

data = []

for line in cech:
    data.append(string.split(line))

cech.close()

n = len(data)

for i in range(0,n):
    for j in range(0,2):
        data[i][j]=string.atof(data[i][j])

x_list=[]
y_list=[]

for i in range(0,n):
    x = data[i][0]
    x_list.append([x**2,x,0,1])
    y_list.append(data[i][1])

data = []

for line in cinan:
    data.append(string.split(line))

cinan.close()

n = len(data)

for i in range(0,n):
    for j in range(0,2):
        data[i][j]=string.atof(data[i][j])

for i in range(0,n):
    x = data[i][0]
    x_list.append([x**2,x,1,0])
    y_list.append(data[i][1])

n = len(x_list)

x = scipy.matrix(x_list)
y = scipy.matrix(y_list)

y = scipy.transpose(y)


v = scipy.transpose(x)*w*x
u = scipy.transpose(x)*w*y
h = np.linalg.inv(v)
b = h*u
b_w = b
h_w = h
r = scipy.transpose(y)*w*y-scipy.transpose(b)*u
sigma = np.sqrt(r/(n-4))
sigma_w = sigma
delta_b = sigma*np.sqrt(np.diagonal(h))

print 'spojene s vahami'

print 'a1: {} +- {}' .format(b[0,0],delta_b[0,0])
print 'a2: {} +- {}' .format(b[1,0],delta_b[0,1])
print 'a3: {} +- {}' .format(b[2,0],delta_b[0,2])
print 'a4: {} +- {}' .format(b[3,0],delta_b[0,3])
print 'sigma: {}' .format(float(sigma))

a1 = b[0,0]
a2 = b[1,0]
t_min = -a2/(2*a1)
dtda1 = a2/(2*a1**2) #derivative tmin / a1
dtda2 = -1/(2*a1) # derivative tmin / a2
dt = np.matrix([dtda1,dtda2,0,0]) #gradient tmin
delta_t_min = sigma*np.sqrt(dt*h*scipy.transpose(dt))

print 't_min: {} +- {}' .format(t_min,float(delta_t_min))

a6 = b[3,0]-b[0,0]*t_min**2
da6 = np.matrix([t_min**2,t_min,0,1]) #gradient a6
delta_a6 = sigma*np.sqrt(da6*h*scipy.transpose(da6))
a5 = b[2,0]-b[0,0]*t_min**2
da5 = np.matrix([t_min**2,t_min,1,0]) #gradient a5
delta_a5 = sigma*np.sqrt(da5*h*scipy.transpose(da5))

print 'a5: {} +- {}' .format(float(a5),float(delta_a5))
print 'a6: {} +- {}' .format(float(a6),float(delta_a6))
print '\n'


################################################################
#######################ERROR BARS#############################
################################################################

x = np.linspace(0,0.8,1000)

f1 = open('out_pasy1.dat','w')
f2 = open('out_pasy2.dat','w')

#calculates errorbars and outputs it into file for graphs to be made
for i in range(0,len(x)):
    g1_list = []
    g2_list = []
    g1_list.append([x[i]**2,x[i],0,1])
    g2_list.append([x[i]**2,x[i],1,0])
    g1 = scipy.matrix(g1_list)
    g2 = scipy.matrix(g2_list)
    y1 = g1*b_w
    y2 = g2*b_w
    tmp = g1*h_w*scipy.transpose(g1)
    delta_y1 = float(sigma_w*np.sqrt(g1*h_w*scipy.transpose(g1)))
    delta_y2 = float(sigma_w*np.sqrt(g2*h_w*scipy.transpose(g2)))
    riadok1 = "%f %f %f\n" %(x[i],y1,delta_y1)
    riadok2 = "%f %f %f\n" %(x[i],y2,delta_y2)
    f1.write(riadok1)
    f2.write(riadok2)

f1.close
f2.close
