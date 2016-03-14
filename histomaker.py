###########################################################
## PROGRAM TO CREATE RADIAL HISTOGRAMS OF THE DATA FILES###
###########################################################

import matplotlib.pyplot as plt
import numpy as np

f = open('2000','r')

vstupne = []

for line in f:
		vstupne.append(line.split())

f.close()

for i in range(0,len(vstupne)):
	for j in range(0,len(vstupne[i])):
		vstupne[i][j] = float(vstupne[i][j])
Rgas = []
Rlessgas = []
for i in range(0,len(vstupne)):
	r=vstupne[i][0]
	if r<=1000:
		Rgas.append(r)
	if r<=100:
		Rlessgas.append(r)

f = open('/home/kubo/skola/diplomka/merger/beta_ver/hvezdy/hq/2000','r')

vstupne = []

for line in f:
		vstupne.append(line.split())

f.close()

for i in range(0,len(vstupne)):
	for j in range(0,len(vstupne[i])):
		vstupne[i][j] = float(vstupne[i][j])
Rstars = []
Rlessstars = []
for i in range(0,len(vstupne)):
	r=vstupne[i][0]
	if r<=1000:
		Rstars.append(r)
	if r<=100:
		Rlessstars.append(r)
#RHq.append(399) #to make HQ data go to 400 and to have the same bin size as PL

g = plt.figure(1)
plt.hist(Rlessgas, bins=100, histtype='stepfilled', color='b', label='stars')
plt.hist(Rlessstars, bins=100, histtype='step', color='r', label='gas')
plt.xlabel("r [kpc]")
plt.ylabel("N")
plt.legend()
g.savefig('histozoom.pdf', bbox_inches='tight')


f = plt.figure(2)
plt.hist(Rstars, bins=10, histtype='stepfilled', color='b', label='stars')
plt.hist(Rgas, bins=10, histtype='step', color='r', label='gas')
plt.xlabel("r [kpc]")
plt.ylabel("N")
plt.legend()
f.savefig('histo.pdf', bbox_inches='tight')
