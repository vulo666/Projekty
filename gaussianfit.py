##############################################
# SCRIPT TO FIT DATA WITH GAUSSIAN FUNCTION ##
##############################################

from scipy.optimize import curve_fit
import string
import numpy as np

def gaussian(x, a, b, c): #defines gaussian
    val = a * np.exp(-(x - b)**2 / c**2)
    return val
# opens data and loads it to MEM to be fitted
vstup = open('names','r')

vstupne = []

for line in vstup:
		vstupne.append(line.split())

vstup.close()

names = []

for i in range(0,len(vstupne)):
	names.append(vstupne[i][0])
	
out1 = open('fit1.dat', 'w')
out2 = open('fit2.dat', 'w')
out3 = open('fit3.dat', 'w')
out4 = open('fit4.dat', 'w')

for name in names:
	imp = open('%s' %name ,'r')
	if len(name) == 16: #lines with peculiar datasets
		time = name[:4]
		shell_number = name[5]
	else:
		time = name[:3] #normal datasets
		shell_number = name[4]
	data = []
	for line in imp:
		data.append(line.split())
	imp.close()
	for i in range(0,len(data)):
		for j in range(0,len(data[i])):
			data[i][j] = float(data[i][j])
	x = []
	y = []
	for i in range(0,len(data)):
		x.append(data[i][0])
		y.append(data[i][1])
	popt, pcov = curve_fit(gaussian, x, y)
	a = popt[0]
	da = np.sqrt(pcov[0,0])
	b = popt[1]
	db = np.sqrt(pcov[1,1])
	c = popt[2]
	dc = np.sqrt(pcov[2,2])
	#print(time, shell_number) JUST TO CHECK!!!
	#print(a, da)
	#print(b, db)
	#print(c, dc)
	riadok = "%s %f %f %f %f %f %f\n" %(time, a, da, b, db, abs(c), dc) #output lines
	if shell_number == '1':
		out1.write(riadok)
	elif shell_number == '2':
		out2.write(riadok)	
	elif shell_number == '3':
		out3.write(riadok)	
	elif shell_number == '4':
		out4.write(riadok)	
	elif shell_number == '5':
		out5.write(riadok)	
	elif shell_number == '6':
		out6.write(riadok)	
	elif shell_number == '7':
		out7.write(riadok)	
	else:
		out8.write(riadok)
	plot = open('%s.gp' %name,'w')
	### OUTPUT READY FOR GNUPLOT
	riadok1 = 'set terminal pdf enhanced\n'
	riadok2 = "set title 'T = %s Myr'\n" %(time)
	riadok3 = "set xlabel '{/Symbol f} [deg]'\n"
	riadok4 = "set ylabel 'N' rotate by 180\n"
	riadok5 = "set output '%s.pdf'\n" %(name)
	riadok6 = "f(x) = %f * exp(-(x-(%f))**2/(%f)**2)\n" %(a, b, c)
	riadok7 = "plot f(x) title 'gaussian fit' lw 2, '%s' title 'data' w boxes" %(name)
	plot.write(riadok1)
	plot.write(riadok2)
	plot.write(riadok3)
	plot.write(riadok4)
	plot.write(riadok5)
	plot.write(riadok6)
	plot.write(riadok7)
	plot.close()
	
out1.close()
out2.close()
out3.close()
out4.close()	
