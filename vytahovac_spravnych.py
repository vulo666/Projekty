import string
import numpy as np

vstup = open('slupky_hernquist','r') #polomery slupek

vstupne = []

for line in vstup:
		vstupne.append(line.split())

vstup.close()

names = []

for i in range(0,len(vstupne)):
	names.append(vstupne[i][0])

for i in range(0,len(vstupne)):
	for j in range(0,len(vstupne[i])):
		vstupne[i][j] = float(vstupne[i][j])

k = 0

for name in names:
	imp = open('%s' %(name),'r')   #otvor spravny subor
	data = []
	for line in imp:
		data.append(line.split())
	imp.close()
	for i in range(0,len(data)):
		for j in range(0,len(data[i])):
			data[i][j] = float(data[i][j])
	for N in range(1,4):
		R = vstupne[k][N]				#hodnota polomeru
		if R != 0:
			print(name, N, k, R)
			delta = 0.1*R			#vhodny interval
			out1 = open('%s_%sxyz.dat' %(name, N),'w') #otvor spravny subor
			out2 = open('%s_%sout.dat' %(name, N),'w')
			print(np.shape(data))
			for i in range(0,len(data)):
				x = data[i][1]
				y = data[i][2]
				z = data[i][3]
				vx = data[i][4]
				vy = data[i][5]
				vz = data[i][6]
				r = data[i][0]
				if r>=(R-delta) and  r<=R:
					if N%2==0:
						if x<0:
							#print(x, y, z)
							riadok1 = '%f %f %f %f %f %f\n' %(x, y, z, vx, vy, vz)
							#riadok = '%f %f %f %f\n' %(data[i][0], data[i][1], data[i][2], data[i][3]) #obsolate
							out1.write(riadok1)
							pomer = y/x
							uhol = np.degrees(np.arctan(pomer))
							riadok2 = '%f\n' %(uhol)
							out2.write(riadok2)
					else:
						if x>0:
							#print(x, y, z)
							riadok1 = '%f %f %f %f %f %f\n' %(x, y, z, vx, vy, vz)
							#riadok = '%f %f %f %f\n' %(data[i][0], data[i][1], data[i][2], data[i][3]) #obsolate
							out1.write(riadok1)
							pomer = y/x
							uhol = np.degrees(np.arctan(pomer))
							riadok2 = '%f\n' %(uhol)
							out2.write(riadok2)
			out1.close()
			out2.close()
#################################################################


			imp = open('%s_%sout.dat' %(name, N),'r')   #otvor spravny subor
			vystup = []
			for line in imp:
				vystup.append(line.split())
			imp.close()
			for i in range(0,len(vystup)):
				for j in range(0,len(vystup[i])):
					vystup[i][j] = float(vystup[i][j])
			maximum = int(np.amax(vystup))
			minimum = int(np.amin(vystup))
			output = open("%s_%svystup.dat" %(name, N),'w')
			n = 1. #velkost intervalu
			i = float(minimum)
			while True:
				suma = 0
				for j in range(0,len(vystup)):
					if vystup[j][0]>=i and vystup[j][0]<(i+n):
						suma += 1
				riadok = "%f %i\n" %(i, suma)
				output.write(riadok)
				i = i + n
				if i == maximum:
					break
			output.close()
		if N==4:
			k += 1

