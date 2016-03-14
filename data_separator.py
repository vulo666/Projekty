##################################################
# SKRIPT KTORY VYBERA DATA V SPRAVNYCH CASOCH ####
##################################################

import string
import numpy as np

for subor in range(1,14): #mnozstvo parcialnych suborov + 1
	inputname = 'formatedHQ.dat'
	imp = open(inputname,'r')
	print(inputname)
	
	data = []
	
	for line in imp:
		data.append(line.split())
		
	imp.close()
	print('done')
	for i in range(0,len(data)):
		for j in range(0,len(data[i])):
			data[i][j] = float(data[i][j])
			
	print('done')
	
	N = 10001                        #pocet riadkov kazdeho output suboru + 1
	nr_of_files = int(len(data)/N)    #pocet suborov
	
	for i in range(0,nr_of_files):
		beg_line = i*N       #pociatocny riadok = riadok s casom
		#print(data[beg_line])
		filename = str(int(data[beg_line][0])) #nazov output suboru s polohami a rychlostami
		f = open(filename,'w')
		histo = filename + 'histo.dat' #nazov output suboru s histogramom na urcenie poloh slupiek
		histogram = open(histo,'w')
		histo_data = np.zeros((1000,2))
		for j in range(0,len(histo_data)):
			histo_data[j][0] = j
		for j in range(beg_line+1,beg_line+N):  #indexuje riadky
			riadok = "%f %f %f %f %f %f %f\n" %(data[j][0], data[j][1], data[j][2], data[j][3], data[j][4], data[j][5], data[j][6])
			f.write(riadok)
			R = int(data[j][0])
			X = data[j][1]
			if R < 1000:
				histo_data[R][1] += 1
		for l in range(0,len(histo_data)):
			riadok = "%i %i\n" %(histo_data[l][0], histo_data[l][1])
			histogram.write(riadok)
		f.close()
		histogram.close()
