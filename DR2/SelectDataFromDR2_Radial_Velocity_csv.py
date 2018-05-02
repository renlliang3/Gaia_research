import numpy as np
import csv

#readresults650 = []
readresults100 = []

def relevantdata(row):
	reducedrow = [row[5],row[6],row[7],row[8],row[9],row[10],row[11],row[12],row[13],row[14],row[15],row[49],row[50],row[54],row[55],row[59],row[60],row[63],row[66],row[67]]
	reducedrow = [float(i) for i in reducedrow]
	reducedrow.insert(0,int(row[2])) #Insert source_id as integer in front of the row
	return reducedrow

#Too much data for low RAM, read in row by row with the list appending and then use that data.
with open('DR2_Radial_Velocity.csv') as csvdata:
	readcsvdata = csv.reader(csvdata, delimiter=',')
	headerinfo = next(readcsvdata) #Save the headers and then go to next row
	for row in readcsvdata:
		if row[9] and row[11] and row[12] and row[50] and row[55] and row[60]: #Check if the wanted parameters have values
			#if (float(row[9]) > 1000/650) and (float(row[11]) >= 5.0) and ((float(row[50])+5*np.log10(float(row[9]))-10) < 2): #d<650pc and parallax_over_error>=5.0 and M_G < 2
				#readresults650.append(relevantdata(row))
			if (float(row[9]) > 1000/100) and (float(row[11]) > 10.0) and (float(row[54]) > 10.0) and (float(row[59]) > 10.0): #d<100pc, parallax_over_error>10, bp_meanfluxovererror>10
				readresults100.append(relevantdata(row))
				
readresults100.insert(0,[headerinfo[2],headerinfo[5],headerinfo[6],headerinfo[7],headerinfo[8],headerinfo[9],headerinfo[10],headerinfo[11],headerinfo[12],headerinfo[13],headerinfo[14],headerinfo[15],headerinfo[49],headerinfo[50],headerinfo[54],headerinfo[55],headerinfo[59],headerinfo[60],headerinfo[63],headerinfo[66],headerinfo[67]])

with open('Reduced_DR2_Radial_Velocity_<10pc', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
    wr.writerows(readresults100)

#d<650
#1000/parallax < 650
#1000/650 < parallax

