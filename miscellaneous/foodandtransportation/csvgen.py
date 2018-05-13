import csv 
import re
import os


finalString = ""
files = os.listdir("foodNetworkOutput")
for textFile in files:
	temp=""
	out = ""
	filePath = "foodNetworkOutput/"+textFile
	for line in open(filePath,'r'):
	    out = out + line
	temp+=out
	temp = temp[temp.index("Graph Stats for Max-Clique:"):]
	allNumbers = [re.search(r'\n\SV(.*?)\n',temp).group().replace('\n','')]
	allNumbers.append(re.search(r'\n\SE(.*?)\n',temp).group().replace('\n',''))
	allNumbers.append(re.search(r'\nd_max(.*?)\n',temp).group().replace('\n',''))
	allNumbers.append(re.search(r'\nd_avg(.*?)\n',temp).group().replace('\n',''))
	allNumbers.append(re.search(r'\n\ST(.*?)\n',temp).group().replace('\n',''))
	allNumbers.append(re.search(r'\nT_avg(.*?)\n',temp).group().replace('\n',''))
	allNumbers.append(re.search(r'\nT_max(.*?)\n',temp).group().replace('\n',''))
	allNumbers.append(re.search(r'cc_avg(.*?)\n',temp).group().replace('\n',''))
	allNumbers.append(re.search(r'cc_global(.*?)\n',temp).group().replace('\n',''))

	allNumbersFinal = []
	for x in allNumbers:
	     allNumbersFinal.append(x[x.index(":")+1:].replace(" ",""))

		
	finalString += str(allNumbersFinal).replace("[","").replace("]","").replace("'","").replace('"',"")+"\n"



file = open("foodnetworks.csv","w")
file.write("|V|,|E|,d_max,d_avg,|T|,T_avg,T_max,cc_avg,cc_global\n")
file.write(finalString)
file.close()