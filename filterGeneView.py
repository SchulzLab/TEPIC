import sys

#Removes genes for which no score good be computed

def isValidAffinity(lineSplit):
	for i in xrange(1,len(lineSplit)):
		if (str(lineSplit[i]) != "0.0"):
			return True	
	return False

def main():
	#Checking Affinity
	infile=open(sys.argv[1],"r")
	output=open(sys.argv[1].replace(".txt","_Filtered.txt"),"w")
	#Copy header line
	output.write(infile.readline())
	#Check individual lines
	for l in infile:
		if (isValidAffinity(l.split())):
			output.write(l)
	infile.close()
	output.close()

main()
