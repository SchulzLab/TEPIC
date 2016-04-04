import sys

#argv1 Affinity file
#argv2 Filtered Affinity file
#argv3 pValue Score file
#argv4 Filtered pValue Score file

def isValidAffinity(lineSplit):
	for i in xrange(1,len(lineSplit)):
		if (str(lineSplit[i]) != "0.0"):
			return True	
	return False

def isValidpValue(lineSplit):
	for i in xrange(1,len(lineSplit)):
		if (str(lineSplit[i]) != "1"):
			return True	
	return False


def main():
	#Checking Affinity
	infile=open(sys.argv[1],"r")
	output=open(sys.argv[2],"w")
	#Copy header line
	output.write(infile.readline())
	#Check individual lines
	for l in infile:
		if (isValidAffinity(l.split())):
			output.write(l)
	infile.close()
	output.close()

	if (len(sys.argv) > 3):
		#Checking pValue
		infile=open(sys.argv[3],"r")
		output=open(sys.argv[4],"w")
		#Copy header line
		output.write(infile.readline())
		#Check individual lines
		for l in infile:
			if (isValidpValue(l.split())):
				output.write(l)
		infile.close()
		output.close()

main()
