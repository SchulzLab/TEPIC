import sys
#Checks whether a region is within the chromosomes 1 to 22 or X

def isValid(prefix):
	prefix=prefix.replace("chr","")
	vset=set(["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"])
	return prefix in vset

def main():
	inp=open(sys.argv[1],"r")
	outp=open(sys.argv[1].split(".bed")[0]+"_Filtered_Regions.bed","w")
	for l in inp:
		if (isValid(l.split()[0])):
			outp.write(l.replace("chr",""))
	inp.close()
	outp.close()


main()
