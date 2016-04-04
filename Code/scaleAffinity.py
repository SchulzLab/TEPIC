import sys
import argparse

def main():
	parser=argparse.ArgumentParser(prog="scaleAffinity.py")
	parser.add_argument("signalScales",nargs=1,help="File containing scaling information")
	parser.add_argument("affinity",nargs=1,help="TAP seperated TF Affinity file")
	args=parser.parse_args() 
	affinities=open(args.affinity[0],"r")
	print("\t"+affinities.readline().strip())
	scales=open(args.signalScales[0],"r")
	for sl in scales:
		if (len(sl.split())!=1):
			factor=float(sl.split()[3])
		else:
			factor=float(sl.split()[0])
		al=affinities.readline()
		als=al.split()
		newstring=str(als[0])
		for i in xrange(1,len(als)):
			newstring+="\t"+str(float(als[i])*factor)
		print newstring	

	affinities.close()
	scales.close()


main()
