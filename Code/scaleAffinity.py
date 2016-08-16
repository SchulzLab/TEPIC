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
	for al in affinities:
		als=al.split()
		while(True):
			sl=scales.readline()
			s=sl.split()
			if ("#" not in sl):
				factor=float(s[3])
				pos=str(s[0])+":"+str(s[1])+"-"+str(s[2])
			else:
				pos="NA"
			if(pos == str(als[0])):
				newstring=str(als[0])
				for i in xrange(1,len(als)):
					newstring+="\t"+str(float(als[i])*factor)
				print newstring	
				break;


	affinities.close()
	scales.close()


main()
