import sys

def main():
	infile=open(sys.argv[1],"r")
	for l in infile:
		if ("#" not in l):
			s=l.split()
			if (s[6]=="+"):
				print(str(s[0].replace("chr",""))+"\t"+str(int(int(s[3])-(int(sys.argv[2])/2.0)))+"\t"+str(int(int(s[3])+(int(sys.argv[2])/2.0))))
			else:
				print(str(s[0].replace("chr",""))+"\t"+str(int(int(s[4])-(int(sys.argv[2])/2.0)))+"\t"+str(int(int(s[4])+(int(sys.argv[2])/2.0))))
	infile.close()


main()
