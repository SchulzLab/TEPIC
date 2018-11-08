import sys

def main():
	infile=open(sys.argv[1],"r")
	counter=0
	tf=""
	for l in infile:
		if (">" in l):
			s=l.split()
			if (tf!=""):
				print(tf+'\t'+str(counter))
				counter=0
			tf=s[1].upper()
		elif ("#" not in l):
			counter+=1
	print(tf+'\t'+str(counter))
	infile.close()


main()
