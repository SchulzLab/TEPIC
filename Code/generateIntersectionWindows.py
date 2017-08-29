import sys

def main():
	if (".gz" in sys.argv[1]):
		print("Gzipped files are not supported")
		return 0;
	infile=open(sys.argv[1],"r")
	identifier="start_codon"
	for l in infile:
		if ("#" not in l):
			s=l.split()
			if (s[2]=="gene"):
				identifier="gene"
				break
	infile.close()
	genebody=sys.argv[3]

	if ((identifier =="start_codon") and (sys.argv[3]!="FALSE")):
		genebody="FALSE"

	infile=open(sys.argv[1],"r")
	for l in infile:
		if ("#" not in l):
			s=l.split()
			if (s[2]==identifier):
				if (s[6]=="+"):
					if (genebody=="FALSE"):
						print(str(s[0].replace("chr",""))+"\t"+str(max(0,int(int(s[3])-int(float(sys.argv[2])/2.0))))+"\t"+str(int(s[3])+int(float(sys.argv[2])/2.0)))
					else:
						print(str(s[0].replace("chr",""))+"\t"+str(max(0,int(int(s[3])-int(int(sys.argv[2])/2.0))))+"\t"+str(s[4]))
					
				else:
					if (genebody=="FALSE"):
						print(str(s[0].replace("chr",""))+"\t"+str(max(0,int(int(s[4])-int(int(sys.argv[2])/2.0))))+"\t"+str(int(s[4])+int(int(sys.argv[2])/2.0)))
					else:
						print(str(s[0].replace("chr",""))+"\t"+str(int(s[3]))+"\t"+str(int(s[4])+int(int(sys.argv[2])/2.0)))						
				
	infile.close()


main()
