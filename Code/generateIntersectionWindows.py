import sys

def main():
	infile=open(sys.argv[1],"r")
	for l in infile:
		if ("#" not in l):
			s=l.split()
			if (s[2]=="gene") or (s[2]=="start_codon"):
				if (s[6]=="+"):
					if (sys.argv[3]=="FALSE"):
						print(str(s[0].replace("chr",""))+"\t"+str(max(0,int(int(s[3])-int(float(sys.argv[2])/2.0))))+"\t"+str(int(s[3])+int(float(sys.argv[2])/2.0)))
					else:
						print(str(s[0].replace("chr",""))+"\t"+str(max(0,int(int(s[3])-int(int(sys.argv[2])/2.0))))+"\t"+str(s[4]))
					
				else:
					if (sys.argv[3]=="FALSE"):
						print(str(s[0].replace("chr",""))+"\t"+str(max(0,int(int(s[4])-int(int(sys.argv[2])/2.0))))+"\t"+str(int(s[4])+int(int(sys.argv[2])/2.0)))
					else:
						print(str(s[0].replace("chr",""))+"\t"+str(int(s[3]))+"\t"+str(int(s[4])+int(int(sys.argv[2])/2.0)))						
				
	infile.close()


main()
