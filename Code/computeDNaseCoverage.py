import sys
import argparse

def normalize(overlap,length):
	for element in overlap:
		if (int(length[element]) != int(0)):
			overlap[element]=overlap[element]/length[element]
	return overlap

def compareChr(dchr,cC):
	if (dchr != "X") and (cC != "X") and  (dchr != "Y") and (cC != "Y"):
		return (int(dchr) > int(cC))
	elif ((dchr == "X") and ((cC != "X") and (cC != "Y"))):
		return True
	elif ((dchr == "Y") and  (cC != "Y")):
		return True
	else:
		return False

def main():
	parser=argparse.ArgumentParser(prog="computeDNaseCoverage.py")
	parser.add_argument("DNase",nargs=1,help="wiggle file of the DNase run. The file has to be sorted.")
	parser.add_argument("Regions",nargs=1,help="bed file containing the regions of interest.")
	args=parser.parse_args()

	overlap=[]
	oD={}
	length={}
	validChr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"}
	dnase=open(args.DNase[0],"r")
	regions=open(args.Regions[0],"r")
	for lr in regions:
		slr=lr.split()
		rchr=slr[0]
		rstart=int(slr[1])
		rend=int(slr[2])
		overlap+=[(rchr,rstart,rend)]
		oD[(rchr,rstart,rend)]=0.0
		length[(rchr,rstart,rend)]=0
	regions.close()
	
	cI=int(0)
	cC=str(overlap[0][0])
	ld=dnase.readline()
	for ld in dnase:
		if ("#" in ld):
			continue
		s=ld.split()
		dchr=str(s[0].replace("chr",""))
		if (dchr in validChr):
			ds=int(s[1])
			de=int(s[2])
			cC=overlap[cI][0]
			while (compareChr(dchr,cC)):
				if (cI == len(overlap)-1):
					break
				cI+=1
				cC=overlap[cI][0]
	
			while (int(overlap[cI][1]) < ds) and (int(overlap[cI][2]) < de) and (ds > int(overlap[cI][2])) and (dchr == cC):
				if (cI == len(overlap)-1):
					break
				cI+=1
				cC=overlap[cI][0]	
			if (dchr == cC):
				if (ds < int(overlap[cI][1])) and (de <= int(overlap[cI][2])) and (de > int(overlap[cI][1])):
					oD[overlap[cI]]+=float(s[3])*(float(de-int(overlap[cI][1]))/(de-ds))
					length[overlap[cI]]+=(de-int(overlap[cI][1]))
				elif (ds < int(overlap[cI][1])) and (de > int(overlap[cI][2])):
					oD[overlap[cI]]+=float(s[3])*(float(int(overlap[cI][2])-int(overlap[cI][1]))/(de-ds))
					length[overlap[cI]]+=(int(overlap[cI][2])-int(overlap[cI][1]))
				elif (ds >= int(overlap[cI][1])) and (de > int(overlap[cI][2])) and (ds < int(overlap[cI][2])):
					oD[overlap[cI]]+=float(s[3])*(float(int(overlap[cI][2])-ds)/(de-ds))
					length[overlap[cI]]+=(int(overlap[cI][2])-ds)
				elif (ds >= int(overlap[cI][1])) and (de <= int(overlap[cI][2])):
					oD[overlap[cI]]+=float(s[3])
					length[overlap[cI]]+=(de-ds)
				else:
					continue
			else:
				continue
	dnase.close()

	oD=normalize(oD,length)
	for element in overlap:
		print(str(element[0])+"\t"+str(element[1])+"\t"+str(element[2])+"\t"+str(oD[element]))
			
main()
