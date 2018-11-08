import sys
import argparse
import os

#This script converts PWMS in JASPAR format to PWMs in TRANSFAC format, such that they can be converted to the TRAP energy format
#arg1 File containing the pwms in Jaspar format

def createHeader(outfile,name,database):
	outfile.write("//\n")
	outfile.write("XX\n")
	outfile.write("ID "+database+"\t"+name+"\n")
	outfile.write("XX\n")
	outfile.write("P0\tA\tC\tG\tT\n")

def storePreviousFactor(outfile,scores):
	motivelength=len(scores)/4
	for i in range(0,motivelength):
		temp=str(i+1)
		for j in range(0,4):
			temp+="\t"+scores[motivelength*j+i]
		outfile.write(temp+"\n")
	outfile.write("XX\n")
			
def main():
	parser=argparse.ArgumentParser(prog="ConvertTrainingDataToMaxAffinityFormat.py")
	parser.add_argument("jaspar",nargs=1,help="File containing the pwms in JASPAR format")
	args=parser.parse_args()
	jasparFormat=open(args.jaspar[0],"r")
	PSCMFormat=open(args.jaspar[0]+".PSCM","w")
	scores=[]
	for l in jasparFormat:
		if (">" in l):
			if (scores!=[]):
				storePreviousFactor(PSCMFormat,scores)
			database=l.split()[0].replace(">","")
			name=l.split()[1]
			createHeader(PSCMFormat,name,database)
			scores=[]
		else:
			scores+=l.split()
	storePreviousFactor(PSCMFormat,scores)
	PSCMFormat.close()
	jasparFormat.close()


main()
