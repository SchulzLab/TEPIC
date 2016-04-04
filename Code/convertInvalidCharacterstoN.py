import sys

def main():
	infasta=open(sys.argv[1],"r")
	outfasta=open(sys.argv[2],"w")

	for line in infasta:
		if (">" in line):
			outfasta.write(line)
		else:
			outfasta.write(line.replace("R","N").replace("Y","N").replace("r","N").replace("y","N").replace("m","N").replace("M","N").replace("w","N").replace("W","N").replace("s","N").replace("S","N").replace("k","N").replace("K","N").replace("v","N").replace("V","N").replace("h","N").replace("H","N").replace("d","N").replace("D","N").replace("b","N").replace("B","N"))

	infasta.close()
	outfasta.close()

main()
