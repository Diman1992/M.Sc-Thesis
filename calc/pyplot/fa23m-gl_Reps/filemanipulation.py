with open(r'fa23mGlReps.txt','r') as infile, open(r'outputGood.txt','w') as outfile:
	data = infile.read()
	data = data.replace("{","")
	data = data.replace(".}, ","\n")
	data = data.replace("}, ","\n")
	data = data.replace(", ","\t")
	data = data.replace("*^","e")
	outfile.write(data)
