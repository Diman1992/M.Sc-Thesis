with open(r'outputpts.txt','r') as infile, open(r'outputGood.txt','w') as outfile:
	data = infile.read()
	data = data.replace("{","")
	data = data.replace(".}, ","\n")
	data = data.replace("}, ","\n")
	data = data.replace(", ","\t")
	data = data.replace("*^","e")
	data = data.replace("}}","")
	outfile.write(data)
