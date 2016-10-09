with open(r'outfa23m-gl/outfa23_m-gl.txt','r') as infile, open(r'outfa23m-gl/outputGood.txt','w') as outfile:
	data = infile.read()
	data = data.replace("{","")
	data = data.replace(".}, ","\n")
	data = data.replace("}, ","\n")
	data = data.replace(", ","\t")
	data = data.replace("*^","e")
	outfile.write(data)
