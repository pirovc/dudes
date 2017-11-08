def parse_names(names_file,taxids):
	file = open(names_file,'r')
	names = {}
	for l in file:
		fields = l.split('\t|\t')
		if fields[3].replace('\t|\n','')=="scientific name":
			if int(fields[0]) in taxids:
				names[int(fields[0])] = fields[1]

	file.close()
	return names