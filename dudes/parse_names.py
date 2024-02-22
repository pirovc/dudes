def parse_names(names_file,taxids):
	"""
	Generate taxid to scientific name map. Only names of taxids in param taxid are parsed.

	:param names_file:
	:param taxids:
	:return: dict, key: taxid, value: scientific name
	"""
	file = open(names_file,'r')
	names = {}
	for l in file:
		fields = l.split('\t|\t')
		if fields[3].replace('\t|\n','')=="scientific name":
			if int(fields[0]) in taxids:
				names[int(fields[0])] = fields[1]

	file.close()
	return names
