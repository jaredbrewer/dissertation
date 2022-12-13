import bibtexparser as btx
import os, random, string
def bibparser(input, output):

	with open(input) as bibtex_file:
		bib = btx.load(bibtex_file)

	chars = [x for x in string.ascii_lowercase]
	dupes = []
	ref_dict = {}

	for entry in bib.entries:
		firstname = entry["author"].split(",", 1)[0]
		if "-" in firstname:
			firstname = firstname.replace("-", "")
		if " " in firstname:
			firstname = firstname.replace(" ", "")
		year = entry["year"]
		entry["ID"] = "".join([firstname, year])

	with open(output, "w") as bibtex:
		btx.dump(bib, bibtex)
