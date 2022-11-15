import re

with open("intro.tex", "r") as file:
	text = file.read()

cites = re.findall(r"\\citep{.*?}", text)

for cite in cites:
	cite_strip = cite.replace("\citep{", "").replace("}", "").replace(" ", "")
	cite_split = cite_strip.split(",")
	cite_split.sort(key = lambda x: int("".join(filter(str.isdigit, x))))
	cite_str = "".join(["\citep{", ", ".join(cite_split), "}"])
	cite_dict[cite] = cite_str

for key, value in cite_dict.items():
	text = text.replace(key, value)

with open("new_intro.tex", "w") as output:
	output.write(text)
