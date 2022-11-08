#!/usr/bin/python3

import sys, argparse


# @click.command("command", )
# @click.option("folder", required = True, help = "Provide the folder you wish to rename. Drag and drop is usually fine.")
# @click.option("--r", help = "Do you want to rename files in the subdirectories as well? If so, provide [r] or [r = True]")

def blindrename(folder, r = False):
	print("Blind Renamer")

def unblind(folder):
	print("Unblinder")

parser = argparse.ArgumentParser(description='Blind rename/unblind a folder of files')
parser.add_argument("function", nargs = 1, action = "store", help = "Which command would you like to run? [blindrename] or [unblind]")
parser.add_argument("folder", nargs = 1, action = "store", help = "Provide the folder you wish to rename. Drag and drop is usually fine.")
parser.add_argument("--r", default = False, required = False, action = "store_true", help = "Do you want to rename files in the subdirectories as well? If so, provide [r] or [r = True]")

args = parser.parse_args()
print(args)
if "blindrename" in args.function:
	blindrename(args.folder, args.r)
if "unblind" in args.function:
	unblind(args.folder)

# if __name__ == '__main__':
# 	try:
# 		sys.argv[3]
# 	except IndexError:
# 		sys.argv.insert(3, False)
# 	globals()[sys.argv[1]](sys.argv[2], sys.argv[3])
# 	if sys.argv[1] == "unblind":
# 		sys.argv.pop(3)
# 	print(sys.argv)
# 	parser = argparse.ArgumentParser(description='Blind rename/unblind a folder of files')
# 	parser.add_argument("command", nargs = 1, choices = ["blindrename", "unblind"], help = "Which command would you like to run? [blindrename] or [unblind]")
# 	parser.add_argument("folder", nargs = 1, help = "Provide the folder you wish to rename. Drag and drop is usually fine.")
# 	parser.add_argument("--r", nargs = 1, default = False, required = False, help = "Do you want to rename files in the subdirectories as well? If so, provide [r] or [r = True]")
#
# 	args = parser.parse_args()
# 	args.command(args.folder)
