#!/usr/bin/python3

import os, csv, string, random, sys, argparse
from os import path

def blindrename(folder, r = False):

  # Root = bad for so many reasons.
  if os.geteuid() == 0:
    exit(f"For safety reasons, you must NOT be root when running renamer.py.  Please become a non-root user, make sure that user has permissions to write to all files in {folder}, and try again.")

  # Needs a folder, if none is provided, exit.
  if not folder:
    exit("You must supply a folder name. I will randomly rename all files in that folder, and create a key file in .CSV format in that folder.")

  # The folder needs to already exist.
  if not path.isdir(folder):
    exit(f"{folder} is not a valid, existing folder.  Try again - typo maybe?")

  # And it presumably needs to have some number of non-hidden files in it. Seems harmless if it is empty, but not expected behavior.

  contents = os.listdir(folder)
  if len(contents) == 0:
    exit("Looks like that folder is empty?")

  # Avoid double renaming.
  if path.exists(os.path.join(folder, "keyfile.csv")):
    exit("Keyfile already exists. Have you already randomized that folder?")

  csvfile = open(path.join(folder, "keyfile.csv"), "w")
  writer = csv.writer(csvfile)
  writer.writerow(["original", "file.path"])

  chars = str(string.ascii_letters + string.digits)

  print("Renaming: ")

  for dirpath, subdirs, files in os.walk(folder):
    subdirs[:] = [d for d in subdirs if not d.startswith(".")]
    files = [f for f in files if not f.startswith(".")]
    if r:
      old_names = [path.join(dirpath, f) for f in files]
    else:
      old_names = [path.join(folder, f) for f in files if dirpath is folder]
    for old_name in old_names:
      if path.isfile(old_name):
        if path.basename(old_name) == "keyfile.csv":
          pass
        else:
          base = path.basename(old_name).rsplit(os.extsep, 1)
          # This ensures that there is a 1 index for files with no extension
          base.append("")
          cloaked_name = "".join(random.choices(chars, k = 5))
          if base[1]:
            new_name = path.join(dirpath, cloaked_name + "." + base[1])
          if not base[1]:
            new_name = path.join(dirpath, cloaked_name)
          print(old_name, " to ", new_name)
          writer.writerow([old_name, new_name])
          os.rename(old_name, new_name)

  print("Finished!")

def unblind(folder):

  if os.geteuid() == 0:
    exit(f"For safety reasons, you must NOT be root when running renamer.py.  Please become a non-root user, make sure that user has permissions to write to all files in {folder}, and try again.")

  # Needs a folder, if none is provided, die
  if not folder:
    exit("You must supply a folder name. I will randomly rename all files in that folder, and create a key file in .CSV format in that folder.")

  # The folder needs to already exist.
  if not path.isdir(folder):
    exit(f"{folder} is not a valid, existing folder.  Try again - typo maybe?")

  # And it presumably needs to have some number of non-hidden files in it. Seems harmless if it is empty, but not expected behavior.

  contents = os.listdir(folder)
  if ".DS_Store" in contents:
    contents.remove(".DS_Store")
  if len(contents) == 0:
    exit("Looks like that folder is empty?")

  # Needs a keyfile to function.
  if not path.exists(os.path.join(folder, "keyfile.csv")):
    exit("No keyfile?")

  csvfile = open(path.join(folder, "keyfile.csv"), "r")
  reader = csv.reader(csvfile)

  print("Renaming: ")

  for file in reader:
    masked = file[1]
    original = file[0]
    if path.exists(masked):
      print(masked + " to " + original)
      os.rename(masked, original)

  print("Finished!")

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Blind rename/unblind a folder of files')
  parser.add_argument("function", nargs = 1, action = "store", help = "Which command would you like to run? [blindrename] or [unblind]")
  parser.add_argument("folder", nargs = 1, action = "store", help = "Provide the folder you wish to rename. Drag and drop is usually fine.")
  parser.add_argument("--r", default = False, required = False, action = "store_true", help = "Do you want to rename files in the subdirectories as well? If so, provide [--r]")

  args = parser.parse_args()
  print(args)
  if "blindrename" in args.function:
    blindrename(args.folder[0], args.r)
  if "unblind" in args.function:
    unblind(args.folder[0])
