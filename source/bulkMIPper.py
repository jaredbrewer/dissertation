from ij import IJ, ImagePlus
from ij.plugin import ZProjector
from ij.gui import GenericDialog
import os, sys, random
from os import path

gui = GenericDialog("File Directory: ")
gui.addDirectoryField("Z-Stacks Location", "~/Documents")
gui.addChoice("File Extension: ", [".czi", ".lif", ".lsm", ".tif", ".tiff", "Other"], ".tif")
gui.addStringField("Custom File Extension: ", "")
gui.showDialog()

basedir = str(gui.getNextString())

inputdir = basedir
outputdir = basedir + "MIPs/"
if not path.isdir(outputdir):
  os.makedirs(outputdir)

filenames = os.listdir(inputdir)

ext = str(gui.getNextChoice())
cusext = str(gui.getNextString())

extension = ""

if ext != "Other":
  extension = ext
elif ext == "Other":
  extension = cusext

rep_files = []

for files in filenames:
  if files.endswith(extension):
    file = os.path.join(inputdir, files)
    rep_files.append(file)

for f in rep_files:
  imp = IJ.openImage(f)
  if "MAX_" + imp.getTitle().rstrip(extension) + ".tif" not in os.listdir(outputdir):
    proj = ZProjector.run(imp, "max all")
    title = proj.getTitle()
    out = path.join(outputdir, title)
    IJ.saveAs(proj, "Tiff", out)
    imp.close()
    proj.close()
  rand = random.randint(1, 100)
  if rand > 85:
    IJ.run("Collect Garbage", "")
  else:
    pass
