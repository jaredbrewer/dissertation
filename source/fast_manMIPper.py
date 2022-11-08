from ij import IJ, ImagePlus
from ij.plugin import ZProjector
from ij.gui import NonBlockingGenericDialog, GenericDialog
import os, sys, random
from os import path

# This is a low-overhead version of the manMIPper that loads and displays images one at a time and only if needed.

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
  extension = ext.lower()
elif ext == "Other":
  extension = cusext.lower()

rep_files = []

for files in filenames:
  if files.lower().endswith(extension) and "_thumb_" not in files.lower():
    file = os.path.join(inputdir, files)
    rep_files.append(file)

for f in rep_files:
  mips = [mip.lower() for mip in os.listdir(outputdir)]
  mip_title = str("MAX_" + os.path.basename(f)).lower()
  if mip_title.replace(extension, ".tif") not in mips:
    imp = IJ.openImage(f)
    imp.show()

    gui = NonBlockingGenericDialog("MIP Ranges")
    gui.addNumericField("First Slice: ", 0, 0)
    gui.addNumericField("Last Slice: ", 0, 0)
    gui.showDialog()

    if gui.wasCanceled():
      IJ.error("Dialog cancelled!")
      exit()

    bot = gui.getNextNumber()
    top = gui.getNextNumber()

    if bot is 0 and top is 0:
      bot = 1
      top = int(imp.getImageStackSize())
    if int(bot) > int(top):
      IJ.error("Bottom is greater than top, invalid combination.")
      exit()
    if int(top) > int(imp.getImageStackSize()):
      IJ.error("Top is greater than the dimensions of the image - try again")
      exit()
    if int(bot) < 1:
      IJ.error("Out of range - bottom number must be greater than or equal to 1")
      exit()

    proj = ZProjector.run(imp, "max", int(bot), int(top))
    title = proj.getTitle()
    out = path.join(outputdir, title)
    IJ.saveAs(proj, "Tiff", out)
    imp.close()
    proj.close()
  else:
    pass
  rand = random.randint(1, 100)
  if rand > 85:
    IJ.run("Collect Garbage", "")
  else:
    pass
