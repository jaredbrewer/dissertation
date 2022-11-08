from ij import IJ, ImagePlus
from ij.io import DirectoryChooser, FileSaver
from ij.plugin import ChannelSplitter, ZProjector
from ij.plugin.frame import RoiManager
from ij.gui import GenericDialog
import os, sys
from os import path

gui = GenericDialog("File Directory: ")
gui.addDirectoryField("Cell Images", "~/Documents")
gui.addChoice("File Extension: ", [".czi", ".lif", ".lsm", ".tif", ".tiff", "Other"], ".tif")
gui.addStringField("Custom File Extension: ", "")
gui.addNumericField("Reference Channel: ", 2, 0)
gui.addNumericField("Measurement Channel: ", 1, 0)
gui.showDialog()

basedir = str(gui.getNextString())

inputdir = basedir
outputdir = basedir + "/Splits/"
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
    if not files.startswith("."):
      file = os.path.join(inputdir, files)
      rep_files.append(file)

for f in rep_files:
  imp = IJ.openImage(f)
  proj = ZProjector.run(imp, "max all")
  imps = ChannelSplitter.split(proj)
  for i in imps:
    title = i.getTitle()
    out = path.join(outputdir, title)
    IJ.saveAs(i, "Tiff", out)
    ImagePlus.close(i)

###

splits = os.listdir(outputdir)

tifs = []

for files in splits:
  if files.endswith(".tif"):
    file = os.path.join(outputdir, files)
    tifs.append(file)

names = []

for tif in tifs:
  base = str(tif.split("-", 1)[1])
  names.append(base)

names = set(names)

rm = RoiManager.getRoiManager()

tifs.sort(reverse = True)

ref = str(int(gui.getNextNumber()))
meas = str(int(gui.getNextNumber()))

for name in names:
  for tif in tifs:
    if tif.startswith(outputdir + "C" + ref) and tif.endswith(name):
      a = IJ.openImage(tif)
      IJ.setAutoThreshold(a, "MinError dark")
      IJ.run(a, "Analyze Particles...", "size=100.00-Infinity clear include add")
      ra = rm.getRoisAsArray()
      rm.reset()
      ImagePlus.close(a)
    if tif.startswith(outputdir + "C" + meas) and tif.endswith(name):
      b = IJ.openImage(tif)
      title = b.getTitle().lstrip(outputdir + "C" + meas + "-")
      out = path.join(outputdir, title)
      b.show()
      for r in ra:
        rm.addRoi(r)
      rm.select(b, 0)
      IJ.run(b, "Clear Outside", "")
      rm.reset()
      IJ.saveAs(b, "Tiff", out)

fin = GenericDialog("All Done!")
fin.showDialog()
