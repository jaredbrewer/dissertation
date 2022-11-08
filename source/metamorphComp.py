from ij import IJ, ImagePlus
from ij.plugin import ZProjector, RGBStackMerge
from ij.gui import GenericDialog
import os, sys, random, re
from os import path

gui = GenericDialog("File Directory: ")
gui.addDirectoryField("Individual Channels: ", "~/Documents")
gui.showDialog()

basedir = str(gui.getNextString())

outputdir = basedir + "/comps/"
if not path.isdir(outputdir):
  os.makedirs(outputdir)

inputs = os.listdir(basedir)
filenames = []

for input in inputs:
  if input.startswith("MAX_"):
    filenames.append(input)

names = []
ends = []

for files in filenames:
  base = str(files.rsplit("_", 1)[0])
  end = str(files.rsplit("_", 1)[1])
  names.append(base)
  ends.append(end)

names = set(names)
# Sorting vastly improved reliability and facilitates static bindings in the merge.
# This entire logic is in support of arbitrary channel number, even though functionally 7 is the most you can do.
ends = sorted(set(ends))
mip_dict = {}

for name in names:
  for end in ends:
    mip = path.join(basedir, name + "_" + end)
    mip_dict["C" + str(int(ends.index(end) + 1))] = mip
  title = name
  # Order is RGBKCMY, K = Gray
  merge = RGBStackMerge.mergeChannels([ImagePlus(mip_dict['C2']), None, ImagePlus(mip_dict['C3']), None, None, None, ImagePlus(mip_dict['C1'])], False)
  out = path.join(outputdir, title)
  del title
  IJ.saveAs(merge, "Tiff", out)

gui = GenericDialog("All Done!")
gui.showDialog()
