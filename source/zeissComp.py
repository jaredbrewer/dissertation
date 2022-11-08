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
chans = []

for files in filenames:
  base = str(files.split("-", 1)[1])
  chan = str(files.split("-", 1)[0])
  names.append(base)
  chans.append(chan)

names = set(names)
chans = sorted(set(chans))
mip_dict = {}

for name in names:
  for chan in chans:
    # Watch the separator - it has played me many times. They all need to be the same
    mip = path.join(basedir, chan + "-" + name)
    mip_dict["C" + str(int(chans.index(chan) + 1))] = mip
  title = name
  merge = RGBStackMerge.mergeChannels([ImagePlus(mip_dict['C3']), None, ImagePlus(mip_dict['C2']), None, ImagePlus(mip_dict['C4']), None, ImagePlus(mip_dict['C1'])], False)
  out = path.join(outputdir, title)
  del title
  IJ.saveAs(merge, "Tiff", out)

gui = GenericDialog("All Done!")
gui.showDialog()
