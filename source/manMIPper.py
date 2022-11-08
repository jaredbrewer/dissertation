from ij import IJ, ImagePlus
from ij.plugin import ZProjector
from ij.gui import NonBlockingGenericDialog
import os, sys, random
from os import path

gui = NonBlockingGenericDialog("MIP Ranges")

gui.addNumericField("First Slice: ", 0, 0)
gui.addNumericField("Last Slice: ", 0, 0)

gui.showDialog()

if gui.wasCanceled():
  IJ.error("Dialog cancelled!")
  exit()
# Everything below this ~should~ be in an else statement, but if it isn't broke, don't fix it.

bot = gui.getNextNumber()
top = gui.getNextNumber()

imp = IJ.getImage()

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

outputdir = str(imp.getOriginalFileInfo().directory)

proj = ZProjector.run(imp, "max", int(bot), int(top))
# This line is to show your MIP - if you want to see it, remove the #
# proj.show()

title = proj.getTitle()
out = path.join(outputdir, title)
IJ.saveAs(proj, "Tiff", out)

#Comment out these lines if you'd like to see your images before you close them.
imp.close()
proj.close()

rand = random.randint(1, 100)

if rand > 85:
  IJ.run("Collect Garbage", "")
else:
  pass
