from ij import IJ, ImagePlus, ImageStack, WindowManager
from ij.gui import GenericDialog
from ij.measure import Calibration
from ij.plugin import Plugin
from ij.process import ImageProcessor
from os import path

# Defaults:
first = 1
last = 9999
inc = 2
imp = IJ.getImage()

def sliceKeeper(imp, first, last, inc):

  imp = IJ.openImage(imp)
  if not imp:
    IJ.noImage()
  stack = imp.getStack()
  if stack.getSize() is 1:
    IJ.error("Stack Required")
  title = imp.getTitle().split(".")[0]
  cal = imp.getCalibration()
  impdir = str(imp.getOriginalFileInfo().directory)

  def keepSlices(stack, first, last, inc):
    if last > stack.getSize():
      last = stack.getSize()
    newstack = ImageStack(stack.getWidth(), stack.getHeight())
    for slice in range(first, last + 1, inc):
      if slice > stack.getSize():
        break
      ip = stack.getProcessor(slice)
      newstack.addSlice(ip)
    imp = ImagePlus(title + "_kept", newstack)
    imp.setCalibration(cal)
    outputdir = path.join(impdir, title)
    IJ.saveAs(imp, "Tiff", outputdir)

  keepSlices(stack, first, last, inc)

def sliceRemover(imp, first, last, inc):

  imp = IJ.openImage(imp)
  if not imp:
    IJ.noImage()
  stack = imp.getStack()
  if stack.getSize() is 1:
    IJ.error("Stack Required")
  title = imp.getTitle().split(".")[0]
  cal = imp.getCalibration()
  impdir = str(imp.getOriginalFileInfo().directory)

  def removeSlices(stack, first, last, inc):
    if last > stack.getSize():
      last = stack.getSize()
    newstack = stack.duplicate()
    for slice in range(first, last + 1, inc):
      if slice > stack.getSize():
        break
      newstack.deleteSlice(slice)
    imp = ImagePlus(title + "_cut", newstack)
    imp.setCalibration(cal)
    outputdir = path.join(impdir, title)
    IJ.saveAs(imp, "Tiff", outputdir)

  removeSlices(stack, first, last, inc)


