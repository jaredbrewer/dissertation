from ij import IJ, ImagePlus
from ij.gui import Roi
from ij.plugin import Duplicator
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor, ImageStatistics
import math, os

def ROIsToMap(imp_path, rois, column_name):

  rm = RoiManager.getRoiManager()
  rm.reset()
  rm.runCommand("Open", rois)

  label = "Measure"
  choices = ["area", "angle", "anglevert", "ar", "circumference", "major", "minor", "mean",
  "median", "mode", "min", "max", "perimeter", "xcenterofmass", "ycenterofmass"]

  if str(column_name).lower() not in choices:
    print("Invalid column_name, try again")
    exit()

  def R2M(imp):

    imp2 = imp.duplicate()
    saveloc = str(imp.getOriginalFileInfo().directory)
    title = imp.getTitle()

    if imp.getBitDepth() == 32:
      IJ.run(imp2, "16-bit", "")

    stackN = imp2.getImageStackSize()
    isStack = bool()

    if stackN > 1:
      isStack = True

    for stack in range(1, stackN + 1):
      imp2.getStack().getProcessor(stack).setValue(0)
      imp2.getStack().getProcessor(stack).fill()

    IJ.run(imp2, "32-bit", "")

    rc = rm.getRoisAsArray()

    filling_value = 0

    for i in range(0, len(rc)):
      if isStack:
        imp.setPosition(rc[i].getPosition())
        imp2.setPosition(rc[i].getPosition())
      ip = imp.getProcessor()
      ip.setRoi(rc[i])

      ip2 = imp2.getProcessor()
      ip2.setRoi(rc[i])

      ip_stats = ip.getStatistics()

      if column_name.lower() == "area":
        filling_value = ip_stats.area
      elif column_name.lower() == "angle":
        filling_value = ip_stats.angle
      elif column_name.lower() == "anglevert":
        filling_value = ip_stats.angle - 90
      elif column_name.lower() == "ar":
        filling_value = ip_stats.major / ip_stats.minor
      elif column_name.lower() == "circumference":
        filling_value = 4 * math.pi * ip_stats.area / math.pow(rx[i].getLength(), 2)
      elif column_name.lower() == "major":
        filling_value = ip_stats.major
      elif column_name.lower() == "minor":
        filling_value = ip_stats.minor
      elif column_name.lower() == "mean":
        filling_value = ip_stats.mean
      elif column_name.lower() == "median":
        filling_value = ip_stats.median
      elif column_name.lower() == "mode":
        filling_value = ip_stats.mode
      elif column_name.lower() == "min":
        filling_value = ip_stats.min
      elif column_name.lower() == "max":
        filling_value = ip_stats.max
      elif column_name.lower() == "perimeter":
        filling_value = rc[i].getLength()
#       elif column_name.lower() = "pattern":
#         roi_name = rc[i].getName()
#         r = Pattern.compile(pattern)
#         m = r.matcher(roi_name)
#         if m.find():
#           group = m.group(1)
#           try:
#             filling_value = Float.parseFloat(group)
      elif column_name.lower() == "xcenterofmass":
        filling_value = ip_stats.xCenterOfMass
      elif column_name.lower() == "ycenterofmass":
        filling_value = ip_stats.yCenterOfMass

      ip2.setValue(filling_value)
      ip2.fill(rc[i])
      imp2.setProcessor(ip2)

    imp2.setTitle(column_name.lower() + "_" + title)
    imp2.show()
    imp2.setRoi(0, 0, imp2.getWidth(), imp2.getHeight())
    IJ.saveAs(imp2, "Tiff", os.path.join(saveloc, imp2.getTitle()))

  R2M(imp)
