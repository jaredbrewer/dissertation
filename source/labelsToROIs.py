from ij import IJ, ImagePlus
from ij.gui import Wand, PolygonRoi, Roi
from ij.plugin import Duplicator
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
import os

def labelsToROIs(imp):

  rm = RoiManager.getRoiManager()

  def L2R(imp):
    ip = imp.getProcessor()
    wand = Wand(ip)

    width = imp.getWidth()
    height = imp.getHeight()

    pixel_width = []

    for num in range(0, width - 1):
        pixel_width.append(num)

    pixel_height = []

    for num in range(0, height - 1):
        pixel_height.append(num)

    ip.setColor(0)

    for y_coord in pixel_height:
      for x_coord in pixel_width:
        if ip.getPixel(x_coord, y_coord) > 0:
          wand.autoOutline(x_coord, y_coord)
          if wand.npoints > 0:
            roi = PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.TRACED_ROI)
            roi.setPosition(imp.getCurrentSlice())
            ip.fill(roi)
            rm.addRoi(roi)
    rx = rm.getRoisAsArray()
    rm.reset()
    rois = []
    for r in rx:
      if r not in rois:
        rois.append(r)
        rm.addRoi(r)

    rm.runCommand(imp, "Show All")

  copy_imp = imp.duplicate()
  rm.reset()

  dimensions = copy_imp.getDimensions()
  nChannels = dimensions[2]
  nSlices = dimensions[3]
  nFrames = dimensions[4]

  if (nChannels > 1 and nSlices > 1) or (nChannels > 1 and nFrames > 1) or (nSlices > 1 and nFrames > 1):
    print(imp.getTitle() + " is a hyperstack (multi c , z or t), please prepare a stack (single c, either z-stack or t-stack) from it.")
    pass
  elif nChannels > 1 or nSlices > 1 or nFrames > 1:
    for stack in range(1, copy_imp.getImageStackSize() + 1):
      copy_imp.setPosition(stack)
      L2R(copy_imp)
  else:
    L2R(copy_imp)


