from ij import IJ, ImagePlus
from ij.gui import Wand, PolygonRoi, Roi
from ij.plugin import Duplicator, Concatenator
from ij.plugin.frame import RoiManager
from ij.process import ImageProcessor
import os, re

def ROIsTolabels(rois, imp_path = None):

  rm = RoiManager.getRoiManager()
  rm.reset()
  rm.runCommand("Open", rois)
  rc = rm.getRoisAsArray()

  if not imp_path:
    pxm = []
    pym = []
    lrs = []
    for r in rc:
      lr = list(r)
      for l in lr:
        lrs.append(str(l).lstrip("java.awt.Point"))
    for lr in lrs:
      px = re.findall("x=\d+", lr)
      for p in px:
        pxm.append(int(p.lstrip("x=")))
      py = re.findall("y=\d+", lr)
      for p in py:
        pym.append(int(p.lstrip("y=")))

    imp_title = os.path.splitext(os.path.basename(rois))[0]
    imp = IJ.createImage(imp_title, "8-bit black", max(pxm), max(pym), 1)
  elif imp_path:
    imp = IJ.openImage(imp_path)

  def R2L(imp):
    label_imp = imp.duplicate()
    IJ.run(label_imp, "Grays", "")
    IJ.run(label_imp, "16-bit", "")
    rc = rm.getRoisAsArray()

    stackN = label_imp.getImageStackSize()
    isStack = bool()
    if stackN > 1:
      isStack = True

    if int(len(rc)) in range(1, 255):
      IJ.run(label_imp, "8-bit", "")
    elif int(len(rc)) in range(256, 65535):
      IJ.run(label_imp, "16-bit", "")
    else:
      IJ.run(label_imp, "32-bit", "")

    for stack in range(1, stackN + 1):
      label_imp.getStack().getProcessor(stack).setValue(0)
      label_imp.getStack().getProcessor(stack).fill()

    for i in range(1, len(rc)):
      if isStack:
        label_imp.setPosition(rc[i].getPosition())
      label_imp.getProcessor().setValue(i + 1)
      label_imp.getProcessor().fill(rc[i])

    label_imp.setTitle("ROIsTolabels_" + imp.getTitle())
    imp_dir = os.path.dirname(rois)
    limp_title = label_imp.getTitle()
    IJ.saveAs(label_imp, "Tiff", os.path.join(imp_dir, limp_title))

  rm.reset()

  dimensions = imp.getDimensions()
  nChannels = dimensions[2]
  nSlices = dimensions[3]
  nFrames = dimensions[4]

  if (nChannels > 1 and nSlices > 1) or (nChannels > 1 and nFrames > 1) or (nSlices > 1 and nFrames > 1):
    print(imp.getTitle() + " is a hyperstack (multi c , z or t), please prepare a stack (single c, either z-stack or t-stack) from it.")
    pass
  elif nChannels > 1 or nSlices > 1 or nFrames > 1:
    for stack in range(1, imp.getImageStackSize() + 1):
      imp.setPosition(stack)
      R2L(imp)
  else:
    R2L(imp)

