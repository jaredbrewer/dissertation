from ij import IJ, ImagePlus
from os import path

outputdir = "" # Insert save directory as a string.
imp = IJ.getImage()
title = imp.getTitle()

out = path.join(outputdir, title)
IJ.saveAs(imp, "Tiff", out)
