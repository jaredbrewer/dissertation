from ij import IJ, ImagePlus
from ij.gui import GenericDialog
from csv import writer
from os import path

gui = GenericDialog("Surface Plot Measure")
gui.addDirectoryField("Save Location: ", "~/Documents")
gui.showDialog()
save_loc = str(gui.getNextString())

imp = IJ.getImage()
title = imp.getTitle()

IJ.run(imp, "8-bit", "")
IJ.setRawThreshold(imp, 1, 255)
IJ.run(imp, "Make Binary", "method=Default")

width = imp.getWidth()
height = imp.getHeight()

pixel_width = []

for num in range(0, width-1):
  pixel_width.append(num)

pixel_height = []

for num in range(0, height-1):
  pixel_height.append(num)

max_coord = []

for x_coord in pixel_width:
  for y_coord in pixel_height:
    if imp.getPixel(x_coord, y_coord)[0] == 255:
      max_coord.append([imp.getTitle(), x_coord, y_coord])

csv_file = path.join(save_loc, title.split(".", 1)[0] + ".csv")

with open(csv_file, "a") as csv:
  writer = writer(csv)
  for coord in max_coord:
    writer.writerow(coord)
