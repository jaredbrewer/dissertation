from ij import IJ, ImagePlus
from ij.gui import GenericDialog
from csv import writer
from os import path

gui = GenericDialog("Surface Plot Measure")
gui.addDirectoryField("File Location: ", "~/Documents")
gui.showDialog()

basedir = str(gui.getNextString())

inputdir = basedir
outputdir = basedir + "csvs/"
if not path.isdir(outputdir):
  os.makedirs(outputdir)

filenames = os.listdir(inputdir)

plot_files = []

for files in filenames:
  if files.startswith("Surface_Plot") and files.endswith(".tif"):
    file = os.path.join(inputdir, files)
    plot_files.append(file)

for f in plot_files:
  imp = IJ.openImage(f)
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

  csv_file = path.join(outputdir, title.split(".", 1)[0] + ".csv")

  with open(csv_file, "a") as csv:
    writer = writer(csv)
    for coord in max_coord:
      writer.writerow(coord)
