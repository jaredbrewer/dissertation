from ij import IJ, ImagePlus
from ij.measure import ResultsTable, Measurements
from ij.plugin.frame import RoiManager
from ij.process import AutoThresholder
from ij.plugin import ZProjector
from ij.plugin.filter import ParticleAnalyzer
from ij.gui import GenericDialog, NonBlockingGenericDialog
import os, random, csv, sys, re
from os import path

def burden(directory, chan, min_threshold, ext, screen_threshold = "Otsu dark", proj_save = False, proj_show = False, imp_show = False, fish_channel = None, outline_threshold = "Triangle dark", brightfield = False, subset = 0, man_psize = 100000):

# directory is a string with path to the files you wish to analyze.
# chan is an integer of the channel # with the bacteria.
# min_threshold is an integer of the minimum threshold value (empirically determined by checking a few images)
# ext is the file extension (".czi", ".lsm", ".tif" are some examples) as a string
# screen_threshold is the threshold algorithm to find the yolk as a string. I found that "Otsu dark" works well, but others may be better in other circumstances.
# proj_save is a boolean that determines whether or not to save the MIPs as you go.
# proj_show is a boolean to determine whether to show the MIPs as they are processed. A bit memory intensive, but valuable for sanity checking.
# imp_show is a boolean for whether or not to show the source image (not recommended except with small datasets)
# fish_channel is an integer for the channel containing something we could use as the basis for an outline of the fish.
# outline_threshold is a string with the thresholding algorithm to find the whole fish (you want something generous, but probably not MinError or similar).
# brightfield is a boolean telling whether or not the selected fish_channel is the brightfield channel. This can be useful and very convenient, but requires some additional processing to be useful.
# subset is to run the function over a subset of the total number of images. Provide an integer value less than the total number of images in the directory. If subset is selected, the projections will automatically show.
# man_psize is a manual particle size value for the "Analyze Particles..." function. Default has worked well for me, but can be changed if needed.

  filenames = os.walk(str(directory))
  bfiles = []

  if ext.startswith("."):
    # Rigid binding of extension per se.
    ext = ext.replace(".", "")
  # This logic is a little over-complex for the task, but essentially I want to pull out the extension and then use it to match the end of the line. Why do it in two steps? It makes the code more flexible for ".tif" vs ".TIFF" and could allow future use of the mext in downstream applications if needed.
  for dirpath, subdir, files in filenames:
    for file in files:
      mext = re.search("\." + ext + "\S*$", file, re.IGNORECASE)
      if mext:
        if file.endswith(mext.group()):
          bfiles.append(path.join(dirpath, file))

  # Make sure all our measurements are set properly.
  IJ.run("Set Measurements...",
  "area mean min limit display redirect=None decimal=3")

  rm = RoiManager.getRoiManager()
  rt = ResultsTable.getResultsTable()
  data = []

  if subset > len(bfiles):
    wrong_subset = NonBlockingGenericDialog("Subset Error")
    wrong_subset.addMessage("The subset is greater than the number of files - running on whole directory.")
    wrong_subset.showDialog()
    subset = 0
  if subset:
    bfiles = bfiles[0:int(subset)]
    proj_show = True
  else:
    subset = 0

  valid_thresholds = [str(i + " dark") for i in AutoThresholder.getMethods()]
  valid_thresholds.append("Try all")
  valid_thresholds.append("None")

  if screen_threshold == "Try all" or outline_threshold == "Try all":
    if screen_threshold == "Try all":
      for f in bfiles:
        imp = IJ.openImage(f)
        proj = ZProjector.run(imp, "max all")
        proj.setC(int(chan))
        IJ.run(proj, "Auto Threshold", "method=[Try all] dark")
    if outline_threshold == "Try all":
      for f in bfiles:
        imp = IJ.openImage(f)
        proj = ZProjector.run(imp, "max all")
        proj.setC(int(fish_channel))
        if brightfield:
          IJ.run(proj, "Find Edges", "slice");
        IJ.run(proj, "Auto Threshold", "method=[Try all] dark")
  else:
    for i, f in enumerate(bfiles):
      imp = IJ.openImage(f)
      if imp_show:
        imp.show()
      if imp.getDimensions()[3] > 1:
        proj = ZProjector.run(imp, "max all")
      else:
        proj = imp.duplicate()
      if proj_show:
        proj.show()
      proj.setC(int(chan))
      if screen_threshold not in valid_thresholds:
        screen_threshold = "Otsu dark"

      if screen_threshold != "None":
        IJ.setAutoThreshold(proj, str(screen_threshold))

        IJ.run(proj, "Analyze Particles...", "size="+str(man_psize)+"-Infinity clear include add slice")

        if rm.getCount() > 1: # This is the dumbest source of a downstream error I have ever encountered. rm.getCount is an attribute and rm.getCount() is a method = two different outcomes, but no evaluation error.
          rm.runCommand(proj, "Select All")
          rm.runCommand(proj, "Combine")
          rm.runCommand(proj, "Delete")
          rm.runCommand(proj, "Add")

        rm.select(proj, 0)

        IJ.run(proj, "Clear", "slice")
        IJ.run(proj, "Select None", "")
        IJ.run(proj, "Remove Overlay", "")
        rm.reset()

      # Build some logic for using the outline of the fish if provided. Brightfield images probably need even more work, but this will *probably* work for fluorescent channels.
      if fish_channel:
        proj.setC(int(fish_channel))
        if brightfield:
          IJ.run(proj, "Find Edges", "slice");
        if outline_threshold not in valid_thresholds:
          outline_threshold = "Triangle dark"
        IJ.setAutoThreshold(proj, outline_threshold)
        IJ.run(proj, "Analyze Particles...", "size="+str(man_psize)+"-Infinity clear include add slice")
        proj.setC(int(chan))
        ra = rm.getRoisAsArray()
        for r in ra:
          rm.addRoi(r)
        if rm.getCount() > 1:
          rm.runCommand(proj, "Select All")
          rm.runCommand(proj, "Combine")
          rm.runCommand(proj, "Delete")
          rm.runCommand(proj, "Add")

        rm.select(proj, 0)

      else:
        IJ.run(proj, "Select All", "")

      # We want this to work for any bit depth image, but due to limitations of the underlying plugin, may only work for 16-bit and lower images. Something weird could happen with 24-bit RGB images, but I can't stop people from doing weird things.

      if proj.getBitDepth() == 24: # This should only be color images - seems a little risky.
        IJ.run(proj, "8-bit", "")
      elif proj.getBitDepth() > 16: # The thresholds don't work for 32-bit images and FIJI doesn't really understand 10 or 12 bit images.
        IJ.run(proj, "16-bit", "")

      IJ.setRawThreshold(proj, int(min_threshold), 2**int(proj.getBitDepth()) - 1)
      IJ.run(proj, "Measure", "")
      row = rt.getRowAsString(0).split("\t")
      row.insert(0, path.join(imp.getOriginalFileInfo().directory, imp.getOriginalFileInfo().fileName))
      if i == 0:
        headings = rt.getColumnHeadings().split("\t")
        headings.insert(0, "file_path")
      data.append(row)

      if not imp_show:
        imp.close()
      if not proj_show:
        proj.close()

      if proj_save:
        out = path.join(directory, proj.getTitle())
        IJ.saveAs(proj, "Tiff", out)

      rand = random.randint(1, 100)
      if rand > 85:
        IJ.run("Collect Garbage", "")
      else:
        pass

    # Write out the results to a CSV file for further processing (probably in R, but also fine for other options.)

    with open(path.join(directory, "burden.csv"), "w") as csvfile:
      writer = csv.writer(csvfile)
      writer.writerow(headings)
      for datum in data:
        writer.writerow(datum)

  gui = GenericDialog("Finished")
  gui.addMessage("All Done! Output is located in " + directory)
  gui.showDialog()

def burden_gui():

  values = [str(x + " dark") for x in AutoThresholder.getMethods()]
  values.append("Try all")
  values.append("None")

  # Build the GUI elements.
  gui = NonBlockingGenericDialog("Bacterial Burden Measurement")
  gui.addDirectoryField("Image Folder: ", " ~/Documents")
  gui.addChoice("File Extension: ", [".czi", ".lif", ".lsm", ".tif", "Other"], ".tif")
  gui.addStringField("Custom File Extension: ", ".jpeg")
  gui.addNumericField("Bacterial Fluorescence Channel #: ", 1)
  gui.addNumericField("Minimum Threshold: ", 200)
  gui.addMessage("'Try all' requires a subset to be selected - the fewer the better. Default = 5")
  gui.addChoice("Background Removal Threshold: ", values, "Otsu dark")
  gui.addCheckbox("Save MIPs?", False)
  gui.addCheckbox("Show MIPs?", False)
  gui.addCheckbox("Show Images?", False)
  gui.addCheckbox("Use Fish Outline?", False)
  gui.addToSameRow()
  gui.addNumericField("Channel # for Fish Outline: ", 2)
  gui.addChoice("Fish Threshold: ", values, "Triangle dark")
  gui.addCheckbox("Brightfield Image?", False)
  gui.addCheckbox("Only do a subset?", False)
  gui.addToSameRow()
  gui.addNumericField("If yes, how many images to try: ", 0)
  gui.addNumericField("Manual Particle Size: ", 100000)
  gui.showDialog()
  if gui.wasCanceled():
    canceled = GenericDialog("Canceled")
    canceled.showDialog()
    sys.exit(0)

  # Parse all of the GUI elements into variables that can be plugged into the burden function.
  directory = str(gui.getNextString())
  ext = str(gui.getNextChoice())
  if ext == "Other":
    ext = str(gui.getNextString()).lower()
  else:
    dump = gui.getNextString()

  chan = int(gui.getNextNumber())
  min_threshold = int(gui.getNextNumber())
  screen_threshold = str(gui.getNextChoice())
  proj_save = gui.getNextBoolean()
  proj_show = gui.getNextBoolean()
  imp_show = gui.getNextBoolean()
  do_fish = gui.getNextBoolean()
  if do_fish:
    fish_channel = int(gui.getNextNumber())
    outline_threshold = str(gui.getNextChoice())
    brightfield = gui.getNextBoolean()
  else:
    dump = gui.getNextNumber()
    dump = gui.getNextChoice()
    dump = gui.getNextBoolean()
    fish_channel = None
    outline_threshold = "Triangle dark"
    brightfield = False
  subsetter = gui.getNextBoolean()
  if subsetter:
    subset = int(gui.getNextNumber())
  else:
    dump = gui.getNextNumber() # Another extraordinarily stupid source of error - every entry /must/ be accessed or you're going to frameshift - was feeding 0 into man_psize when subset = False and getting a correspondingly bad result (all points included).
    subset = 0
  if screen_threshold == "Try all":
    subsetter = True
    if not subset:
      subset = 5
  if outline_threshold == "Try all":
    subsetter = True
    if not subset:
      subset = 5
  man_psize = int(gui.getNextNumber())

  burden(directory = directory,
    chan = chan,
    min_threshold = min_threshold,
    ext = ext,
    screen_threshold = screen_threshold,
    proj_save = proj_save,
    proj_show = proj_show,
    imp_show = imp_show,
    fish_channel = fish_channel,
    outline_threshold = outline_threshold,
    brightfield = brightfield,
    subset = subset,
    man_psize = man_psize)

burden_gui()
