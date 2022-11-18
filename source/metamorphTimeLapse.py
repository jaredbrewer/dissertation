from ij import IJ, ImagePlus
from ij.plugin import ZProjector, RGBStackMerge, Concatenator, Duplicator
from ij.gui import GenericDialog, NonBlockingGenericDialog
import os, sys, random, re
from os import path

gui = GenericDialog("Metamorph File Compiler")
gui.addMessage("Version 0.1a, 17 November 2022")
gui.addDirectoryField("Files Directory: ", "~/Documents")
gui.addCheckbox("Multiple Scenes?", True)
gui.addCheckbox("Multiple Time Points?", True)
gui.addCheckbox("Show Merged Image?", False)
gui.addCheckbox("Use Metadata from .nd File?", False)

gui.showDialog()

basedir = str(gui.getNextString())
multiscene = gui.getNextBoolean()
timelapse = gui.getNextBoolean()
opener = gui.getNextBoolean()
nder = gui.getNextBoolean()

outputdir = path.join(basedir, "comps")
if not path.isdir(outputdir):
  os.makedirs(outputdir)

inputs = os.listdir(basedir)
filenames = []
ext = []

for input in inputs:
  if "thumb" not in input and input.lower().endswith(".tif") and not input.startswith(".") and not input.startswith("_"):
    filenames.append(input)
    ext.append("." + input.rsplit(".", 1)[1])

ext = set(ext)

parts = []

for files in filenames:
  part = files.rsplit("_", 3)
  parts.append(part)

names = []
chans = []
scenes = []
times = []

for p in parts:
  names.append(p[0])
  if multiscene and timelapse:
    chans.append(p[1])
    scenes.append(int(p[2].replace("s", "")))
    times.append(int(p[3].split(".")[0].replace("t", "")))
  elif multiscene:
    chans.append(p[1])
    scenes.append(int(p[2].split(".")[0].replace("s", "")))
  elif timelapse:
    chans.append(p[1])
    times.append(int(p[2].split(".")[0].replace("t", "")))
  else:
    chans.append(p[1]).split(".")[0]

names = list(set(names))
chans = list(set(chans))
scenes = list(sorted(set(scenes)))
times = list(sorted(set(times)))

chan_dict = {}
chan_pairs = ""
for chan in chans:
  chan_dict["C" + str(int(chans.index(chan) + 1))] = chan
  chan_pairs += "C" + str(int(chans.index(chan) + 1)) + ": " + chan + "\n"

chan_list = ["C" + str(int(chans.index(chan) + 1)) for chan in chans]
chan_list.append("None")

gui = NonBlockingGenericDialog("Color Matcher")
gui.addMessage(chan_pairs)
gui.addChoice("Red: ", chan_list, "None")
gui.addChoice("Green: ", chan_list, "None")
gui.addChoice("Blue: ", chan_list, "None")
gui.addChoice("Grey: ", chan_list, "None")
gui.addChoice("Cyan: ", chan_list, "None")
gui.addChoice("Magenta: ", chan_list, "None")
gui.addChoice("Yellow: ", chan_list, "None")

gui.showDialog()

red_pick = str(gui.getNextChoice())
green_pick = str(gui.getNextChoice())
blue_pick = str(gui.getNextChoice())
grey_pick = str(gui.getNextChoice())
cyan_pick = str(gui.getNextChoice())
magenta_pick = str(gui.getNextChoice())
yellow_pick = str(gui.getNextChoice())

colors = [red_pick, green_pick, blue_pick, grey_pick, cyan_pick, magenta_pick, yellow_pick]
for key, value in chan_dict.items():
  if key not in colors:
    del chan_dict[key]

chans = chan_dict

def metamorpher(names, chans, scenes, times, ext, colors, timelapse = True, multiscene = True, opener = False, nder = False):

  def titler(nder, name, scene):
    if nder:
      if path.isfile(path.join(basedir, name + ".nd")):
        pos_dict = {}
        with open(path.join(basedir, name + ".nd"), "r") as nd_file:
          for line in nd_file.readlines():
            if line.startswith("\"Stage"):
              pos_dict[line.split(",")[0].strip("\"\n").replace("Stage", "")] = line.split(",")[1].strip("\"\n\r ")
        if str(scene) in pos_dict:
          scene_name = pos_dict.get(str(scene))
          title = "_".join([name, scene_name])
    else:
      title = "_".join([name, str(scene)])
    return title

  def stiller(timelapse, combo_dict, name, chan, times, scene = "scene"):
    stills = []
    if timelapse:
      for time in times:
        still = path.join(basedir, "_".join([name, chan, "s" + str(scene), "t" + str(time) + "".join(ext)]))
        if path.isfile(still):
          still = ImagePlus(still)
          if still.getNSlices() > 1:
            still = ZProjector.run(still, "max")
            stills.append(still)
          else:
            stills.append(still)
      if stills:
        combo_dict[key] = Concatenator.run(stills)
    else:
      still = path.join(basedir, "_".join([name, chan, "s" + str(scene) + "".join(ext)]))
      if path.isfile(still):
        still = ImagePlus(still)
        if still.getNSlices() > 1:
          still = ZProjector.run(still, "max")
          stills.append(still)
        else:
          stills.append(still)
      if stills:
        combo_dict[key] = ImagePlus(still)

  def merge_and_save(combo_dict, colors, title, outputdir, opener):
    color_picker = []
    for color in colors:
      if color != "None":
        color_picker.append(combo_dict[color])
      elif color == "None":
        color_picker.append(None)
    merge = RGBStackMerge.mergeChannels(color_picker, False)
    if opener:
      merge.show()
    out = path.join(outputdir, title)
    del title
    IJ.saveAs(merge, "Tiff", out)
    rand = random.randint(1, 100)
    if rand > 85:
      IJ.run("Collect Garbage", "")
    else:
      pass

  for name in names:
    if timelapse and multiscene:
      for scene in scenes:
        title = titler(nder, name, scene)
        print(title)
        combo_dict = {}
        for key, chan in chans.items():
          stiller(timelapse, combo_dict, name, chan, times, scene)
        if combo_dict:
          merge_and_save(combo_dict, colors, title, outputdir, opener)
    elif timelapse:
      title = name
      combo_dict = {}
      for key, chan in chans.items():
        stiller(timelapse, combo_dict, name, chan, times, scene)
      if combo_dict:
        merge_and_save(combo_dict, colors, title, outputdir, opener)
    elif multiscene:
      # If you have a multiple scene, non-timelapse image (I'm not sure how or why it would have been set up like that?), there are definitely better options out there, but this should work?
      for scene in scenes:
        title = titler(nder, name, scene)
        combo_dict = {}
        for key, chan in chans.items():
          stiller(timelapse, combo_dict, name, chan, times, scene)
        if combo_dict:
          merge_and_save(combo_dict, colors, title, outputdir, opener)
    else:
      title = name
      combo_dict = {}
      num_stacks = []
      for key, chan in chans.items():
        stiller(timelapse, combo_dict, name, chan, times, scene)
      if combo_dict:
        merge_and_save(combo_dict, colors, title, outputdir, opener)

metamorpher(names, chans, scenes, times, ext, colors, timelapse, multiscene, opener, nder)

gui = GenericDialog("All Done!")
gui.showDialog()
