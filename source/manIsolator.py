from ij import IJ, ImagePlus
from ij.plugin import ZProjector, RGBStackMerge, ChannelSplitter
from ij.plugin.frame import RoiManager
from csv import writer

imp = IJ.getImage()
proj = ZProjector.run(imp, "max all")

imps = ChannelSplitter.split(proj)
rm = RoiManager.getRoiManager()

imp.close()

for i in imps:
	# You will need to change the "C#" to whatever channel numbers are appropriate for your use case!
	if i.getTitle().startswith("C2"):
		IJ.getImage()
		i.show()
		IJ.setAutoThreshold(i, "MinError dark")
		IJ.run(i, "Analyze Particles...", "size=100.00-Infinity clear include add slice")
		i.close()
	if i.getTitle().startswith("C1"):
		IJ.getImage()
		i.show()
		rm.select(0)
		IJ.run(i, "Clear Outside", "slice")

