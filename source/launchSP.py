from ij import IJ, ImagePlus

imp = IJ.getImage()

IJ.run(imp, "3D Surface Plot", 
	'''plotType=1 
	colorType=3 
	drawAxes=0 
	drawLines=0 
	drawText=0 
	grid=256 
	drawLegend=0 
	smooth=8.5 
	backgroundColor=000000 
	windowHeight=600 
	windowWidth=720''')
