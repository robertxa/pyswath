######!/usr/bin/env python
# -*- coding: utf-8 -*-

# Do divisions with Reals, not with integers
# Must be at the beginning of the file
from __future__ import division

# To do :
#	- test other projections, I have not tested it
#   - mask the frequency outside of min/max
#   - remove holes from the dem... TO BE TESTED
#   - remove holes from the min file
# End To Do

###### CHANGES
#    -2015/08/07 : Correct the stepping if two types of profiles (with or without multipoints) are asked in the same run
#                  Add options to choose the size of the graphic
#                  Add options to force the density bar color scale
#    -2015/08/06 : Save a shapefile with the points of the profil
#    -2015/08/04 : Multipoint profils
#    -2015/07/30 : correct bug in function findline to take in account if B>A
#	 -2015/07/30 : remove NaN from graphs
#    -2015/09/15 : Add swath on synthetic DEM
#
######  End CHANGES


############################################################################

## Import Python modules
## I have problems to install rasterio : it does not find gdal libraries... from kingchaos
##modulesNames = ['sys', 'math', 'os', 'utm', 'warnings', 'rasterstats', 'shapely', 'copy', 'time', 'rasterio']
modulesNames = ['sys', 'os', 'warnings']
for module in modulesNames:
	try:
		# because we want to import using a variable, do it this way
		module_obj = __import__(module)
		# create a global object containging our module
		globals()[module] = module_obj
	except ImportError:
		sys.exit(u"ERROR : Module " + module + " not present. \n\n Please, install it \
			      \n\n Edit the source code for more information")
from os import path, access, R_OK, mkdir         # W_OK for write permission.
try:
	import numpy as np                               # need version 1.7 or higher
except ImportError:
	sys.exit(u"ERROR : Module Numpy not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	from osgeo import gdal, gdalnumeric, ogr, osr    # For GIS operations
	from osgeo.gdalconst import *
except ImportError:
	sys.exit(u"ERROR : Module osgeo/gdal not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")

from raster_tools import *
from profiles import *
from plotgraph import *

############################################################################


def checkfiles(rasterfnme, A, B, xsteps, boxwidths, shpbox, title, 
				Coord, synthetic,
				multipoints, remNoData):	
	"""
	
	INPUTS
		- rasterfnme, 
		- A, B, 
		- xsteps, 
		- boxwidths, 
		- shpbox, 
		- title, 
		- Coord, 
		- synthetic,
		- multipoints, 
		- remNoData
	OUTPUTS
		-
		
	"""
	# Check if the DEM exists
	if not path.isfile(rasterfnme) or not access(rasterfnme, R_OK):
		# if not, exist
		raise NameError(u'ERROR: raster {FileNa} dem does not exist'.format(FileNa = rasterfnme))
		
	# Change list to numpy array
	a = np.array(A)
	b = np.array(B)
	a_utm = a
	b_utm = b
	#xstepsa = np.array(xsteps)
	#boxwidthsa = np.array(boxwidths)
	
	test_N = 1
	projdone = False
	
	# Check if the number of points A are the same of point B
	if len(xsteps) != 1 or len(boxwidths) != 1:
		if a.shape != b.shape \
		   or (a.shape[0] != len(xsteps) and len(xsteps) != 1) \
		   or (a.shape[0] != len(multipoints) and len(multipoints) != 1) \
		   or (a.shape[0] != len(boxwidths) and len(boxwidths) != 1) \
		   or (b.shape[0] != len(xsteps) and len(xsteps) != 1) \
		   or (b.shape[0] != len(boxwidths) and len(boxwidths) != 1) \
		   or (b.shape[0] != len(multipoints) and len(multipoints) != 1) \
		   or (len(xsteps) != len(boxwidths) and len(xsteps) != 1 and len(boxwidths) != 1 and len(multipoints) != 1):
			print(u'A =', A)
			print(u'B =', B)
			print(u'xsteps = ', xsteps)
			print(u'boxwidths = ', boxwidths)
			print(u'multipoints', multipoints)
			# If not, exist...
			raise NameError(u"ERROR : the number of points A and points B are different. \n" 
				     u"        Please, check your entries A, B, xsteps, boxwidth and multipoints! \n \n")
	
	# Change the name of the shapefile
	if shpbox[-4:] == '.shp':
		shpbox = 'SHP_' + title + '/' + shpbox[0:-4] + '_' + title + '.shp'
	else:
		shpbox = 'SHP_' + title + '/' + shpbox + '_' + title + '.shp'
	# Erase all the preexisting shapefiles of the project
	#if path.exists("SHP_" + title) == True and path.isfile(shpbox[0:-4] + '1.shp'):
	if path.exists(u"SHP_" + title) == True:
		# erase the preexisting shapefiles
		files = os.listdir("SHP_" + title)
		for f in files:
			full_path = os.path.join("SHP_" + title,f)
			if not os.path.isdir(full_path):
				os.remove(full_path)
	if path.exists("SHP_" + title) == False:
		os.mkdir("SHP_" + title) 
	if path.exists("TMP") == False:
		os.mkdir("TMP")
	if path.exists("TMP"):
		# erase the preexisting shapefile
		files = os.listdir("TMP")
		for f in files:
			full_path = os.path.join("TMP",f)
			if not os.path.isdir(full_path):
				os.remove(full_path) 
	if path.exists("Graphs") == False:
		os.mkdir("Graphs")
	if path.exists("Data") == False:
		os.mkdir("Data")
	
	dst_filename = rasterfnme
	
	print('  ')
	# Check if the xstep is greater than the distance between 2 pixels
	# Load the DEM
	print(u'Reading the DEM...')
	if remNoData:
		# Check if there are no data values
		srcband = rasterdata.GetRasterBand(1)
		#print srcband.GetNoDataValue()
		###### with rasterio if working on your machine
		#with rasterio.drivers():
		#	with rasterio.open(rasterfnme) as src:
		#		
		#		kwargs = src.meta
		#		#print kwargs
		#		#kwargs['nodata'] = nodatav
		#		kwargs.update(nodata = nodatav)
		#		
		#		src = None
		#		#with rasterio.open(rasterfnme, 'w', **kwargs) as dst:
		#		#	dst.write()
		#		
		#		#print src.meta['nodata']
		#		#print src.nodata
		###### with gdal... always working, by tricky... Pb today !!!	
		#if srcband.GetNoDataValue() == None:
		#	rasterdata = None
		#	with rasterio.drivers():
		#		with rasterio.open(rasterfnme) as src:
		#			kwargs = src.meta
		#			kwargs['nodata'] = nodatav
		#	#rasterdata = gdal.Open(rasterfnme, GA_Update)
		#	#srcband = rasterdata.GetRasterBand(1)
		#	#srcband.SetNoDataValue(nodatav)
		#	#rasterdata = None
		#	rasterdata = gdal.Open(rasterfnme, GA_ReadOnly)
		# Get the caracteristics of the DEM
		dst_fnme = rasterfnme[0:-4] + '_fill.tif'
		if path.isfile(dst_fnme):
			print(u'file %s already exists...' % dst_fnme)
			rasterdata = gdal.Open(dst_fnme, GA_ReadOnly)
		else:
			rasterdataor = gdal.Open(rasterfnme, GA_ReadOnly)
			rasterdata = removeNoData(rasterdataor, dst_fnme)
	else:
		rasterdata = gdal.Open(rasterfnme, GA_ReadOnly)
	
	geotransform = rasterdata.GetGeoTransform()
	# Check if xteps is compatible with the DEM resolution
	for iii in range (0, len(xsteps)):
		if min(geotransform[1], geotransform[5]) >= xsteps[iii]:
			warnings.warnings(u'Xstep smaller that the pixel size of the DEM ! \n \
	    		              Xstep is set to 3 times the pixel size')
			xsteps[iii] = max(geotransform[1], geotransform[5]) * 3
	# xtep is now OK for processing
	
	# Check if A and B are in the DEM region or not
	# Calcul the extent of the DEM
	(ulx, uly) = (geotransform[0], geotransform[3])
	(lrx, lry) = (geotransform[0] + geotransform[1] * rasterdata.RasterXSize,
		          geotransform[3] + geotransform[5] * rasterdata.RasterYSize)	
	for i in range (0, a.shape[0]):
		if  not (ulx < A[i][0]< lrx):
			raise NameError(u"ERROR : A%s is outside the DEM border %s" % (A[i], (ulx, lrx)))
		if  not (lry < A[i][1]< uly):
			raise NameError(u"ERROR : A%s is outside the DEM border %s" % (A[i], (lry, uly)))
		if  not (ulx < B[i][0]< lrx):
			raise NameError(u"ERROR : B%s is outside the DEM border %s" % (A[i], (ulx, lrx)))
		if  not (lry < B[i][1]< uly):
			raise NameError(u"ERROR : B%s is outside the DEM border %s" % (A[i], (lry, uly)))
	# Check if the raster is in UTM
	# get projection of the raster source	
	src_proj = rasterdata.GetProjection()
	srs = osr.SpatialReference(wkt = src_proj)
	# To find the projection : srs.GetAttrValue('geogcs')
	if not synthetic and ('UTM' not in src_proj or (Coord[0:3] != 'utm' and Coord[0:3] != 'UTM')):
		projdone = True
		src_geotrans = geotransform
		print(u'The Dem is not projected.')
		print(u'Doing the projection to UTM...')
		print(u'It could be long...')            
		#xstep, boxwidth, Coord, srs, dst_filename, a_utm, b_utm, test_N = \
		xsteps, boxwidths, Coord, srs, dst_filename, a_utm, b_utm, test_N = \
		       project_raster(rasterfnme, 
		                      rasterdata, 
		                      src_geotrans, 
		                      srs, 
		                      A, 
		                      B, 
		                      Coord, 
		                      xsteps, 
		                      boxwidths)
	elif synthetic:
		print(u'Working on a synthetic DEM...')
		text_N = 1
	else:
		if 'south' in srs.ExportToProj4():
			text_N = 0
		else:
			text_N = 1	
	
	return xsteps, boxwidths, Coord, srs, dst_filename, a, b, a_utm, b_utm, test_N, shpbox, ulx, lrx, lry, uly, projdone
	
	