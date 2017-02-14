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
modulesNames = ['sys', 'math', 'os', 'utm', 'warnings', 'shapely']
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
from rasterstats import raster_stats, zonal_stats             # For stats on rasters
from shapely.geometry import Polygon, LineString, shape, box, MultiPolygon
from shapely import affinity
from shapely import speedups
if speedups.available:
	speedups.enable()
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
try:
	from progress.bar import Bar                     # For a progress bar in the terminal output
except ImportError:
	sys.exit(u"ERROR : Module progress not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")

from raster_tools import *
from plotgraph import *
############################################################################

def makeshape(A, B, boxwidth, shp, srs, xstep = 'NULL', nx = 'NULL' , ny = 'NULL', beta = 'NULL'):
	"""
	Create a shapefile around a profile defined the points A and B, 
	with the boxwidth given
	
	INPUTS:
	   A(Xa,Ya) and B(Xb,Yb) = Points defining the line around which we build the shapefile
	   boxwidth = width of the box around the profile A-B
	   shp = path+name of the output shapefile
	   srs = coordinate system of the shapefile to create
	   xstep = size of the x-stepping in the unit of the DEM/shapefiles
	   shp = name or the shapefile that will be produced
	   nx = number of points in the direction normal to the main profile
	   ny = number of points in the direction of the main profile
	   beta = angle between the parallel and the transect A-B  
	OUTPUTS:
	   nx = number of points in the direction normal to the main profile
	   ny = number of points in the direction of the main profile
	   beta = angle between the parallel and the transect A-B
	   shapefile shp
	USAGE:
	  nsteps, ny, beta = makeshape(A, B, boxwidth, shpbox2, srs, 
	                               xstep = xsteps, [nx], [ny], [beta])
	  or:
	  makeshape(A, B, boxwidth, shpbox2, srs, 
	            xstep = 'NULL', [nx], [ny], [beta])
	
	"""
	
	# Calcul the length of the box
	xdistm = calc_length(A,B)
	if xstep != 'NULL':
		xstep = xstep
		nx = np.int(xdistm / xstep)
		ny = np.int(boxwidth / xstep)
		# beta is the angle between the parallel and the transect
		if beta == 'NULL':
			beta = math.atan((B[1] - A[1]) / (B[0] - A[0]))
	
	# create the corner points of the shapefile
	aa = np.empty(2)
	bb = np.empty(2)
	cc = np.empty(2)
	dd = np.empty(2)
	
	# Recalcul the coordinates of the corner :
	xxx = boxwidth / 2 * np.sin(beta)
	yyy = boxwidth / 2 * np.cos(beta)
	
	aa[0] = A[0] - xxx
	aa[1] = A[1] + yyy
	dd[0] = A[0] + xxx
	dd[1] = A[1] - yyy
	cc[0] = B[0] + xxx
	cc[1] = B[1] - yyy
	bb[0] = B[0] - xxx
	bb[1] = B[1] + yyy
	
	# from http://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
	# Create the polygon to be used for the shapefile
	poly = Polygon([(aa[0],aa[1]), (bb[0],bb[1]), (cc[0],cc[1]), (dd[0],dd[1]), (aa[0],aa[1])])
	
	# Now convert it to a shapefile with OGR    
	driver = ogr.GetDriverByName('Esri Shapefile')
	ds = driver.CreateDataSource(shp)
	# Set the Reference
	layer = ds.CreateLayer('zoup', srs, ogr.wkbPolygon)
	# Add one attribute
	layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
	defn = layer.GetLayerDefn()
	
	## If there are multiple geometries, put the "for" loop here
	# Create a new feature (attribute and geometry)
	feat = ogr.Feature(defn)
	feat.SetField('id', 123)
	# Make a geometry, from Shapely object
	geom = ogr.CreateGeometryFromWkb(poly.wkb)
	feat.SetGeometry(geom)
	layer.CreateFeature(feat)
	feat = geom = None  # destroy these
	# Save and close everything
	ds = layer = feat = geom = None
	
	if xstep != 'NULL':
		#return nx, ny, x, y
		return nx, ny, beta
	else:
		return

############################################################################

def save_shape(A, B, shp, srs, iii, kkk = None):
	"""
	Create a shapefile around a profile defined the points A and B, 
	with the boxwidth given
	
	INPUTS:
	   A(Xa,Ya) and B(Xb,Yb) = Points defining the line around which we build the shapefile
	   shp = path+name of the output shapefile
	   srs = coordinate system of the shapefile to create
	   iii = profile number
	   kkk = sub-profile number
	OUTPUTS:
	   shapefile shp
	USAGE:
	   save_shape(A, B, shp, srs)
	
	"""

	## Create the shapefile points
	shp2 = shp[0:-4] + '_points.shp'
	print(u'   Saving points defining the profil in the shapefile %s' % shp2)
	# Now convert it to a shapefile with OGR    
	driver = ogr.GetDriverByName('Esri Shapefile')
	ds = driver.CreateDataSource(shp2)
	# Set the Reference
	layer = ds.CreateLayer('zoup', srs, ogr.wkbPoint)
	# Add one attribute
	layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
	layer.CreateField(ogr.FieldDefn('Name', ogr.OFTString))
	
	defn = layer.GetLayerDefn()
	
	# create point A
	point = ogr.Geometry(ogr.wkbPoint)
	point.SetPoint(0, A[0], A[1])
	# Create feature
	feature = ogr.Feature(defn)
	feature.SetGeometry(point)
	feature.SetFID(1)
	feature.SetField('id', 1)
	if kkk == None:
		feature.SetField('Name', 'A_' + str(iii))
	else:
		feature.SetField('Name', 'A_' + str(iii) + '_' + str(kkk))
	# Save feature
	layer.CreateFeature(feature)
	# Cleanup
	point.Destroy()
	feature.Destroy()

	# create point B
	point = ogr.Geometry(ogr.wkbPoint)
	point.SetPoint(0, B[0], B[1])
	# Create feature
	feature = ogr.Feature(defn)
	feature.SetGeometry(point)
	feature.SetFID(2)
	feature.SetField('id', 2)
	if kkk == None:
		feature.SetField('Name', 'B_' + str(iii))
	else:
		feature.SetField('Name', 'B_' + str(iii) + '_' + str(kkk))
	# Save feature
	layer.CreateFeature(feature)
	# Cleanup
	point.Destroy()
	feature.Destroy()

	# Save and close everything
	ds = layer = None
	
	return

############################################################################
		
def findline(A, B, Xi, Xii, inc, i, nsteps, dist):
	"""
	Function to find the sub-line to build the sub-shapefile
	
	INPUTS:
	   A = A coordinates (xai, yai)
	   B = B coordinates (xbi, ybi)
	   Xi = coordinates of the first point of the line
	   Xii = coordinates of the second point of the line
	   inc = distance between two points
	   i = step's number
	   nsteps = total number of steps
	   dist = nsteps-Array of distance along the profile.
	OUTPUTS:
	   Xi = coordinates of the first point of the line updated for the next step
	   Xii = coordinates of the second point of the line updated for the next step
	   dist = nsteps-Array of distance along the profile updated
	USAGE:
	   Xi, Xii, dist = findline(A, B, Xi, Xii, inc, i, nsteps, dist)
	
	Written by X. Robert, june 2014
	"""
	# beta is the angle between the parallel and the transect
	#beta = np.pi / 2 - math.radians(bearing(A[0], A[1], B[0], B[1]))
	# The previous line is true if coordinates in lat/long (deg)
	# For utm in meters, the true value is the next one
	beta = np.arctan((B[1]-A[1]) / (B[0]-A[0]))
	if B[0] > A[0]:
		beta = -beta
	if B[1] > A[1]:
		beta = -beta
	incx = inc
	incy = inc
	if A[0] > B[0]:
		incx = -incx
	if	A[1] > B[1]:
		incy = -incy
	if i == 0 :
		Xi = [A[0],A[1]]
		# 15/06/2015; XR; inc/2 replaced by inc	
		Xii = [Xi[0] + incx * np.cos(beta), Xi[1] + incy * np.sin(beta)]
	else:
		Xi = Xii
		Xii = [Xi[0] + incx * np.cos(beta), Xi[1] + incy * np.sin(beta)]
	if i > 0: dist[i] = dist[i-1] + calc_length((Xi[0], Xi[1]), (Xii[0], Xii[1]))
	
	return Xi, Xii, dist	

############################################################################

def main_xa(A, B, xsteps, boxwidth, binsize, title, shpbox, rasterfnme, srs, factor, 
            sizeplotx = None, sizeploty = None, densitymin = None, densitymax = None,
            corrnan = False, synthetic = False, iii = None, 
            multipoints = None, nbpointsint = None, kkk = None):
	"""
	Main function that extract the data and plot the graph for each sets of A and B
	
	INPUTS:
	   A = A coordinates (xai,yai)
	   B = B coordinates (xbi,ybi)
	   xsteps = steps for the distance stepping
	   boxwidth = width of the box inside which we extract the data
	   binsize = altitudes binsize to build the histogram/frequency of the altitudes
	   title = Title of the graphs/transects
	   shpbox = Name of the shapefile that define the box around the profile to extract
	   rasterfnme = Name of the raster from which we extract the profile
	   srs = Coordinate system of the input dem
	   factor = factor to apply to path from the projected dem measure system to km 
	            (generally = 1000)
	   corrnan = Flag to correct (True) or no (False) the graph from the Nan values
	   synthetic = True if synthetic dem from a SPEM
	   iii = transect's number, i.e. rank in the main list A or B
	   multipoints = flag to set a profile along a multiline with different points
	   nbpointsmulti = if multipoints, number of intermediary poitns
	   kkk = rank of the intermediary points
	OUTPUTS:
	   no variables
	USAGE:
	   main_xa(A, B, xsteps, boxwidth, binsize, title, shpbox, rasterfnme, srs, factor, [iii])
	
	Written by X. Robert - June 2014
	"""
	# build the line A-B
	xdist = calc_length(A,B) / factor
	
	print(u'   Profile between A(' + str(A[0]) + u', ' + str(A[1]) + u') and B(' + str(B[0]) + u', ' + str(B[1]) + u')')
	#print('   Length of the profile : ' + str(xdist) + ' km')
	print(u'   Build the box shapefile...  ')
		
	# Build the line shapefile
	lineshp = 'SHP_' + title + '/line_' + title + '_' + str(iii+1) + '.shp'
	if multipoints:
		lineshp = lineshp[0:-4] + '_sub_' + str(kkk+1) + '.shp'
	# from http://gis.stackexchange.com/questions/52705/how-to-write-shapely-geometries-to-shapefiles
	# Create the polygon to be used for the shapefile
	polyl = LineString([(A[0],A[1]), (B[0],B[1])])
	# Now convert it to a shapefile with OGR    
	driver = ogr.GetDriverByName('Esri Shapefile')
	ds = driver.CreateDataSource(lineshp)
	
	# Set the Reference
	# This is in srs_new
	layer = ds.CreateLayer('zoup', srs, ogr.wkbLineString)
	# Add one attribute
	layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
	defn = layer.GetLayerDefn()
	## If there are multiple geometries, put the "for" loop here
	# Create a new feature (attribute and geometry)
	feat = ogr.Feature(defn)
	feat.SetField('id', 123)
	# Make a geometry, from Shapely object
	geom = ogr.CreateGeometryFromWkb(polyl.wkb)
	feat.SetGeometry(geom)
	layer.CreateFeature(feat)
	feat = geom = None  # destroy these
	# Save and close everything
	ds = layer = feat = geom = None
	# Build the box shapefile around the main line
	shpbox2 = shpbox[0:-4] + str(iii + 1) + '.shp'	
	# build the general shapefile from the [AB] line
	# and calcul the number of steps between A and B :
	nsteps, ny, beta = makeshape(A, B, boxwidth, shpbox2, srs, xstep = xsteps)
	# Save the points in a shapefile
	save_shape(A, B, shpbox2, srs, iii, kkk)
	
	# extract the statistics data for the whole box:
	stats = zonal_stats(shpbox2,
    	            	rasterfnme,
        	            stats=['min', 'max'])                         

	# initiate the tables for the heat map
	dist = np.empty(nsteps)
	inc = xsteps
	dist = np.zeros(nsteps)
	dist[0] = 0
	altimin = [f['min'] for f in stats][0]
	altimax = [f['max'] for f in stats][0]
	Xi = [0,0]
	Xii = [0,0]
	# set minimum altitude to 0 if it is inferior
	if altimin < 0:
		altimin = 0
	if altimax == None or altimax < altimin:
		raise NameError(u"   ERROR : There is a problem with the values in the DEM. \n"
		         u"   Check it! ")
	altirg = np.arange(altimin, altimax + binsize, binsize)	
	falti = np.zeros((dist.shape[0], altirg.shape[0]-1))
	statslines = np.empty((dist.shape[0],5))
	
	print(u'   Data extraction...')
	if multipoints:
		sub = 'sub-'
	else:
		sub = ''
	# Define the progress-bar that is output in the terminal window
	#bar = Bar('   Processing %sprofile %i' % sub, (iii+1), max = nsteps, suffix = '%(index)d / %(max)d ')
	if multipoints:
		bar = Bar(u'   Processing profile %i.%i' % ((iii+1), (kkk+1)), max = nsteps+1, suffix = '%(index)d / %(max)d ')
	else:
		bar = Bar(u'   Processing profile %i' % (iii+1), max = nsteps+1, suffix = '%(index)d / %(max)d ')
	# Stepping along the line
	for i in range (0,nsteps):
		
		# first add a step in the progress-bar
		bar.next()

		# build the shapefile between A+i*istep and A+(i+1)*istep
		# find the line:
		Xi, Xii, dist = findline(A, B, Xi, Xii, inc, i, nsteps, dist)
		
		# Store the sub-shapefile in TMP/
		# If they are needed, comment the end of the code in the function clean() 
		#                     that erase them when the run is finished
		shp = 'TMP/shp' + str(i) + '.shp'
		# build the sub-shapefile between Xi and Xii
		junk = makeshape(Xi, Xii, boxwidth, shp, srs, beta = beta)
		
		# calcul the distance from A
		statslines[i,0] = calc_length(A, Xi)
		
		# calcul the stats
		# Zonal_stats replace raster_stats.
		stats = zonal_stats(vectors = shp,
			                raster = rasterfnme, 
			                nodata_value = -9999,
			                raster_out = True,
		                    stats = ['min', 'max', 'median', 'mean', 'count'])
		# report it in a distance table statslines as followed:
		#	- row 0 = distance
		#	- row 1 = min altitude value of the shapefile
		#	- row 2 = max altitude value of the shapefile
		#	- row 3 = median altitude value of the shapefile
		#	- row 4 = mean altitude value of the shapefile
		statslines[i,1] = [f['min'] for f in stats][0]
		statslines[i,2] = [f['max'] for f in stats][0]
		statslines[i,3] = [f['median'] for f in stats][0]
		statslines[i,4] = [f['mean'] for f in stats][0]
		# Compute the altitude histogram
		try:
			#falti[i,:], alt = np.histogram([f['mini_raster_array'] for f in stats][0], altirg, density = True)
			falti[i,:], alt = np.histogram([f['mini_raster_array'] for f in stats][0].compressed(), altirg, density = True)
		except KeyError, e:
			print(repr(e))
			raise NameError(u'Well... This error could come from a mismatch between one point of the profil and the region of the DEM \n'
			          u'Please, check if the ongoing profile is well defined (Plot few of the TMP/shapefile on the DEM to check that)')
			
		# norm the altitude density
		falti[i,:] = falti[i,:] * np.diff(alt)
		
	bar.next()
	# Finish the progress-bar
	bar.finish()
	
	# Change the altitude def for the plot : 
	# choose the middle of the interval between the step i and the step j+1
	alti = [(altirg[hh] + altirg[hh+1])/2 for hh in range (0, altirg.shape[0]-1)]
	# initiate data:
	data = np.empty([falti.size,3])
	# change the falti table to a table with 3 columns : dist, altitude, frequency
	# to be able to plot it as a heatmap (see the plt.hexbin plot).
	data[:,0] = np.repeat(dist, falti.shape[1])
	data[:,1] = np.tile(alti, falti.shape[0])
	data[:,2] = falti.flatten()
	# OLD:
	#for i in range (0,falti.shape[0]):		 
		#dist = np.vstack((dist,dist))
		#alti = np.vstack((alti,alti))
		#for j in range (0,falti.shape[1]):		
			#data[j+i*falti.shape[1],0] = dist[i]
			#data[j+i*falti.shape[1],1] = (altirg[j] + altirg[j+1])/2
			#data[j+i*falti.shape[1],2] = falti[i,j]
					
	if multipoints:
		return statslines, data, xdist
	else:
		# mask the data where it is out of the bounds statslines[:,1] and statslines[:,2]
		# mask the data where it is null (set null value to nan)
		datamask1 = np.ma.masked_where(data[:,2] == 0, data[:,2])
		# mask the data where it is closed to 1 to get a colorbar not overprinted by only 1 value
		datamask = np.ma.masked_where(datamask1 == 1, datamask1)
		# Apply the mask to the data
		for dou in range (0,2):
			data[:,dou] = np.ma.masked_where(np.ma.getmask(datamask), data[:,dou])
	
		# save data in the folder Data/
		# Be carefull, there is no header...	
		print(u'   Saving the data in Data/ ...')
		np.savetxt('Data/data_' + title + '_' + str(iii) + '.txt',data)
		np.savetxt('Data/datamask_' + title + '_' + str(iii) + '.txt',datamask)
		np.savetxt('Data/falti_' + title + '_' + str(iii) + '.txt',falti)
		np.savetxt('Data/statslines_' + title + '_' + str(iii) + '.txt',statslines)
		
		# Plot the graph
		plot_graph(data = data, 
	    	       datamask = datamask, 
		           statslines = statslines, 
		           title = title, 
	    	       xdist = xdist,
	    	       xstep = xsteps,
		           boxwidth = boxwidth, 
		           factor = factor,
		           iii = iii,
		           corrnan = corrnan,
		           sizeplotx = sizeplotx,
		           sizeploty = sizeploty,
		           densitymin = densitymin,
		           densitymax = densitymax,
	    	       synthetic = synthetic)
	# Clean the workspace
	clean()
	
	return

############################################################################

def clean():
	"""
	Function to clean the workspace after the plot of each profile
	It removes the temporary shapefiles from the TMP/ directory
	"""	
	# cleaning
	# erase the temporary shapefiles
	print(u'   Cleaning....')
	files = os.listdir("TMP")
	for f in files:
		full_path = os.path.join("TMP",f)
		if not os.path.isdir(full_path):
			os.remove(full_path) 
	os.rmdir("TMP")
	os.mkdir("TMP")
	print(' ')
	return

