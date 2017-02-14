######!/usr/bin/env python
# -*- coding: utf-8 -*-

# Do divisions with Reals, not with integers
# Must be at the beginning of the file
from __future__ import division

## Import Python modules
## I have problems to install rasterio : it does not find gdal libraries... from kingchaos
modulesNames = ['sys', 'math', 'os', 'utm', 'warnings', 'rasterstats', 'shapely']
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
###############################################################################

def calc_length(A,B):       
	"""
	function to calcule the distance
	between two points A and B given by their projected coordinates (in m),   
	The result is in meters !
	
	INPUTS:
	   A = Coordinates of the point A
	   B = Coordinates of the point B
	OUTPUTS:
	   calc_lengthm = length between A and B in m
	USAGE:
	   calc_lengthm = calc_length(A,B)
	"""
	
	calc_lengthm = np.sqrt((A[0] - B[0]) ** 2 + (A[1] - B[1]) ** 2)
	
	return calc_lengthm

############################################################################
       
def bearing(XA,YA,XB,YB):       
	"""
	Function to calcule the bearing on a sphere
	between two points A and B given by their long/lat coordinates
	
	INPUTS:
	   XA, YA = A coordinates in lat-long
	   XB, YB = B coordinates in lat-long
	OUTPUTS:
	   angle = bearing from A to B
	USAGE:
	   angle = bearing(XA,YA,XB,YB)
	   
	Written by X. Robert, june 2014
    """   
	# Not used anymore
	
	R = 6371  # Rayon de la terre en km
	
	# Convert in Radians
	lon1, lat1, lon2, lat2 = map(math.radians, [XA, YA, XB, YB])
	
	dLon = lon2 - lon1
	y = np.sin(dLon) * np.cos(lat2)
	x = np.cos(lat1) * np.sin(lat2) \
		- np.sin(lat1) * np.cos(lat2) * np.cos(dLon)
	
	return math.degrees(math.atan2(y, x))

############################################################################

def project_points(A, B, srs, A_utm, B_utm, test_N, iii, ggg):
	"""
	Function to project A and B from Lat-Long to UTM coordinates
	
	INPUTS:
	   A, B = the two points to project in utm
	   srs = Coordinate system of A & B
	   A_utm, B_utm = the points where the value of the projection will be recorded
	   test_N = Boolean that tells if the points are north or south
	   iii = iteration number (rank in the A list)
	   ggg = iteration number (rank in the B list)
	OUTPUTS:
	   A_utm, B_utm = Output projected coordinates
	   srs_new = new projection system
	USAGE:
	   A_utm, B_utm, srs_new = project_points(A, B, A_utm, B_utm, test_N, iii, ggg)
	
	Written by X. Robert, june 2015
	"""
	
	# Check if the points of the profile are in the same UTM zone
	if utm.from_latlon(A[iii ,1], A[iii ,0])[2] != utm.from_latlon(B[ggg, 1], B[ggg, 0])[2] \
	    or utm.from_latlon(A[0, 1], A[0, 0])[2] != utm.from_latlon(A[ggg,1], A[ggg, 0])[2]:
		# If not, Raise Warning, and choose only one UTM zone
		if utm.from_latlon(A[iii ,1], A[iii ,0])[0] < utm.from_latlon(B[iii, 1], B[iii, 0])[0]:
			utm_zone = utm.from_latlon(A[0,1], A[0,0])[2]
		else:
			utm_zone = utm.from_latlon(B[0,1], B[0,0])[2]
		if iii == 0:
			print(u'WARNING : The profile is over two UTM zones ! \n'
			      u'          I am choosing one of them : UTM%i' % utm_zone)	
			print(u'          If you get this warning, you may check if on a GIS software if '
			      u'          the computed profile is at the right place on the DEM. \n'
			      u'          Do not be surprised if you get an error in the calcul of statistics, '
			      u'          that means that your projected profile is out of the bounds of the projected DEM. \n'
			      u'          To solve the problem, use as entry :\n'
			      u'               - The projected DEM, \n'
			      u'               - Coordinates of points that you extract on the projected DEM!')
		# Define the new spatial reference
		srs_new = osr.SpatialReference()
		# Set the projection to UTM with the zone that correspond to the datapoints
		srs_new.SetUTM(utm_zone, test_N)
		transformxx = osr.CoordinateTransformation(srs, srs_new)
		# Find coordinates of the new raster
		(A_utm [iii, 0], A_utm [iii, 1], junk) = transformxx.TransformPoint(A[iii,0], A[iii,1])
		(B_utm [ggg, 0], B_utm [ggg, 1], junk) = transformxx.TransformPoint(B[ggg,0], B[ggg,1])
	else:
		# If in the same UTM zone, calculate the new coordinates in UTM
		A_utm [iii, 0] = utm.from_latlon(A[iii,1], A[iii, 0])[0]
		A_utm [iii, 1] = utm.from_latlon(A[iii,1], A[iii, 0])[1]
		B_utm [ggg, 0] = utm.from_latlon(B[ggg,1], B[ggg, 0])[0]
		B_utm [ggg, 1] = utm.from_latlon(B[ggg,1], B[ggg, 0])[1] 		
		# Get the UTM zone
		utm_zone = utm.from_latlon(A[0,1], A[0,0])[2]
		# Define the new spatial reference
		srs_new = osr.SpatialReference()
		# Set the projection to UTM with the zone that correspond to the datapoints
		srs_new.SetUTM(utm_zone, test_N)
		
	#return A_utm[iii, 0], A_utm[iii, 1], B_utm[iii, 0], B_utm[iii, 1], srs_new
	return A_utm, B_utm, srs_new

############################################################################

def raster_projection(rasterfnme, rasterdata, src_geotrans,  srs, dst_filename, srs_new, kmdeg, mem = 'NULL'):
	"""
	Project raster. Return a memory raster if mem != 'NULL' or a tiff file if mem == 'NULL'
	
	INPUTS:
	   rasterfnme = input raster name
	   rasterdata = input raster data
	   src_geotrans = Coordinate transform vector
	   srs = input coordinate system
	   dst_filename = output raster name
	   srs_new = output coordinate system into to project the raster
	   kmdeg = factor to convert from degrees to km
	   mem = tells if the output raster should be stored in a file ( = 'NULL')
	          or in memory (= 'MEMORY')
	OUTPUTS:
	   if mem = 'NULL', no variable outpût, only a projected raster is recorded in a new file
	   if mem = 'MEMORY',
	          dst = Raster reprojected in Lat Long, stored in memory
	          out_srs  = coordinate system
	USAGE:
	   to project a raster in a new tiff file:
	       raster_projection(rasterfnme, rasterdata, srs, dst_filename, srs_new, mem = 'NULL')
	   to project a raster and to store it in Memory:
	        dst_tmp = raster_projection(rasterfnme, 
		                                rasterdata, 
		                                srs, 
		                                dst_filename,  
		                                latlong_srs, 
		                                mem = 'MEMORY')
	   
	Written by X. Robert, june 2015
	"""
	
	if mem == 'NULL':
		# Test if the input raster is in lat-long
		if int(srs.GetAttrValue(u"AUTHORITY", 1)) != 4326:
			raise NameError(u'ERROR : the input raster for the raster_projection function '
			         u'is not in geographic coordinates (i.e. EPSG != 4326)')
		# Test if raster already exists
		if path.isfile(dst_filename):
			# If yes, do not do the projection, skip it
			print(u'The raster %s has already be prejected, I will not do it another time, '
			      u'I am thus using %s' % (rasterfnme, dst_filename))
		else:
			#print rasterdata.GetMetadata()
			#print rasterdata.GetRasterBand(1).GetStatistics(True, True)
			# Define the new raster
			#band = rasterdata.GetRasterBand(1)
			dst = gdal.GetDriverByName('GTiff').Create(dst_filename,
													   rasterdata.RasterXSize,
													   rasterdata.RasterYSize,
													   1,
													   GDT_Float32)
			#    	                                   band.DataType)
			# Set up the transform (from http://pangea.stanford.edu/~samuelj/musings/dems-in-python-pt-2-projecting-dems.html)
			# Create the coordinate transform object
			transformxx = osr.CoordinateTransformation(srs, srs_new)
			# Find coordinates of the new raster
			(ulx, uly, ulz) = transformxx.TransformPoint(src_geotrans[0], src_geotrans[3])
			(lrx, lry, lrz) = transformxx.TransformPoint(src_geotrans[0] + src_geotrans[1] * rasterdata.RasterXSize,
														 src_geotrans[3] + src_geotrans[5] * rasterdata.RasterYSize)
			# Set the new geotransform
			geotransform_proj = (ulx, src_geotrans[1] * 1000 / kmdeg, src_geotrans[2], uly, src_geotrans[4], -src_geotrans[1] * 1000 / kmdeg)
			# Set the new metadata
			dst.SetMetadata(rasterdata.GetMetadata())
			# Give the same extension to the new raster (REM : is it good ????)
			dst.SetGeoTransform(geotransform_proj)
			# Apply the projection		
			dst.SetProjection(srs_new.ExportToWkt())		
			# Do the reprojection
			gdal.ReprojectImage(rasterdata, # InRaster
								dst,        # OutRaster
								srs.ExportToWkt(), # src_wkt
								srs_new.ExportToWkt(), # dst_wkt
								#gdal.GRA_CubicSpline)
								#GRA_CubicSpline)
								#GRA_NearestNeighbour)
								GRA_Bilinear)    # eResampleAlg		
			# clean the memory
			dst = None
		return
	else:
		# Define the new raster in Memory
		dst = gdal.GetDriverByName('MEM').Create('',
												rasterdata.RasterXSize,
												rasterdata.RasterYSize,
												1,
												GDT_Float32)
		# Set up the transform (from http://pangea.stanford.edu/~samuelj/musings/dems-in-python-pt-2-projecting-dems.html)
		# Create the coordinate transform object
		transformxx = osr.CoordinateTransformation(srs, srs_new)
		# Find coordinates of the new raster
		(ulx, uly, ulz) = transformxx.TransformPoint(src_geotrans[0], src_geotrans[3])
		(lrx, lry, lrz) = transformxx.TransformPoint(src_geotrans[0] + src_geotrans[1] * rasterdata.RasterXSize,
													 src_geotrans[3] + src_geotrans[5] * rasterdata.RasterYSize)
		# Set the new geotransform
		geotransform_proj = (ulx, src_geotrans[1] * 1000 / kmdeg, src_geotrans[2], uly, src_geotrans[4], -src_geotrans[1] * 1000 / kmdeg)
		# Set the new metadata
		dst.SetMetadata(rasterdata.GetMetadata())
		# Give the same extension to the new raster (REM : is it good ????)
		dst.SetGeoTransform(geotransform_proj)
		# Apply the projection		
		dst.SetProjection(srs_new.ExportToWkt())		
		# Do the reprojection
		gdal.ReprojectImage(rasterdata, # InRaster
							dst,        # OutRaster
							src.ExportToWkt(), # src_wkt
							srs_new.ExportToWkt(), # dst_wkt
							GRA_Bilinear)    # eResampleAlg		
		return dst
			
############################################################################

def project_raster(rasterfnme, rasterdata, src_geotrans, srs, A, B, Coord, xstep, boxwidth):
	"""
	Function to project the raster and the lines
	
 	INPUT:
 		- rasterfnme : name of the raster to project (with the extension)
 		- rasterdata : 
 		- src_geotrans : data Geotransform of the raster to project
 		- srs : coordinate system of the original raster
 		- A and B : Coordinates that define the different lines
 		- Coord : Set the coordinate system in which are given A and B
 		- xstep = size of the x stepping in the original raster metric
 		- boxwidth : y-width in the original raster metric of the box
 	OUTPUT:
 		- xstep = size of the x stepping in the projected raster metric
 		- boxwidth = y-width in the projected raster metric of the box
 		- Coord = new coordinate system
 		- srs = new projected system
 		- dst_filename
 		- a_utm, b_utm = Coordinates of the points A and B transformed in UTM
 	USAGE:
 	    xstep, boxwidth, Coord, srs, dst_filename, a_utm, b_utm =
 	       project_raster(rasterfnme, rasterdata, src_geotrans, srs, A, B, Coord, xstep, boxwidth)
 	
 	Written by X. Robert, June 2015
 	"""
 	
	A = np.array(A)
	B = np.array(B)
	A_utm = np.zeros(A.shape)
	B_utm = np.zeros(B.shape)
	
	# find new projection
	# If data in Lat-long:
	# Calcul the factor to convert from degrees to km
	R = 6378
	kmdeg = 180 /(math.pi * R * np.cos(np.radians(A[0,1])))
	if srs.GetAttrValue('projcs') == None:
	#if Coord == 'latlong':# or srs.GetAttrValue('projcs') == None:
		# Check if north or South
		if A[0,1] < 0:
			test_N = 0
		else:
			test_N = 1
		# Do the lat-long to UTM conversion
		for iii in range (0, A.shape[0]):
			A_utm, B_utm, srs_new = project_points(A, B, srs, A_utm, B_utm, test_N, iii, iii)
			
		# Set the output
		# New raster file name
		dst_filename = rasterfnme[0:-4] + '_projUTM.tif'
		
		raster_projection(rasterfnme, rasterdata, src_geotrans, srs, dst_filename, srs_new, kmdeg, mem = 'NULL')
		
	else:
		# set the projection reference for lat long
		latlong_srs = osr.SpatialReference()
		latlong_srs.ImportFromEPSG(4326)
		dst_filename = ''
		# convert to lat-long (in memory)
		dst_tmp = raster_projection(rasterfnme, 
		                            rasterdata,
		                            src_geotrans, 
		                            srs, 
		                            dst_filename,  
		                            latlong_srs,
		                            kmdeg, 
		                            mem = 'MEMORY')
		# Project A and B
		A_lat = np.zeros(A.shape)
		B_lat = np.zeros(B.shape)
		transformxx = osr.CoordinateTransformation(srs, latlong_srs)
		for iii in range (0, A.shape[0]):
			# Find coordinates of the new raster
			(A_lat [iii, 0], A_lat [iii, 1], junk) = transformxx.TransformPoint(A[iii,0], A[iii,1])
			(B_lat [iii, 0], B_lat [iii, 1], junk) = transformxx.TransformPoint(B[iii,0], B[iii,1])
	    # find if in northern or southern hemisphere
		if A_lat[0,1] < 0:
			test_N = 0
		else:
			test_N = 1
		for iii in range (0, A.shape[0]):
			A_utm, B_utm, srs_new = project_points(A_lat, B_lat, srs, A_utm, B_utm, test_N, iii, iii)
	    # Set the output
		# New raster file name
		dst_filename = rasterfnme[0:-4] + '_projUTM.tif'
		# Project the lat-long memory raster to UTM
		raster_projection(rasterfnme, dst_tmp, src_geotrans, srs, dst_filename, srs_new, kmdeg, mem = 'NULL')
		dst_tmp = None
	
	# change to numpy array
	a_utm = np.array(A_utm)
	b_utm = np.array(B_utm)	
	# Convert the parameters to meters
	xstep[:] = [zz * 1000 / kmdeg for zz in xstep]
	boxwidth[:] = [zz * 1000 / kmdeg for zz in boxwidth]
	# Change the coordinates system variable
	srs = srs_new
	# Change the coordinates description
	Coord = 'utm'
	
	return xstep, boxwidth, Coord, srs, dst_filename, a_utm, b_utm, test_N

############################################################################
def CopyBand(srcband, dstband):
	"""
	from http://www.mit.edu/afs.new/athena/software/fwtools/current/FWTools-linux-i686-3.1.0/usr/bin/gdal_fillnodata.py
	"""
	for line in range(srcband.YSize):
		line_data = srcband.ReadRaster( 0, line, srcband.XSize, 1 )
		dstband.WriteRaster( 0, line, srcband.XSize, 1, line_data,
		                     buf_type = srcband.DataType )

############################################################################

def removeNoData(inraster, dst_filename = 'None', mask = 'default', NoData = 'None'):
	"""
	function to remove NoData values from a raster
	
	Modification of http://www.mit.edu/afs.new/athena/software/fwtools/current/FWTools-linux-i686-3.1.0/usr/bin/gdal_fillnodata.py
	
	INPUTS:
	   inraster = raster data to clean
	   dst_filename = name of the raster to save
	   mask = see ogr.FillNodata.py
	   NoData = Value of the NoData. 
	            If 'None', I consider that the NoData value is under or equal to 0.0
	
	OUTPUTS:
		newraster = cleaned raster data
	"""
	
	print(u'Cleanning the DEM, removing holes...')
	
	# Options for gdal
	max_distance = 100
	smoothing_iterations = 0
	options = []
	quiet_flag = 0
	#src_filename = None
	src_band = 1

	# 	Verify we have next gen bindings with the sievefilter method.
	try:
		gdal.FillNodata
	except:
		print(u'')
		print(u'gdal.FillNodata() not available.  You are likely using "old gen"')
		print(u'bindings or an older version of the next gen bindings.')
		print(u'')
		sys.exit(1)
		
	srcband = inraster.GetRasterBand(src_band)
	
	if mask is 'default':
		maskband = srcband.GetMaskBand()
	elif mask is 'none':
		maskband = None
	else:
		mask_ds = gdal.Open( mask )
		maskband = mask_ds.GetRasterBand(1)
	
	#   Create output file if one is specified.
	if dst_filename is not None:
		
		drv = gdal.GetDriverByName('GTiff')
		dst_ds = drv.Create( dst_filename,inraster.RasterXSize, inraster.RasterYSize,1,
		                 srcband.DataType )
		wkt = inraster.GetProjection()
		if wkt != '':
			inraster.SetProjection( wkt )
		inraster.SetGeoTransform( inraster.GetGeoTransform() )
		
		dstband = inraster.GetRasterBand(1)
		CopyBand( srcband, dstband )	
	else:
		dstband = srcband
	
	#	Invoke algorithm.
	if quiet_flag:
		prog_func = None
	else:
		prog_func = gdal.TermProgress
	
	result = gdal.FillNodata(dstband, 
	                         maskband,
                             max_distance, 
                             smoothing_iterations, 
                             options,
                             callback = prog_func )
	# close rasters
	src_ds = None
	dst_ds = None
	mask_ds = None
	
	newraster = gdal.Open(dst_filename, gdal.GA_ReadOnly)
	
	return newraster






