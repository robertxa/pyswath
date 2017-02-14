######!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
##########################################################
#                                                        #  
#          Slope profile extraction script               #
#                                                        #  
#                 By Xavier Robert                       #
#           From Peter Van der Beek codes                #
#                Montreal, 09/10/18                      #
#                  Lima, june 2014                       #
#                  Lima, june 2015                       #
#                                                        #  
##########################################################

Written by Xavier Robert, june 2015

xavier.robert@univ-grenoble-alpes.fr

(c) licence CCby-nc : http://creativecommons.org/licenses/by-nc/3.0/ 2015

"""

# Do divisions with Reals, not with integers
# Must be at the beginning of the file
from __future__ import division

# Import Python modules
# I have problems to install rasterio : it does not find gdal libraries... from kingchaos
#modulesNames = ['sys', 'math', 'os', 'utm', 'warnings', 'rasterstats', 'shapely', 'copy', 'time', 'rasterio']
modulesNames = ['sys', 'warnings', 'time', 'copy']
for module in modulesNames:
	try:
		# because we want to import using a variable, do it this way
		module_obj = __import__(module)
		# create a global object containging our module
		globals()[module] = module_obj
	except ImportError:
		sys.exit("ERROR : Module " + module + " not present. \n\n Please, install it \
			      \n\n Edit the source code for more information")
try:
	import numpy as np                               # need version 1.7 or higher
except ImportError:
	sys.exit("ERROR : Module Numpy not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	from osgeo import gdal, gdalnumeric, ogr, osr    # For GIS operations
	from osgeo.gdalconst import *
except ImportError:
	sys.exit("ERROR : Module osgeo/gdal not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	from progress.bar import Bar                     # For a progress bar in the terminal output
except ImportError:
	sys.exit("ERROR : Module progress not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")

from raster_tools import *
from profiles import *
from plotgraph import *
from checks import *

###############################################################################

def swathp(rasterfnme = None, A = None, B = None, Coord = 'utm', factor = 1000,
			xsteps = [5000], boxwidths = [20000], binsize = 20, title = 'Swath_profile',
			shpbox = 'shpbox.shp', sizeplotx = None, sizeploty = None,
			densitymin = None, densitymax = None,
			remNoData = False, corrnan = False, nodatav = 0.0,
			multipoints = [False], nbpointsint = [0], pointsdic = None, printpointsi = False,
			C = None, D = None, E = None, F = None, G = None, H = None, I = None, J = None,
			idensity = False):
	"""
	
	This code is to built topographic swath profile through a region, 
	with a plot of the density of elevation in the box defined around the profile.
	
	It uses a DEM that is projected in UTM (best input), 
	but will also work with geographic coordinates (lat-long).
	
	USAGE
	Inside a (i)python environnement:
	To import the module:
		>>> from swath import swathp
	To plot a swath profile [A,B] through the raster 'dem.tif':
		>>> swathp()

	INPUTS
	To use options or inputs, you need to set them as
	swath(option_name = option_value, [...])
	
	1. rasterfnme: name of the raster to work with
					Add the full path to the raster. Personally, I like to store my rasters in a DEM/folder
					Be aware that if the input raster is not projected in UTM, the code will create a projected raster in the same folder.
					ex: rasterfnme = 'Dem/Dem_Fusion-Peru_projUTM.tif'
					Default = None
	2. A, B: Coordinates of the 2 points A and B defining the whole profile in the same projection system than the DEM
				If multipoints profile (see further), these two points define the two extreme points.
				Be careful, do no take oceans, otherwise, statistics will be biased
				If different profiles:
				ex: A = [(-78.255,-9.713),(,),...]
					B = [(-77.255,-9.218),(,),...]
				each column corresponding to one different profile
				Default: A = None, B = None
	3. Coord: Units of the points coordinates. 
			- If Lat/Long : Coord = 'latlong'
			- If UTM : Coord = 'utmZONE'
			- If the dem is a synthetic dem (This is to avoid the problem of projections
			  If it is set to True, change the value of Factor if needed!) : Coord = 'synthetic'
			- If else, give a projection name that is NOT 'latlong' or 'utmZONE' or 'synthetic'
			ex: Coord = 'utm'
				Coord = 'latlong'
				Coord = 'synthetic'
			Default = 'utm'
	4. factor: Factor is to convert to km (generaly 1000)
				if the unit of the DEM is 'meters', factor = 1000
				if the unit of the DEM is 'kilometers', factor = 1
				if the DEM is synthetic, units are arbitrary, so it could be factor = 0.001 (Test it !)
				Default: factor = 1000
				(ex: factor = 0.01)
	5. xsteps: Stepping along the profile in the same projection/coordinates system than the DEM
				If more than one profile with different profiles: xsteps = [5000, 2000,...], each column corresponding to one different profile
				If all the profiles have the same xsteps, just use one column 
				ex: xsteps = [5000]
				Default = [5000]
	6. boxwidths: with of the box around the profile from where are extracted the stats in the same units than the DEM (m if m; km if km; deg if deg)
				if several profiles with different profiles: boxwidths = [20000, 15000,...], each column corresponding to one different profile
				If all the profiles have the same boxwidth, just use one column 
				ex: boxwidths = [20000]
				Default = [20000]
	7. binsize: altitude binsize (for the altitude frequency plot) in the same units than the DEM (m if m; km if km; deg if deg)
				ex: binsize = 20
				Default = 20
	8. title: title of the graphic
			The name will also be used to define the name:
				- in which the shapefiles are stored
				- of the output files
			ex: title = 'Synth-Essai'
			Default = 'Swath_profile'
	9. shpbox: Name of the shapefile in which we extract the profile
			Default: shpbox = 'shpbox.shp'
	10. sizeplotx, sizeploty: size of the plot.
							Standard size is sizeplotx = 8 and sizeploty = 6
							If you want to use the default/automatic setting, just give the value None to the variables
							Default: sizeplotx = None
									 sizeploty = None
	11. densitymin, densitymax: set the density colorbar limits (between [0,1]).
								Set it to None, if you want to keep the automatic settings
								Default: densitymin = None
										 densitymax = None
	12. remNoData: Flag to remove (True) or not (False) the NoData values from a DEM
					2015/08 : does not work very well, avoid it for the moment.
					Default: remNoData = False
	13. corrnan: Flag to correct (True) or no (False) the graph from the Nan values
				Be careful, it replace the NaN values with the min value of the frequency
				Default: corrnan = False
	14. nodatav: value of the NoData
				Default: nodatav = 0.0
	15. multipoints: Multipoints section : 
					Flag to set a profil with multipoints
					[False] = only two points (Default)
					[True] = more than two points
					[True, False,...] if several  transects
					ex: multipoints = [False]
	16. nbpointsint: Multipoints section :
					number of intermediary points in the profile
					I choose to limit to 9 intermediary points.
					If different profiles : nbpointsint = [...,3,2,1]
					Choose the order of the profiles with a decreasing number of intermediary points to avoid error in the code
					ex: nbpointsint = [0]
	17. pointsdic: Multipoints section :
					dictionnary to assign a number to the different points. It should contain the same number of lines than the number of points
				ex: pointsdic = {1 : C,
								2 : D,
								3 : E,
								4 : F,
								... : ...
								}
	18. printpointsi: Multipoints section :
					Flag to print (True) or not (False) the position of the intermediary points on the profile
					ex: printpointsi = True
					Default = False
	19. C,D,...,J: Multipoints section :
				intermediary points in the profile, given from A to B
				Be aware of the order !
				I choose to limit to 9 intermediary points C,D,E,F,G,H,I,J
				Give the name C for the 1st intermediary point (C = [(-78.255,-9.713),(,),...])
							  D for the 2nd intermediary point (D = [(-78.255,-9.713),(,),...])
							  E for the 3rd intermediary point
							  ...
	20. idensity: Flag to plot the density (True) or not (Default, False)
		
	OUPUTS
	Inside the working directory, the code build several folders :
		- Data/: For each profile, the code outputs XXXXX files in Data/:
			+ data_title_Nbprofile.txt:
			+ datamask_title_Nbprofile.txt:
			+ falti_title_Nbprofile.txt:
			+ statslines_title_Nbprofile.txt:	
		- Graphs/: for each profile, the code outputs here the graphs in pdf
		- shpbox/ (defined in the Variable declaration): In this directory, for each profile, the code outputs:
			+ a shapefile defining the line between the two points of the profile
	    	+ a shapefile the define the box in which the transect is extracted

	CONTACT
		If needed, do not hesitate to contact the author. 
		Please, use https://isterre.fr/spip.php?page=contact&id_auteur=303

	LICENCE
	This package is licenced with `CCby-nc` (https://creativecommons.org/licenses/by-nc/2.0/)
	"""
	start_time = time.time()
	
	print(' ')
	print(' ')
	print('******************************************')
	print('Program to extract swath profiles from DEM')
	print('     Written by X. Robert, ISTerre')
	print('           June 2014-2015      ')
	print('******************************************')
	
	# If the minima input data are not given print the help file
	if rasterfnme == None or A == None or B == None: help(swathp)
	
	# Check the files... 
	if Coord == 'synthetic': 
		synthetic = True
	else:
		synthetic = False
	xsteps, boxwidths, Coord, srs, dst_filename, a, b, a_utm, b_utm, test_N, shpbox, ulx, lrx, lry, uly, projdone = checkfiles(rasterfnme, 
						A, B, xsteps, boxwidths, shpbox, title, 
						Coord, synthetic,
						multipoints, remNoData)
	
	print(u'  ')
	print(u'Number of profiles : ' + str(a.shape[0]))
	print(u'With xsteps = %s and boxwidths = %s' % (xsteps, boxwidths)) 
	
	# Begin the stepping on the number of profiles (a.shape[0])
	for iii in range (0, a.shape[0]):
		print(u' ')
		print (u'Doing swath profile ' + str(iii + 1) + u'/' + str(a.shape[0]))
		if len(xsteps) != 1:
			xstep = xsteps[iii]
		else:
			xstep = xsteps[0]
		if len(boxwidths) != 1:
			boxwidth = boxwidths[iii]
		else:
			boxwidth = boxwidths[0]
		aa = (a_utm[iii,0], a_utm[iii,1])
		bb = (b_utm[iii,0], b_utm[iii,1])
		
		multip = False
		if len(multipoints) != 1:
			if multipoints[iii]:
				ggg = -1
				for indexg in range (0,iii+1):
					if multipoints[indexg]:
						ggg += 1
			multip = multipoints[iii]
		else:
			ggg = 0
			multip = multipoints[0]

		if multip:
			print(u'   profile with %i intermediary points' % nbpointsint[ggg])
			# initiate the first point with A
			aa1 = aa
			A = np.array(A)
			xdisttot = 0
			if printpointsi:
				xpoint = []
			#cc = np.zeros((nbpointsint[iii], 2))
			for kkk in range (0, nbpointsint[ggg]+1):
				D_utm = np.zeros(A.shape)
				if kkk != nbpointsint[ggg]:
					# if this is not the end of the profile, set the second point to cc[kkk]											
					c = pointsdic[kkk+1]
					c = eval(c)
					#ggg = -1
					#for indexg in range (0,iii+1):
					#	if multipoints[indexg]:
					#		ggg += 1
					
					# check if the point is in the profile
					if  not (ulx < c[ggg][0]< lrx):
						raise NameError(u"ERROR : Intermediary point%s is outside the DEM border %s" % (c, (ulx, lrx)))
					if  not (lry < c[ggg][1]< uly):
						raise NameError(u"ERROR : Intermediary point%s is outside the DEM border %s" % (c, (lry, uly)))
					
					if projdone:
						C2_utm = np.zeros(c.shape)
						junk = D_utm
						# do the projection if needed
						junk, c, junk2 = project_points(A, c, D_utm, C2_utm, test_N, iii, ggg)
					bb1 = c[ggg]
				else:
					# if this is the end of the profile, set the second point to B
					bb1 = bb
					
				# call the main code for each sub-profile
				if kkk == 0:
					shpbox = shpbox[0:-4] + '_sub_' + str(kkk+1) + '.shp'
				else:
					shpbox = shpbox[0:-5] + str(kkk+1) + '.shp'
				time_st = time.time()
				statslines, data, xdist = main_xa(A = aa1, B = bb1,
			                                      xsteps = xstep, 
			                                      boxwidth = boxwidth, 
			                                      binsize = binsize, 
			                                      title = title,
			                                      shpbox = shpbox, 
			                                      rasterfnme = dst_filename, 
			                                      srs = srs, 
			                                      factor = factor,
			                                      corrnan = corrnan,
			                                      synthetic = synthetic, 
			                                      iii = iii,
			                                      multipoints = multip,
			                                      nbpointsint = nbpointsint[ggg],
			                                      kkk = kkk)
				# update the first point of the sub-profile
				aa1 = c[ggg]
				xdisttot = xdisttot + xdist
				# create backup files
				if kkk == 0:
					statslinestot = copy.copy(statslines)
					datatot = copy.copy(data)
				else:
					print(u'   Concatenating profiles... \n'
					      u'       It can be long... Be patient')
					try: 
						statslinestotbackup = copy.copy(statslinestot)
						del statslinestot
					except NameError:
						raise NameError(u'ERROR: Statslinetot does not exists...')
					try:
						datatotbackup = datatot
						del datatot
					except NameError:
						raise NameError(u'ERROR: datatot does not exists...')	
					statslinestot = np.empty((statslinestotbackup.shape[0] + statslines.shape[0],5))
					datatot = np.empty((datatotbackup.shape[0] + data.shape[0]-1,3))
					if printpointsi:
						xpoint.append(max(statslinestotbackup[:,0]))
					# Stacking the files
					statslinestot = np.concatenate((statslinestotbackup, statslines), axis = 0)
					for lll in range (statslinestotbackup.shape[0],statslinestot.shape[0]):
						statslinestot[lll,0] = max(statslinestotbackup[:,0]) \
							                     + statslines[lll - statslinestotbackup.shape[0],0]
					datatot = np.concatenate((datatotbackup, data), axis = 0)
					for lll in range (datatotbackup.shape[0], datatot.shape[0]):
						# This is this loop that is long...
						datatot[lll,0] = max(datatotbackup[:,0]) + data[lll - datatotbackup.shape[0],0]
					time_int = time.time() - time_st
					print(u'       This subprofil needed %i seconds to run' % time_int)
				# Clean the workspace
				clean()
			# mask the data where it is out of the bounds statslines[:,1] and statslines[:,2]
			# mask the data where it is null (set null value to nan)
			datamask1 = np.ma.masked_where(datatot[:,2] == 0, datatot[:,2])
			# mask the data where it is closed to 1 to get a colorbar not overprinted by only 1 value
			datamask = np.ma.masked_where(datamask1 == 1, datamask1)
	
			# save data in the folder Data/
			# Be carefull, there is no header...	
			print(u'   Saving the data in Data/ ...')
			np.savetxt('Data/data_' + title + '_' + str(iii) + '.txt',datatot)
			np.savetxt('Data/datamask_' + title + '_' + str(iii) + '.txt',datamask)
			#np.savetxt('Data/falti_' + title + '_' + str(iii) + '.txt',falti)
			np.savetxt('Data/statslines_' + title + '_' + str(iii) + '.txt',statslinestot)
			
			# Plot the graph
			plot_graph(data = datatot, 
				       datamask = datamask, 
		    	       statslines = statslinestot, 
			           title = title, 
		    	       xdist = xdisttot,
		    	       xstep = xstep,
			           boxwidth = boxwidth, 
				       factor = factor,
				       iii = iii,
				       corrnan = corrnan,
				       sizeplotx = sizeplotx,
				       sizeploty = sizeploty,
				       densitymin = densitymin,
				       densitymax = densitymax,
				       synthetic = synthetic,
				       xpoint = xpoint)	
		else:
			# call the main code
			main_xa(A = aa, B = bb,
			        xsteps = xstep, 
			        boxwidth = boxwidth, 
			        binsize = binsize, 
			        title = title, 
			        shpbox = shpbox, 
			        rasterfnme = dst_filename, 
			        srs = srs, 
			        factor = factor,
			        sizeplotx = sizeplotx,
			        sizeploty = sizeploty,
			        densitymin = densitymin,
			        densitymax = densitymax,
			        corrnan = corrnan,
			        synthetic = synthetic, 
			        iii = iii)

	print('  ')
	elapsed_time = time.time() - start_time
	m, s = divmod(elapsed_time, 60)
	h, m = divmod(m, 60)
	print(u'Finished in %d hours, %d minutes and %d seconds...' %  (h, m, s))
	print(u'_______________________________________________________________________________')
	print(u' ')
	print(u' ')


if __name__ == u'__main__':
	
	
	# initiate variables
	
	# run the transformation  
	swathp()

# END of program
# Pffffiouuuu