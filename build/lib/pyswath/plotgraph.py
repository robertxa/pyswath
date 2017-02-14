######!/usr/bin/env python
# -*- coding: utf-8 -*-

# Do divisions with Reals, not with integers
# Must be at the beginning of the file
from __future__ import division

## Import Python modules
## I have problems to install rasterio : it does not find gdal libraries... from kingchaos
modulesNames = ['sys']
for module in modulesNames:
	try:
		# because we want to import using a variable, do it this way
		module_obj = __import__(module)
		# create a global object containging our module
		globals()[module] = module_obj
	except ImportError:
		sys.exit(u"ERROR : Module " + module + " not present. \n\n Please, install it \
			      \n\n Edit the source code for more information")
try:
	import numpy as np                               # need version 1.7 or higher
except ImportError:
	sys.exit(u"ERROR : Module Numpy not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")
try:
	import matplotlib.pyplot as plt                  # module to plot figures
	from matplotlib import cm
except ImportError:
	sys.exit(u"ERROR : Module matplotlib not present. \n\n Please, install it \
		      \n\n Edit the source code for more information")

############################################################################
	
def plot_graph(data, datamask, statslines, title, xdist, xstep, boxwidth, factor, iii, corrnan, \
               sizeplotx = None, sizeploty = None, densitymin = None, densitymax = None,
               synthetic = False, xpoint = None):
	"""
	Function to plot swath profile
	
	INPUTS:
	   data = frequency of the altitude
	   datamask = mask the data where it is null (set null values to nan)
	              used to build the frequency range
	   statslines = min, mac, mean, median data to plot
	   title = title of the graph
	   xdist = length of the transect
	   factor = factor to apply to plot the x-axes in km
	   corrnan = Flag to correct (True) or no (False) the graph from the Nan values
	   synthetic = True if synthetic dem
	   xpoint = distance along profile of the intermediary points if it exists
	OUTPUTS:
	   Graph in pdf
	USAGE:
	   plot_graph(data, datamask, statslines, title, xdist, factor)
	
	"""	
	# Do the plot:
	print(u'   Plot the transect...')	
	
	linewidth = 1
	
	plt.clf()
	# Change the settings of the graph
	if sizeplotx != None and sizeploty != None:
		plt.figure(num = None, figsize=(sizeplotx,sizeploty))
	
	## plot the altitude density plot in red
	plt.hexbin(data[:,0] / factor, data[:,1], C = datamask, cmap = cm.Reds, bins = None)
	## make a color palette with the altitude frequency
	# to force the boundaries of the color palette :
	if densitymin != None and densitymax != None:
		plt.clim(densitymin, densitymax)
	cd = plt.colorbar()
	cd.set_label(u'Altitude frequency')
	
	#if synthetic:
	#	factor = 1 / factor
	
	# plot the max, the mean, and the min altitude along the profile
	#     You can change from mean to median altitude by commenting the plt.plot mean block
	#     and de-commenting the plt.plot Median plot
	plt.plot(statslines[:,0] / factor, statslines[:,2], 
	         color = 'blue', 
	         linestyle = '-',
	         linewidth = linewidth, 
	         label = u'Max alt')
	#plt.plot(statslines[:,0], statslines[:,3], 
	#         color = 'black', 
	#         linestyle = '-', 
	#         linewidth = 2
	#         label = 'Median alt')
	
	plt.plot(statslines[:,0] / factor, statslines[:,4], 
	         color = 'black', 
	         linewidth = 2, 
	         linestyle = '-', 
	         label = u'Mean alt')
	
	# remove NoData (valid only for the min ; I assume that NaN data are data inf. or equal to 0)
	# TO BE CHANGED	
	statslines[statslines[:,1] <= 0.0, 1] = np.nan
	
	#altifreq = np.ma.masked_where(data[:,2] == 0, data[:,1])
	altifreq1 = data
	altifreq1[altifreq1[:,2] == 0.0, 1] = np.nan
	altifreq1[altifreq1[:,1] == 0.0, 1] = np.nan
	altifreq = altifreq1[:,1]
	
	if corrnan:
		ploty = np.zeros(statslines.shape[0])
		for k in range (0, ploty.shape[0]):
			for jj in range (int((k-1) * altifreq.shape[0]/statslines.shape[0]), 
			                 int(k * altifreq.shape[0]/statslines.shape[0])):
				if np.isnan(altifreq[jj]):
					altifreq[jj] = altifreq[jj+1]
			
			altifreqmin = np.nanmin(altifreq[k * altifreq.shape[0]/statslines.shape[0] : 
			                           (k +1) * altifreq.shape[0]/statslines.shape[0]])
			if not np.isnan(statslines[k,1]):
				ploty[k] = min(statslines[k,1], altifreqmin)
			else:
				ploty[k] = altifreqmin
		
		plt.plot(statslines[:,0] / factor, ploty, 
		#plt.plot(statslines[:,0] / factor, statslines[:,1], 
	    	     color = 'green', 
	        	 linestyle = '-', 
	        	 linewidth = linewidth,
		         label = u'Min alt')
	else:
		 plt.plot(statslines[:,0] / factor, statslines[:,1], 
	    	     color = 'green', 
	        	 linestyle = '-', 
	        	 linewidth = linewidth,
		         label = u'Min alt')
	if xpoint != None:
		for jjj in range(0,len(xpoint)):
			plt.axvline(x = xpoint[jjj]/factor, ls = '--', linewidth = 1, color = '0.75')#, label = 'STDS position')
	
	#plt.ylim(0,statslines[:,2].max() + 500)
	plt.ylim(0,statslines[:,2].max())
	plt.xlim(0,xdist)
	
	plt.xlabel(u'Distance along profile (km)')
	plt.ylabel(u'Altitude (m)')
	plt.legend(loc='best', numpoints = 1)
	plt.title(title + ' ' + str(iii + 1) + u' (xstep = ' + str(round(xstep,1)) + u' m; boxwidth = ' + str(round(boxwidth,1)) + u' m)')
	plt.savefig("Graphs/" + title + '_transect_' + str(iii + 1) + '.pdf')

	plt.close()
	
	print(u'   Saved in Graphs/' + title + u'_transect_' + str(iii + 1) + u'.pdf')

	return
	