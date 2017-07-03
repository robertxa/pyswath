pyswath
========

This code is to built topographic swath profiles through a region, with a plot of the density of elevation in the box defined around the profile.

It uses a DEM that is projected in UTM (best input), but will also work with geographic coordinates (lat-long).

This module may contain numerous bugs, and could probably be ameliorate and/or optimized. If you have any comments, do not hesitate to add a new branch or to contact the author.
To know the history, please check the file History.txt

Install
-------

To install it :

.. code-block:: bash

	pip install pyswath

To update it :

.. code-block:: bash

	pip install -U pyswath

If during the update you get a problem with the update of a dependency, you may try :

.. code-block:: bash

	pip install -U --no-deps pyswath


The module has been written and tested with Python 2.7, but not with Python 3.

Dependencies
------------
This code needs the following python modules and their dependencies, you may install them before the installation of the **pyswath** module:
	- math
	- copy
	- utm
	- rasterstats>=0.11
	- rasterio
	- shapely
	- numpy
	- matplotlib
	- osgeo/gdal
	- progress
	- time

Usage
-----

Inside a (i)python environnement:

To import the module:

.. code-block:: python

>>> from pyswath import swathp
	
To plot a swath profile [A,B] through the raster 'DEM/dem.tif':

.. code-block:: python

    >>> swathp(rasterfnme = 'DEM/dem.tif',A = [(-78.4,-9.3)], B = [(-77.5,-8.5)],Coord = 'latlong', xsteps = [0.02], boxwidths = [0.2], binsize = 20,title = 'CB')

.. danger:: 
	
	Be careful, the units used for the points definition should be the same than for the DEM, and also for xsteps and boxwidth !! 
	Everything in degrees if Coord = LatLong, and everything in meters if Coord = utm !!

Options/inputs
--------------

To use options or inputs, you need to set them as	

.. code-block:: python

    >>> swathp(option_name = option_value, [...])
	
Options/inputs are (option_names):
	1. rasterfnme: name of the raster to work with
					Add the full path to the raster. Personally, I like to store my rasters in a DEM/folder
					Be aware that if the input raster is not projected in UTM, the code will create a projected raster in the same folder.
						
						ex: ``rasterfnme = 'Dem/Dem_Fusion-Peru_projUTM.tif'``
					
					Default = ``None``
	2. A, B: Coordinates of the 2 points A and B defining the whole profile in the same projection system than the DEM
				If multipoints profile (see further), these two points define the two extreme points.
				Be careful, do no take oceans, otherwise, statistics will be biased
				If different profiles:
				
					ex: ``A = [(-78.255,-9.713),(,),...]``
					
					``B = [(-77.255,-9.218),(,),...]``
				
				each column corresponding to one different profile
				
				Default: ``A = None, B = None``
	3. Coord: Units of the points coordinates. 
			- If Lat/Long : latlong
			- If UTM : utmZONE
			- If the dem is a synthetic dem (This is to avoid the problem of projections
			  If it is set to True, change the value of Factor if needed!) : Coord = 'synthetic'
			- If else, give a projection name that is NOT 'latlong' or 'utmZONE' or 'synthetic'
			
				ex: ``Coord = 'utm'``
				
				``Coord = 'latlong'``
				
				``Coord = 'synthetic'``
			
			Default ``Coord = 'utm'``
	4. factor: Factor is to convert to km (generaly 1000)
				if the unit of the DEM is 'meters', factor = 1000
				
				if the unit of the DEM is 'kilometers', factor = 1
				
				if the DEM is synthetic, units are arbitrary, so it could be factor = 0.001 (Test it !)
				
				Default: ``factor = 1000``

	5. xstep: Stepping along the profile in the same projection/coordinates system than the DEM
				If more than one profile with different profiles: xsteps = [5000, 2000,...], each column corresponding to one different profile
				
				If all the profiles have the same xsteps, just use one column 
				
					ex: ``xsteps = [500]`` (here in m)
				
				Default ``xsteps = [5000]``
	6. boxwidth: with of the box around the profile from where are extracted the stats in the same units than the DEM (m if m; km if km; deg if deg)
				If several profiles with different profiles: boxwidths = [20000, 15000,...], each column corresponding to one different profile
				
				If all the profiles have the same boxwidth, just use one column 
				
					ex: ``boxwidths = [2000]`` (here in m)
				
				Default ``boxwidths = [20000]``
	
	7. binsize: altitude binsize (for the altitude frequency plot) in the same units than the DEM (m if m; km if km; deg if deg)
				
		ex: ``binsize = 100``(in m)
				
		Default ``binsize = 20``
	
	8. title: title of the graphic
			The name will also be used to define the name:
				- in which the shapefiles are stored
				- of the output files
			
					ex: ``title = 'Synth-Essai'``
			
			Default ``title = 'Swath_profile'``
	
	9. shpbox: Name of the shapefile in which we extract the profile
			Default: ``shpbox = 'shpbox.shp'``
	
	10. sizeplotx, sizeploty: size of the plot.
							Standard size is ``sizeplotx = 8`` and ``sizeploty = 6``
							If you want to use the default/automatic setting, just give the value None to the variables
								
							Default: ``sizeplotx = None``
								
							Default: ``sizeploty = None``
								
	11. densitymin, densitymax: set the density colorbar limits (between [0,1]).
								Set it to None, if you want to keep the automatic settings
								
								Default: ``densitymin = None``

								Default: ``densitymax = None``
									
	12. ylimmin, ylimmax: set the range of the Y axis in the graph.
						  Set it None, if you want to keep the automatic settings
						  
						  Default: ``ylimmin = None``
						  
						  Default: ``ylimmax = None``
	
	
	13. remNoData: Flag to remove (``True``) or not (``False``) the NoData values from a DEM
					2015/08 : does not work very well, avoid it for the moment.
					
					Default: ``remNoData = False``
					
	14. corrnan: Flag to correct (``True``) or no (``False``) the graph from the Nan values
				Be careful, it replace the NaN values with the min value of the frequency
					
				Default: ``corrnan = False``
					
	15. nodatav: value of the NoData
				Default: ``nodatav = 0.0``
	
	16. multipoints: Multipoints section, Flag to set a profil with multipoints
						- [False] = only two points (Default)
						- [True] = more than two points
						- [True, False,...] if several transects, i.e. one boolean for each transect
					
							ex: ``multipoints = [False]``
					
	17. nbpointsint: Multipoints section, number of intermediary points in the profile
					If different profiles : nbpointsint = [...,3,2,1]
					Choose the order of the profiles with a decreasing number of intermediary points to avoid error in the code
					
						ex: ``nbpointsint = [0]``
					
	18. pointint: C,D,...: Multipoints section, intermediary points in the profile, given from A to B
				Be aware of the order !
				
					Give the name C for the 1st intermediary point (C = [(-78.255,-9.713),(,),...])
					
					Give the name D for the 2nd intermediary point (D = [(-78.255,-9.713),(,),...])
					
					Give the name E for the 3rd intermediary point
					
					...
					
	19. pointsdic: Multipoints section :
					dictionnary to assign a number to the different points. It should contain the same number of lines than the number of points
					
						ex: ``pointsdic = {1 : C, 2 : D, 3 : E, 4 : F, ... : ...}``
					
	20. printpointsi: Multipoints section, Flag to print (True) or not (False) the position of the intermediary points on the profile
						
					ex: ``printpointsi = True``
						
				Default ``printpointsi = False``
					
	21. idensity: Flag to plot the density (``True``) or not (Default, ``False``)


Help files
----------

To get help in your (i)python environnement:

.. code-block:: python

	>>> help(swathp)

or simply:

.. code-block:: python

	>>> swathp()

Examples
--------

To plot a swath profile [A,B] through the raster 'DEM/dem.tif' that is in lat-long (not projected):

.. code-block:: python

    >>> swathp(rasterfnme = 'DEM/dem.tif',A = [(-78.4,-9.3)], B = [(-77.5,-8.5)], Coord = 'latlong', xsteps = [0.02], boxwidths = [0.2], binsize = 20, title = 'CB')

The previous line permits to build the graph :

.. image:: https://github.com/robertxa/pyswath/tree/master/pyswath/DEMswath2.png?raw=true
   :scale: 100 %
   :align: center

To plot a swath profile through the raster 'DEM/Nperu_proj.tif' that is projected to UTM zone 18S:

.. code-block:: python
    >>> swathp(rasterfnme = 'DEM/Nperu_proj.tif',A = [(162374,9299742)], B = [(321829,9399929)], Coord = 'utm',xsteps = [10000], boxwidths = [20000], binsize = 20, title = 'NPeru')
	
To plot 2 swath profiles though the raster 'DEM/dem.tif' that is in lat-long (not projected):

.. code-block:: python

    >>> swathp(rasterfnme = 'DEM/dem.tif',A = [(-78.4,-9.3),(-78.4,-8.0)], B = [(-77.5,-8.5),(-76.0,-9.2)], Coord = 'latlong', xsteps = [0.02], boxwidths = [0.2], binsize = 20, title = 'CB')

To plot 1 swath profile with an intermediary point (kink) through the raster 'DEM/NPeru_proj.tif' that is in Lat-Long:

.. code-block:: python

    >>> swathp(rasterfnme = 'DEM/Nperu_proj.tif',A = [(162374,9299742)], B = [(321829,9399929)], Coord = 'utm', xsteps = [10000], boxwidths = [20000], binsize = 20, title = 'NPeru', multipoints = [True], nbpointsint = [1], pointsdic = {1 : 'C'}, printpointsi = True, C = [(217433,9383481)])
			
Outputs
-------

Inside the working directory, the code build several folders :
	- Data/: For each profile, the code outputs XXXXX files in Data/:
		+ data_title_Nbprofile.txt: 
			* Column 1 = Distance along the profile
			* Column 2 = Altitude
			* Column 3 = Altitude frequency
		+ datamask_title_Nbprofile.txt
		+ falti_title_Nbprofile.txt: altitude frequency
		+ statslines_title_Nbprofile.txt: 
			* Column 1 = Distance along profile
			* Column 2 = Min altitude
			* Colunm 3 = Max altitude
			* Column 4 = Median altitude
			* Column 5 = Mean altitude
	- Graphs/: for each profile, the code outputs here the graphs in pdf
	- shpbox/ (defined in the Variable declaration): In this directory, for each profile (or sub-profile if there are intermediary points), the code outputs:
		+ a shapefile defining the line between the two points of the profile
		+ a shapefile that define the box in which the transect is extracted

Contact
-------

If needed, do not hesitate to add a new branch or to contact the author. 
Please, use `https://isterre.fr/spip.php?page=contact&id_auteur=303 <https://isterre.fr/spip.php?page=contact&id_auteur=303>`_

Licence
-------

This package is licenced with `CCby-nc <https://creativecommons.org/licenses/by-nc/2.0/>`_
