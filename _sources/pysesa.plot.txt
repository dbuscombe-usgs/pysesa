.. _pysesa.plot:

pysesa.plot module
======================

Initialise the pysesa plot class
Takes a file output from pysesa.process and allows a number
of different 2d or 3d plots of the outputs

Syntax
----------
p = pysesa.plot()
p = pysesa.plot('/home/my_pysesa_output_file.xyz')

Parameters
-----------
pysesa_file : str
	pysesa::process output file

If no arguments given, it prompts you to choose a pysesa::process output file

Returns
----------
self : instance
      pysesa.plot instance

DATA functions
---------------
pc = p.get_pc() : ndarray
      NxM contents of pysesa_file
  
xyz = p.get_xyz() : ndarray        
      Nx3 contents of raw point cloud (the file processed by pysesa::process)

vars = p.parse_pc_vars() : dict
      NxM contents of pysesa_file parsed into dict object
      1 key per variable in p.get_pc()
      vars.keys() returns list of variables in dict

2D plotting functions (mayavi not required)
--------------------------------------------

p.grd_xyz()
------------- 
   produces 2d plot of the gridded [x,y,z] surface made from decimated point cloud,
   as returned by parse_pc_vars()

   Syntax
   ----------
   [] = p.grd_xyz(res, azimuth, altitude, zf, cmap, dpi, alpha, ticksize, labelsize)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

   Optional Parameters
   --------------------
   res : float, *optional* [default = 0.1]
   	   grid resolution

   azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see `here <http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm>`_

   altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see `here <http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm>`_

   zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

   cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented `here <http://matplotlib.org/examples/color/colormaps_reference.html>`_
   	   
   dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

   alpha : float, *optional* [default = 0.5]
           transparency, between 0.0 and 1.0

   ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

   labelsize : int, *optional* [default = 6]
           size of x and y axes labels

p.grd_var()
------------      
   produces 2d plot of the gridded surface made from 1 output variable in p.parse_pc_vars()
   e.g. p.grd_var('sigma')

   Syntax
   ----------
      [] = p.grd_var(var, res, azimuth, altitude, zf, cmap, dpi, log_scale, smooth, filtsz, alpha, ticksize, labelsize)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

   var : str
   	   name of variable in p.parse_pc_vars() that will be plotted
           e.g. p.grd_var('sigma')


   Optional Parameters
   --------------------
   res : float, *optional* [default = 0.1]
   	   grid resolution

   azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see `here <http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm>`_

   altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see `here <http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm>`_

   zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

   cmap : str, *optional* [default = 'hot']
 	   colormap
           possible colormaps are documented `here <http://matplotlib.org/examples/color/colormaps_reference.html>`_
   	   
   dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

   log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

   smooth : bool, *optional* [default = True]
   	   if True, will smooth plotted dependent variable
           with a median filter and window size specified by filtsz (below)

   filtsz : int, *optional* [default = 3]
   	   size of filter (pixels) if smooth==1

   alpha : float, *optional* [default = 0.5]
           transparency, between 0.0 and 1.0

   ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

   labelsize : int, *optional* [default = 6]
           size of x and y axes labels


p.grd_vars()
-------------
   produces a 2d plot of the gridded surface made from each output variable in p.parse_pc_vars()

   Syntax
   ----------
   [] = p.grd_vars(res, azimuth, altitude, zf, cmap, dpi, log_scale, smooth, filtsz, alpha, ticksize, labelsize)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

   Optional Parameters
   --------------------
   res : float, *optional* [default = 0.1]
   	   grid resolution

   azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see `here <http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm>`_

   altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see `here <http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm>`_

   zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

   cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented `here <http://matplotlib.org/examples/color/colormaps_reference.html>`_
   	   
   dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

   log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

   smooth : bool, *optional* [default = True]
   	   if True, will smooth plotted dependent variable
           with a median filter and window size specified by filtsz (below)

   filtsz : int, *optional* [default = 3]
   	   size of filter (pixels) if smooth==1

   alpha : float, *optional* [default = 0.5]
           transparency, between 0.0 and 1.0

   ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

   labelsize : int, *optional* [default = 6]
           size of x and y axes labels


3D plotting functions (mayavi not required)
--------------------------------------------

p.plt_xyz()
------------ 
   produces 3d plot of Nx3 contents of raw point cloud, as returned by p.get_xyz()

   Syntax
   ----------
   [] = p.plt_xyz(elev, azim, markersize, dpi, ticksize, labelsize)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

   Optional Parameters
   --------------------
   elev : float, *optional* [default = 65]
           the elevation angle in the z plane

   azimuth : float, *optional* [default = -115]
           azimuth angle in the x,y plane 

   markersize : float, *optional* [default = 0.01]
           marker size in x and y axes units

   dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

   ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

   labelsize : int, *optional* [default = 6]
           size of x and y axes labels       


p.plt_xy_var()  
--------------- 
   produces 3d plot of 1 output variable in p.parse_pc_vars(), e.g. p.grd_var('sigma')

   Syntax
   ----------
   [] = p.plt_xy_var(var, log_scale, dpi, markersize, ticksize, labelsize, elev, azim)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

   var : str
   	   name of variable in p.parse_pc_vars() that will be plotted
           e.g. p.grd_var('sigma')

   Optional Parameters
   --------------------
   log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

   dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

   markersize : float, *optional* [default = 5]
           marker size in points^2. See `here <http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter>`_

   ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

   labelsize : int, *optional* [default = 6]
           size of x and y axes labels

   elev : float, *optional* [default = 65]
           the elevation angle in the z plane

   azimuth : float, *optional* [default = -115]
           azimuth angle in the x,y plane          


p.plt_xy_vars()
----------------
   produces a 3d plot of each output variable in p.parse_pc_vars() 

   Syntax
   ----------
   [] = p.plt_xy_vars(log_scale, dpi, markersize, ticksize, labelsize, elev, azim)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

   Optional Parameters
   --------------------
   log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

   dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

   markersize : float, *optional* [default = 5]
           marker size in points^2. See `here <http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter>`_

   ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

   labelsize : int, *optional* [default = 6]
           size of x and y axes labels

   elev : float, *optional* [default = 65]
           the elevation angle in the z plane

   azimuth : float, *optional* [default = -115]
           azimuth angle in the x,y plane        
   

3D plotting functions (requires mayavi) 
----------------------------------------

p.grd_xyz3d() 
-------------  
   produces 3d plot of the gridded surface made from the Nx3 contents of raw point cloud,
   as returned by p.get_xyz()

   Syntax
   ----------
   [] = p.grd_xyz3d(res, cmap, pitch, azimuth, distance, xsize, ysize)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()


   Optional Parameters
   --------------------
   res : float, *optional* [default = 0.1]
   	   grid resolution

   cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here `here <http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html>`_

   pitch : float, *optional* [default = 10]
   	   rotates the camera. see `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch>`_

   azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360). See `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view>`_

   distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame. See `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view>`_

   xsize : int, *optional* [default = 2000]
   	   size (number of pixels) of output image in x dimension

   ysize : int, *optional* [default = 1000]
   	   size (number of pixels) of output image in y dimension
   
p.grd_var_3d()  
---------------
   produces 3d plot of the gridded surface made from 1 output variable in p.parse_pc_vars()

   Syntax
   ----------
   [] = p.grd_var_3d(var, res, cmap, pitch, azimuth, distance, log_scale, smooth, filtsz, xsize, ysize)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

   var : str
   	   name of variable in p.parse_pc_vars() that will be plotted
           e.g. p.grd_var_3d('sigma')


   Optional Parameters
   --------------------
   res : float, *optional* [default = 0.1]
   	   grid resolution

   cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented `here <http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html>`_

   pitch : float, *optional* [default = 10]
   	   rotates the camera. see `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch>`_

   azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360). See `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view>`_

   distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame. See `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view>`_

   log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

   smooth : bool, *optional* [default = True]
   	   if True, will smooth plotted dependent variable
           with a median filter and window size specified by filtsz (below)

   filtsz : int, *optional* [default = 3]
   	   size of filter (pixels) if smooth==1

   xsize : int, *optional* [default = 2000]
   	   size (number of pixels) of output image in x dimension

   ysize : int, *optional* [default = 1000]
   	   size (number of pixels) of output image in y dimension
          
p.grd_vars_3d() 
----------------
   produces a 3d plot of the gridded surface made from each output variable in p.parse_pc_vars()

   Syntax
   ----------
   [] = p.grd_vars_3d(res, cmap, pitch, azimuth, distance, log_scale, smooth, filtsz, xsize, ysize)

   Parameters
   ------------
   p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

   Optional Parameters
   --------------------
   res : float, *optional* [default = 0.1]
   	   grid resolution

   cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented `here <http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html>`_

   pitch : float, *optional* [default = 10]
   	   rotates the camera. see `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch>`_

   azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360). See `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view>`_

   distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame. See `here <http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view>`_

   log_scale : bool, *optional* [default = False]
   	   if True, will log scale plotted dependent variable

   smooth : bool, *optional* [default = True]
   	   if True, will smooth plotted dependent variable
           with a median filter and window size specified by filtsz (below)

   filtsz : int, *optional* [default = 3]
   	   size of filter (pixels) if smooth==1

   xsize : int, *optional* [default = 2000]
   	   size (number of pixels) of output image in x dimension

   ysize : int, *optional* [default = 1000]
   	   size (number of pixels) of output image in y dimension


  .. image:: _static/pysesa_colour.jpg


