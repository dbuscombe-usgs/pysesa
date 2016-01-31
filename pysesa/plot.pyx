## PySESA (Python program for Spatially Explicit Spectral Analysis) 
## has been developed at the Grand Canyon Monitorinf & Research Center,
## U.S. Geological Survey
##
## Author: Daniel Buscombe
## Project homepage: <https://github.com/dbuscombe-usgs/pysesa>
##
##This software is in the public domain because it contains materials that originally came from 
##the United States Geological Survey, an agency of the United States Department of Interior. 
##For more information, see the official USGS copyright policy at 
##http://www.usgs.gov/visual-id/credit_usgs.html#copyright
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

##"""
## ___      ___ ___ ___   _     _   _ 
##| _ \_  _/ __| __/ __| /_\   (_) (_)
##|  _/ || \__ \ _|\__ \/ _ \   _   _ 
##|_|  \_, |___/___|___/_/ \_\ (_) (_)
##     |__/                           
##           __      __ 
##    ____  / /___  / /_
##   / __ \/ / __ \/ __/
##  / /_/ / / /_/ / /_  
## / .___/_/\____/\__/  
##/_/  
##
###+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
##|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
##+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
##|d|b|u|s|c|o|m|b|e|@|u|s|g|s|.|g|o|v|
##+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
##+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
##|U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
##+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+

##"""


# =========================================================
# =============== import libraries ======================
# =========================================================
from __future__ import division
import os, sys, getopt
import matplotlib.pylab as plt
import numpy as np
from Tkinter import Tk
from tkFileDialog import askopenfilename

# =========================================================
# ===================== class ======================
# =========================================================

# =========================================================
class plot:
   '''
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
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
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
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
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
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
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
           marker size in points^2
           http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter

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
           marker size in points^2
           http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter

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
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

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
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

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
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

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

   '''

   def __init__(self, pysesa_file=None): 
      
      '''
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
      	   produces 2d plot of the gridded [x,y,z] surface made from decimated point cloud,
           as returned by parse_pc_vars()
               
      p.grd_var()   
   	   produces 2d plot of the gridded surface made from 1 output variable in p.parse_pc_vars()
           e.g. p.grd_var('sigma')
   
      p.grd_vars()
   	   produces a 2d plot of the gridded surface made from each output variable in p.parse_pc_vars()


      3D plotting functions (mayavi not required)
      --------------------------------------------
      p.plt_xyz() 
   	   produces 3d plot of Nx3 contents of raw point cloud, as returned by p.get_xyz()
       
      p.plt_xy_var()   
   	   produces 3d plot of 1 output variable in p.parse_pc_vars(), e.g. p.grd_var('sigma')
         
      p.plt_xy_vars()
   	   produces a 3d plot of each output variable in p.parse_pc_vars()      
          

      3D plotting functions (requires mayavi) 
      ----------------------------------------
      p.grd_xyz3d()   
   	   produces 3d plot of the gridded surface made from the Nx3 contents of raw point cloud,
           as returned by p.get_xyz()
   
      p.grd_var_3d()  
   	   produces 3d plot of the gridded surface made from 1 output variable in p.parse_pc_vars()
           e.g. p.grd_var('sigma')
          
      p.grd_vars_3d() 
   	   produces a 3d plot of the gridded surface made from each output variable in p.parse_pc_vars()

      '''

      if pysesa_file is None:
         Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
         pysesa_file = askopenfilename(filetypes=[("Pysesa output file","*.xyz")], multiple=True) 

      if type(pysesa_file) is str:
         pysesa_file = (pysesa_file,)
         
      self._outfiles = pysesa_file   

      self._proctype = int(pysesa_file[0].split('proctype')[1].split('_')[0])

      self._plt = plt

   # =========================================================
   # ===================== 3d plotting functions ======================
   # =========================================================

   #==================================================        
   def grd_vars_3d(self, res=0.1, thresdist = 0, nn=1, cmap='hot', pitch=10, azimuth=-200, distance=50, log_scale=False, smooth=True, filtsz=3, xsize=2000, ysize=1000):
      '''
      produces a 3d plot of the gridded surface made from each output variable in p.parse_pc_vars()

      Syntax
      ----------
      [] = p.grd_vars_3d(res, thresdist, nn, cmap, pitch, azimuth, distance, log_scale, smooth, filtsz, xsize, ysize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution
   	   
      thresdist : float, *optional* [default = 0]
   	   maximum interpolation distance 
   	   if 0, thresdist calculated as sqrt(1/res)
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding 
   	   
      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

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

      '''

      try:
         from mayavi import mlab

         pcdict = self.parse_pc_vars()

         x = pcdict['x']
         y = pcdict['y']
         z = pcdict['z_mean']

         dat, grid_x, grid_y = self._get_grid(x,y,z,res,thresdist,nn)

         for key in pcdict.keys():
            if key not in ['x','y','z_mean']:
               print "printing %s" % (key)

               if log_scale is True:
                  datz = self._grid_var(x, y, np.log(pcdict[key]), grid_x, grid_y, res)
               else:
                  datz = self._grid_var(x, y, pcdict[key], grid_x, grid_y, res)

               if smooth is True:
                  from scipy.ndimage.filters import median_filter
                  datz = median_filter(datz,(filtsz,filtsz))

               fig = self._get_3d_fig()
               surf = mlab.mesh(grid_x,grid_y,dat,scalars=datz,colormap=cmap) 
               mlab.pitch(pitch)
               mlab.view(azimuth=azimuth)
               mlab.view(distance=distance)

               mlab.savefig(self._outfiles[0].split('.')[-1-1]+'grd_pc_'+str(key)+'_3d.png',figure=fig,size=(xsize,ysize))              
               mlab.close()
               del fig
      except:
         print "error: mayavi is required for 3d plots"


   #==================================================        
   def grd_var_3d(self, var, res=0.1, thresdist = 0, nn=1, cmap='hot', pitch=10, azimuth=-200, distance=50, log_scale=False, smooth=True, filtsz=3, xsize=2000, ysize=1000):
      '''
      produces 3d plot of the gridded surface made from 1 output variable in p.parse_pc_vars()

      Syntax
      ----------
      [] = p.grd_var_3d(var, res, thresdist, nn, cmap, pitch, azimuth, distance, log_scale, smooth, filtsz, xsize, ysize)

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

      thresdist : float, *optional* [default = 0]
   	   maximum interpolation distance 
   	   if 0, thresdist calculated as sqrt(1/res)
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding 
   	      	   
      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

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

      '''

      try:
         from mayavi import mlab

         pcdict = self.parse_pc_vars()

         x = pcdict['x']
         y = pcdict['y']
         z = pcdict['z_mean']

         dat, grid_x, grid_y = self._get_grid(x,y,z,res,thresdist,nn)

         for key in pcdict.keys():
            if key in [var]:
               print "printing %s" % (key)

               if log_scale is True:
                  datz = self._grid_var(x, y, np.log(pcdict[key]), grid_x, grid_y, res)
               else:
                  datz = self._grid_var(x, y, pcdict[key], grid_x, grid_y, res)

               if smooth is True:
                  from scipy.ndimage.filters import median_filter
                  datz = median_filter(datz,(filtsz,filtsz))

               fig = self._get_3d_fig()
               surf = mlab.mesh(grid_x,grid_y,dat,scalars=datz,colormap=cmap) 
               mlab.pitch(pitch)
               mlab.view(azimuth=azimuth)
               mlab.view(distance=distance)

               mlab.savefig(self._outfiles[0].split('.')[-1-1]+'grd_pc_'+str(key)+'_3d.png',figure=fig,size=(xsize,ysize))              
               mlab.close()
               del fig
      except:
         print "error: mayavi is required for 3d plots"


   #==================================================        
   def grd_xyz3d(self, res=0.1, thresdist = 0, nn=1, cmap='hot', pitch=10, azimuth=-200, distance=50, xsize=2000, ysize=1000):
      '''
      produces 3d plot of the gridded surface made from the Nx3 contents of raw point cloud,
      as returned by p.get_xyz()

      Syntax
      ----------
      [] = p.grd_xyz3d(res, thresdist, nn, cmap, pitch, azimuth, distance, xsize, ysize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()


      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      thresdist : float, *optional* [default = 0]
   	   maximum interpolation distance 
   	   if 0, thresdist calculated as sqrt(1/res)
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding 
   	      	   
      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html

      pitch : float, *optional* [default = 10]
   	   rotates the camera. see:
           http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=pitch#mayavi.mlab.pitch

      azimuth : float, *optional* [default = -200]
           The azimuthal angle (in degrees, 0-360)
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      distance : float or 'auto', *optional* [default = 50]
           A positive floating point number representing the distance from the focal point to place the camera.
           if ‘auto’ is passed, the distance is computed to have a best fit of objects in the frame.
   	   http://docs.enthought.com/mayavi/mayavi/auto/mlab_camera.html?highlight=view#mayavi.mlab.view

      xsize : int, *optional* [default = 2000]
   	   size (number of pixels) of output image in x dimension

      ysize : int, *optional* [default = 1000]
   	   size (number of pixels) of output image in y dimension

      '''

      try:
         from mayavi import mlab
         pcdict = self.parse_pc_vars()

         x = pcdict['x']
         y = pcdict['y']
         z = pcdict['z_mean']

         dat, grid_x, grid_y = self._get_grid(x,y,z,res,thresdist,nn)

         fig = self._get_3d_fig()
         surf = mlab.mesh(grid_x,grid_y,dat,scalars=dat,colormap=cmap) 
         mlab.pitch(pitch)
         mlab.view(azimuth=azimuth)
         mlab.view(distance=distance)

         mlab.savefig(self._outfiles[0].split('.')[-1-1]+'grd_pc_3d.png',figure=fig,size=(xsize,ysize))        
         mlab.close()
         del fig
      except:
         print "error: mayavi is required for 3d plots"


   # =========================================================
   # ===================== 2d plotting functions ======================
   # =========================================================

   #==================================================        
   def grd_var(self, var, res=0.1, thresdist = 0, nn=1, azimuth=315, altitude=45, zf=1, cmap='hot', dpi=300, log_scale=False, smooth=True, filtsz=3, alpha=0.5, ticksize=4, labelsize=6):
      '''
      produces 2d plot of the gridded surface made from 1 output variable in p.parse_pc_vars()
      e.g. p.grd_var('sigma')

      Syntax
      ----------
      [] = p.grd_var(var, res, thresdist, nn, azimuth, altitude, zf, cmap, dpi, log_scale, smooth, filtsz, alpha, ticksize, labelsize)

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

      thresdist : float, *optional* [default = 0]
   	   maximum interpolation distance 
   	   if 0, thresdist calculated as sqrt(1/res)
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding 
   	      	   
      azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
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

      '''

      pcdict = self.parse_pc_vars()

      x = pcdict['x']
      y = pcdict['y']
      z = pcdict['z_mean']

      dat, grid_x, grid_y = self._get_grid(x,y,z,res,thresdist,nn)

      h = self._hillshade(dat,res,res,azimuth,altitude,zf)
      del dat, z

      for key in pcdict.keys():
         if key in [var]:
            print "printing %s" % (key)
            if log_scale is True:
               dat = self._grid_var(x, y, np.log(pcdict[key]), grid_x, grid_y, res)
            else:
               dat = self._grid_var(x, y, pcdict[key], grid_x, grid_y, res)

            if smooth is True:
               from scipy.ndimage.filters import median_filter
               dat = median_filter(dat,(filtsz,filtsz))

            fig = self._plt.figure()
            ax = fig.add_subplot(2,2,1)   
            self._plt.pcolor(grid_x-np.min(grid_x), grid_y-np.min(grid_y), np.ma.masked_invalid(h), cmap='gray');
            im=self._plt.pcolormesh(grid_x-np.min(grid_x), grid_y-np.min(grid_y), np.ma.masked_invalid(dat), cmap=cmap, alpha=alpha);    

            ax = self._set_tick_size(ax, ticksize)
            ax = self._make_xy_axislabels(ax, labelsize)

            self._divide_colorbar(ax, im, str(key))

            self._plt.savefig(self._outfiles[0].split('.')[-1-1]+'grd_pc_'+key,bbox_inches='tight',dpi=dpi)            
            del fig
            self._plt.close()

   #==================================================        
   def grd_vars(self, res=0.1, thresdist = 0, nn=1, azimuth=315, altitude=45, zf=1, cmap='hot', dpi=300, log_scale=False, smooth=True, filtsz=3, alpha=0.5, ticksize=4, labelsize=6):
      '''
      produces a 2d plot of the gridded surface made from each output variable in p.parse_pc_vars()

      Syntax
      ----------
      [] = p.grd_vars(res, thresdist, nn, azimuth, altitude, zf, cmap, dpi, log_scale, smooth, filtsz, alpha, ticksize, labelsize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      thresdist : float, *optional* [default = 0]
   	   maximum interpolation distance 
   	   if 0, thresdist calculated as sqrt(1/res)
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding 
   	      	   
      azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
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

      '''

      pcdict = self.parse_pc_vars()

      x = pcdict['x']
      y = pcdict['y']
      z = pcdict['z_mean']

      dat, grid_x, grid_y = self._get_grid(x,y,z,res,thresdist,nn)

      h = self._hillshade(dat,res,res,azimuth,altitude,zf)
      del dat, z

      for key in pcdict.keys():
         if key not in ['x','y','z_mean']:
            print "printing %s" % (key)
            if log_scale is True:
               dat = self._grid_var(x, y, np.log(pcdict[key]), grid_x, grid_y, res, thresdist)
            else:
               dat = self._grid_var(x, y, pcdict[key], grid_x, grid_y, res, thresdist)
               
            if smooth is True:
               from scipy.ndimage.filters import median_filter
               dat = median_filter(dat,(filtsz,filtsz))

            fig = self._plt.figure()
            ax = fig.add_subplot(2,2,1)   
            self._plt.pcolor(grid_x-np.min(grid_x), grid_y-np.min(grid_y), np.ma.masked_invalid(h), cmap='gray');
            im=self._plt.pcolormesh(grid_x-np.min(grid_x), grid_y-np.min(grid_y), np.ma.masked_invalid(dat), cmap=cmap, alpha=alpha);    

            ax = self._set_tick_size(ax, ticksize)
            ax = self._make_xy_axislabels(ax, labelsize)

            self._divide_colorbar(ax, im, str(key))

            self._plt.savefig(self._outfiles[0].split('.')[-1-1]+'grd_pc_'+key,bbox_inches='tight',dpi=dpi)            
            del fig
            self._plt.close()


   #==================================================        
   def grd_xyz(self, res=0.1, thresdist = 0, nn = 1, azimuth=315, altitude=45, zf=1, cmap='hot', dpi=300, alpha=0.5, ticksize=4, labelsize=6):
      '''
      produces 2d plot of the gridded [x,y,z] surface made from decimated point cloud,
      as returned by parse_pc_vars()

      Syntax
      ----------
      [] = p.grd_xyz(res, thresdist, nn, azimuth, altitude, zf, cmap, dpi, alpha, ticksize, labelsize)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      thresdist : float, *optional* [default = 0]
   	   maximum interpolation distance 
   	   if 0, thresdist calculated as sqrt(1/res)
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding 
   	      	   
      azimuth : float, *optional* [default = 315]
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float, *optional* [default = 45]
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float, *optional* [default = 1]
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      cmap : str, *optional* [default = 'hot']
   	   colormap
           possible colormaps are documented here:
           http://matplotlib.org/examples/color/colormaps_reference.html
   	   
      dpi : int, *optional* [default = 300]
           figure resolution in dots per inch

      alpha : float, *optional* [default = 0.5]
           transparency, between 0.0 and 1.0

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels

      '''
      pcdict = self.parse_pc_vars()

      x = pcdict['x']
      y = pcdict['y']
      z = pcdict['z_mean']

      dat, grid_x, grid_y = self._get_grid(x,y,z,res,thresdist,nn)

      h = self._hillshade(dat,res,res,azimuth,altitude,zf)

      fig = self._plt.figure()
      ax = fig.add_subplot(2,2,1)   
      self._plt.pcolor(grid_x-np.min(grid_x), grid_y-np.min(grid_y), np.ma.masked_invalid(h), cmap='gray');
      im=self._plt.contourf(grid_x-np.min(grid_x), grid_y-np.min(grid_y), np.ma.masked_invalid(dat), cmap=cmap, alpha=alpha);    

      ax = self._set_tick_size(ax, ticksize)
      ax = self._make_xy_axislabels(ax, labelsize)

      self._divide_colorbar(ax, im, r"Amplitude")

      plt.savefig(self._outfiles[0].split('.')[-1-1]+'grd_pc',bbox_inches='tight',dpi=dpi)
      del fig
      self._plt.close()
      
   #==================================================        
   def plt_xyz(self,elev=65, azim=-115, markersize=.01, dpi=300, ticksize=4, labelsize=6):
      '''
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

      '''
      from mpl_toolkits.mplot3d import Axes3D

      xyz = self.get_xyz()

      fig = self._plt.figure(dpi=dpi)
      fig.subplots_adjust(wspace = -0.15, hspace=-0.15)
      ax = fig.add_subplot(221, projection='3d', axisbg='white')
      ax.view_init(elev=elev, azim=azim)
      ax.plot(xyz[:,0]-np.min(xyz[:,0]),xyz[:,1]-np.min(xyz[:,1]),xyz[:,2]-np.min(xyz[:,2]), linestyle="None", color='k',  marker=".", markersize=markersize)
      ax.grid(False)
      self._plt.axis('tight')

      ax = self._set_tick_size(ax, ticksize)
      ax = self._set_ztick_size(ax, ticksize)
      ax = self._make_axislabels(ax, labelsize)

      self._plt.savefig(self._outfiles[0].split('.')[-1-1]+'xyp_3d_pc',bbox_inches='tight',dpi=dpi)
      del fig
      self._plt.close()

   #==================================================        
   def plt_xy_var(self, var, log_scale=True, dpi=300, markersize=5, ticksize=4, labelsize=6, elev=65, azim=-115):
      '''
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
           marker size in points^2
           http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels

      elev : float, *optional* [default = 65]
           the elevation angle in the z plane

      azimuth : float, *optional* [default = -115]
           azimuth angle in the x,y plane 

      '''
      from mpl_toolkits.mplot3d import Axes3D

      pcdict = self.parse_pc_vars()

      x = pcdict['x']
      y = pcdict['y']
      z = pcdict['z_mean']
      for key in pcdict.keys():
         if key in [var]:
            print "printing %s" % (key)
            fig = self._plt.figure(dpi=dpi)
            fig.subplots_adjust(wspace = -0.15, hspace=-0.15)
            ax = fig.add_subplot(221, projection='3d', axisbg='white')
            ax.view_init(elev=elev, azim=azim)
            if log_scale is True:
               im=ax.scatter(x-np.min(x),y-np.min(y), z, s=markersize, c=pcdict[key], linewidth=0)
            else:
               im=ax.scatter(x-np.min(x),y-np.min(y), z, s=markersize, c=np.log(pcdict[key]), linewidth=0)
            ax.grid(False)
            self._plt.axis('tight')
            self._plt.colorbar(im)

            ax = self._set_tick_size(ax, ticksize)
            ax = self._set_ztick_size(ax, ticksize)
            ax = self._make_axislabels(ax, labelsize)

            self._plt.savefig(self._outfiles[0].split('.')[-1-1]+'xyp_3d_'+key,bbox_inches='tight',dpi=dpi)            
            del fig
            self._plt.close()


   #==================================================        
   def plt_xy_vars(self, log_scale=True, dpi=300, markersize=5, ticksize=4, labelsize=6, elev=65, azim=-115):
      '''
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
           marker size in points^2
           http://matplotlib.org/mpl_toolkits/mplot3d/api.html#mpl_toolkits.mplot3d.axes3d.Axes3D.scatter

      ticksize : int, *optional* [default = 4]
           size of x, y, and z tick labels

      labelsize : int, *optional* [default = 6]
           size of x and y axes labels

      elev : float, *optional* [default = 65]
           the elevation angle in the z plane

      azimuth : float, *optional* [default = -115]
           azimuth angle in the x,y plane 

      '''
      from mpl_toolkits.mplot3d import Axes3D

      pcdict = self.parse_pc_vars()

      x = pcdict['x']
      y = pcdict['y']
      z = pcdict['z_mean']
      for key in pcdict.keys():
         if key not in ['x','y','z_mean']:
            print "printing %s" % (key)
            fig = self._plt.figure(dpi=dpi)
            fig.subplots_adjust(wspace = -0.15, hspace=-0.15)
            ax = fig.add_subplot(221, projection='3d', axisbg='white')
            ax.view_init(elev=elev, azim=azim)
            if log_scale is True:
               im=ax.scatter(x-np.min(x),y-np.min(y), z, s=markersize, c=pcdict[key], linewidth=0)
            else:
               im=ax.scatter(x-np.min(x),y-np.min(y), z, s=markersize, c=np.log(pcdict[key]), linewidth=0)
            ax.grid(False)
            self._plt.axis('tight')
            self._plt.colorbar(im)

            ax = self._set_tick_size(ax, ticksize)
            ax = self._set_ztick_size(ax, ticksize)
            ax = self._make_axislabels(ax, labelsize)

            self._plt.savefig(self._outfiles[0].split('.')[-1-1]+'xyp_3d_'+key,bbox_inches='tight',dpi=dpi)            
            del fig
            self._plt.close()


   # =========================================================
   # ===================== data reading functions ======================
   # =========================================================

   #==================================================        
   def get_pc(self):
      '''
      parse an NxM contents of pysesa_file

      Syntax
      ----------
      pc = p.get_pc()

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Returns
      ----------
      pc : ndarray
   	   list of M ndarrays, each containing a column
           output from the pysesa_file
      '''
      pc = []
      for outfile in self._outfiles:
         pc.append(np.genfromtxt(outfile, delimiter=','))

      pc = np.vstack(pc).astype('float32')
      return pc

   #==================================================        
   def get_xyz(self):
      '''
      parse the Nx3 contents of raw point cloud (the file processed by pysesa::process)

      Syntax
      ----------
      xyz = p.get_xyz()

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Returns
      ----------
      xyz : ndarray
   	   Nx3 array point cloud input
      '''
      xyz = []
      for outfile in self._outfiles:
         infile = outfile.split('_zstat')[0]
         xyz.append(self._txtread(infile))

      xyz = np.vstack(xyz).astype('float32')
      return xyz

   #==================================================        
   def _txtread(self,infile):
      '''
      anonymous fast(er than numpy's genfromtxt) function to read a Nx3 txt file to numpy array

      Syntax
      ----------
      xyz = p._txtread()

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Returns
      ----------
      xyz : ndarray
   	   Nx3 array point cloud input
      '''
      f = open(infile, 'rb'); data = f.read(); f.close()
      data = data.splitlines()
      if len(data[0].split(','))>1: # the file is comma delimited
        return np.array([x.split(',',2)[0:3] for x in data], dtype='float32')
      else: 
         try: # space delimited
            return np.array([x.split(' ',2)[0:3] for x in data], dtype='float32')
         except: #tab delimited
            return np.array([x.split('\t',2)[0:3] for x in data], dtype='float32')

   #==================================================        
   def parse_pc_vars(self):
      '''
      NxM contents of pysesa_file parsed into dict object
      1 key per variable in p.get_pc()
      vars.keys() returns list of variables in dict

      Syntax
      ----------
      [] = p.parse_pc_vars()

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Returns
      ----------
      pcdict : dict
   	   dict object of M, Nx1 arrays point cloud outputs
      '''
      pc = self.get_pc()
      pcdict = {}

      ## proctype = 1 or 2
      if (self._proctype==1) or (self._proctype==2):
         pcdict['x'] = pc[:,0]
         pcdict['y'] = pc[:,1]
         pcdict['z_mean'] = pc[:,2]
         pcdict['slope'] = pc[:,3]
         pcdict['intercept'] = pc[:,4]
         pcdict['r_value'] = pc[:,5]
         pcdict['p_value'] = pc[:,6]
         pcdict['std_err'] = pc[:,7]
         pcdict['d'] = pc[:,8]
         pcdict['l'] = pc[:,9]
         pcdict['wmax'] = pc[:,10]
         pcdict['wmean'] = pc[:,11]
         pcdict['rms1'] = pc[:,12]
         pcdict['rms2'] = pc[:,13]
         pcdict['Z'] = pc[:,14]
         pcdict['E'] = pc[:,15]
         pcdict['sigma'] = pc[:,16] 
         pcdict['T0_1'] = pc[:,17]
         pcdict['T0_2'] = pc[:,18] 
         pcdict['sw1'] = pc[:,19]
         pcdict['sw2'] = pc[:,20]
         pcdict['m0'] = pc[:,21]
         pcdict['m1'] = pc[:,22] 
         pcdict['m2'] = pc[:,23]
         pcdict['m3'] = pc[:,24]
         pcdict['m4'] = pc[:,25]
         pcdict['phi'] = pc[:,26]

      # proctype = 3
      if (self._proctype==3):
         pcdict['x'] = pc[:,0]
         pcdict['y'] = pc[:,1]
         pcdict['z_mean'] = pc[:,2]
         pcdict['z_max'] = pc[:,3]
         pcdict['z_min'] = pc[:,4]
         pcdict['z_range'] = pc[:,5]
         pcdict['sigma'] = pc[:,6]
         pcdict['skewness'] = pc[:,7]
         pcdict['kurtosis'] = pc[:,8]
         pcdict['n'] = pc[:,9]

      # proctype = 4 or 5
      if (self._proctype==4) or (self._proctype==5):
         pcdict['x'] = pc[:,0]
         pcdict['y'] = pc[:,1]
         pcdict['z_mean'] = pc[:,2]
         pcdict['z_max'] = pc[:,3]
         pcdict['z_min'] = pc[:,4]
         pcdict['z_range'] = pc[:,5]
         pcdict['sigma'] = pc[:,6]
         pcdict['skewness'] = pc[:,7]
         pcdict['kurtosis'] = pc[:,8]
         pcdict['n'] =pc[:,9]
         pcdict['slope'] = pc[:,10]
         pcdict['intercept'] = pc[:,11]
         pcdict['r_value'] = pc[:,12]
         pcdict['p_value'] = pc[:,13]
         pcdict['std_err'] = pc[:,14]
         pcdict['d'] = pc[:,15]
         pcdict['l'] = pc[:,16]
         pcdict['wmax'] = pc[:,17]
         pcdict['wmean'] = pc[:,18]
         pcdict['rms1'] = pc[:,19]
         pcdict['rms2'] = pc[:,20]
         pcdict['Z'] = pc[:,21]
         pcdict['E'] = pc[:,22]
         pcdict['sigma'] = pc[:,23]
         pcdict['T0_1'] = pc[:,24]
         pcdict['T0_2'] = pc[:,25]
         pcdict['sw1'] = pc[:,26]
         pcdict['sw2'] = pc[:,27]
         pcdict['m0'] = pc[:,28]
         pcdict['m1'] = pc[:,29]
         pcdict['m2'] = pc[:,30]
         pcdict['m3'] = pc[:,31]
         pcdict['m4'] = pc[:,32]
         pcdict['phi'] = pc[:,33]

      return pcdict


   # =========================================================
   # ===================== numerical functions ======================
   # =========================================================

   #==================================================        
   def _grid_var(self, x, y, z, grid_x, grid_y, res=0.1, thresdist=0, nn=1):
      '''
      anonymous function to produce a gridded surface

      Syntax
      ----------
      [] = p._grid_var(x, y, z, grid_x, grid_y, res, thresdist, nn)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      x : ndarray
   	   array of x values

      y : ndarray
   	   array of y values

      z : ndarray
   	   array of amplitude values

      grid_x : ndarray
   	   mesh of [x,y]

      grid_y : ndarray
   	   mesh of [x,y].T

      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      thresdist : float, *optional* [default = 0]
   	   maximum interpolation distance 
   	   if 0, thresdist calculated as sqrt(1/res)
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding    
   	   	   
      '''
      #from scipy.interpolate import griddata

      #dat = griddata(np.c_[x,y], z, (grid_x, grid_y), method='nearest')

      dat, dist, res = self._get_dist(x, y, z, grid_x, grid_y, res, nn)
      ## mask
      if thresdist == 0:
         thresdist = np.sqrt(1/res)     
 
      dat[dist>thresdist] = np.nan
      #dat[dist>(np.sqrt(1/res)) ] = np.nan
      dat[dat<np.min(z)] = np.nan

      return dat


   # =========================================================
   def _getmesh(self, minX, maxX, minY, maxY, res):

      complete=0
      while complete==0:
         try:
            grid_x, grid_y = np.meshgrid( np.arange(minX, maxX, res), np.arange(minY, maxY, res) )
            if 'grid_x' in locals(): 
               complete=1 
         except:
            print "memory error: trying grid resolution of %s" % (str(res*2))
            res = res*2
         
      return grid_x, grid_y, res


   #==================================================        
   def _get_grid(self, x, y, z, res=0.1, thresdist = 0, nn=1):
      '''
      anonymous function to produce a gridded surface

      Syntax
      ----------
      [] = p._get_grid(x, y, z, res, nn)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      x : ndarray
   	   array of x values

      y : ndarray
   	   array of y values

      z : ndarray
   	   array of amplitude values

      Optional Parameters
      --------------------
      res : float, *optional* [default = 0.1]
   	   grid resolution

      thresdist : float, *optional* [default = 0]
   	   maximum interpolation distance 
   	   if 0, thresdist calculated as sqrt(1/res)
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding 
   	   
      '''
      #from scipy.interpolate import griddata
      
      # create grid for gridding
      #grid_x, grid_y = np.meshgrid( np.arange(np.min(x), np.max(x), res), np.arange(np.min(y), np.max(y), res) )  
      grid_x, grid_y, res = self._getmesh(np.min(x), np.max(x), np.min(y), np.max(y), res)
      #xyz = self.get_xyz()

      # create grid
      #dat = griddata(np.c_[xyz[:,0],xyz[:,1]], xyz[:,2], (grid_x, grid_y), method='nearest')
      
      dat, dist, res = self._get_dist(x, y, z, grid_x, grid_y, res, nn)

      ## mask
      if thresdist == 0:
         thresdist = np.sqrt(1/res)
      
      dat[dist>thresdist] = np.nan
      #dat[dist>(np.sqrt(1/res)) ] = np.nan
      #dat[dat<np.min(xyz[:,2])] = np.nan
      return dat, grid_x, grid_y

   # =========================================================
   def _hillshade(self, dem, dx, dy, azimuth, altitude, zf):
      '''
      anonymous function to produce a shaded relief of a gridded surface
      using the ESRI algorithm

      Syntax
      ----------
      [] = p._hillshade(dem, dx, dy, azimuth, altitude, zf)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      dem : ndarray
   	   2d array of gridded surface values

      dx : float
   	   difference in x dimension

      dy : float
   	   difference in y dimension

      azimuth : float
           Lighting azimuthal angle (in degrees, 0-360)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      altitude : float
           Lighting zenith angle (in degrees, 0-90)
           for hillshade calculation, see:
           http://edndoc.esri.com/arcobjects/9.2/net/shared/geoprocessing/spatial_analyst_tools/how_hillshade_works.htm

      zf : float
           Vertical exaggeration factor
           1=no exaggeration, <1 minimizes, >1 exaggerates 

      '''
      # lighting azimuth
      azimuth = 360.0-azimuth+90 #convert to mathematic unit
      if azimuth>360 or azimuth==360:
         azimuth = azimuth-360
      azimuth = azimuth * (np.pi/180) #  convert to radians

      # lighting altitude
      altitude = (90-altitude) * (np.pi/180) # convert to zenith angle in radians

      # calc slope and aspect (radians)
      fx,fy = np.gradient(dem,dx) # uses simple, unweighted gradient of immediate $
      [asp,grad] = self._cart2pol(fy,fx) # % convert to carthesian coordinates

      grad = np.arctan(zf*grad) #steepest slope
      # convert asp
      asp[asp<np.pi] = asp[asp<np.pi]+(np.pi/2)
      asp[asp<0] = asp[asp<0]+(2*np.pi)

      ## hillshade calculation
      h = 255.0*( (np.cos(altitude)*np.cos(grad) ) + ( np.sin(altitude)*np.sin(grad)*np.cos(azimuth - asp) ))
      h[h<0]=0 # % set hillshade values to min of 0.

      return h

   # =========================================================
   # ===================== auxiliary numerical functions ======================
   # =========================================================

   #==================================================        
   def _get_dist(self, x, y, z, grid_x, grid_y, res, nn):
      '''
      anonymous function to create a nearest neighbour distance matrix for dem masking

      Syntax
      ----------
      [] = p._get_dist(x, y, z, grid_x, grid_y, res, nn)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      x : ndarray
   	   array of x values

      y : ndarray
   	   array of y values

      z : ndarray
   	   array of z value
   	   
      grid_x : ndarray
   	   mesh of [x,y]

      grid_y : ndarray
   	   mesh of [x,y].T 
   	   
      res : float
   	   grid resolution  
   	   
      nn : int, *optional* [default = 1]
   	   number of nearest neighbours used in the gridding 
      '''
      #from scipy.spatial import cKDTree as KDTree
      try:
         from pykdtree.kdtree import KDTree
         pykdtree=1   
      except:
         print "install pykdtree for faster kd-tree operations: https://github.com/storpipfugl/pykdtree"
         from scipy.spatial import cKDTree as KDTree
         pykdtree=0   

      ## create mask for where the data is not
      tree = KDTree(np.c_[x,y])         
            
      complete=0
      while complete==0:
         try:
            dist, inds = tree.query(np.c_[grid_x.ravel(), grid_y.ravel()], k=nn)
            
            if nn==1:
               dat = z.flatten()[inds].reshape(np.shape(grid_x))    
            else:
               w = 1.0 / dist**2
               dat = np.sum(w * z.flatten()[inds,2], axis=1) / np.sum(w, axis=1)
               dat.shape = grid_x.shape
            if 'dat' in locals(): 
               complete=1 
         except:  
            grid_x, grid_y, res = self._getmesh(np.min(x), np.max(x), np.min(y), np.max(y), res*2)          
         
      return dat, dist.reshape(grid_x.shape), res

   # =========================================================
   def _cart2pol(self, x, y):
      '''
      anonymous function to convert cartesian to polar coordinates

      Syntax
      ----------
      [] = p._cart2pol(x, y)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      x : ndarray
   	   array of x values

      y : ndarray
   	   array of y values
      '''
      theta = np.arctan2(y, x)
      rho = np.sqrt(x**2 + y**2)
      return (theta, rho) 

   # =========================================================
   # ===================== auxiliary plotting functions ======================
   # =========================================================

   # =========================================================
   def _dolabels(self, rot=30): 
      '''
      anonymous function to create x and y labels with less density than default and rotate them

      Syntax
      ----------
      [] = p._dolabels(rot)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      rot : float, *optional* [default = 30]
   	   angle of rotation
      '''
      locs, labels = self._plt.xticks()
      self._plt.xticks(locs[::3])
      locs, labels = self._plt.xticks()
      self._plt.setp(labels, rotation=rot)
   
      locs, labels = self._plt.yticks()
      self._plt.yticks(locs[::3])
      locs, labels = self._plt.yticks()
      self._plt.setp(labels, rotation=rot)

   # =========================================================
   def _set_tick_size(self, ax, size=4):
      '''
      anonymous function to set tick size of x and y ticks

      Syntax
      ----------
      [] = p._set_tick_size(ax, size)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      ax : object
   	   matplotlib axes handle

      Optional Parameters
      --------------------
      size : float, *optional* [default = 4]
   	   font size (points)
      '''
      for tick in ax.xaxis.get_major_ticks():
         tick.label.set_fontsize(size)
      for tick in ax.yaxis.get_major_ticks():
         tick.label.set_fontsize(size)
      return ax

   # =========================================================
   def _set_ztick_size(self, ax, size=4):
      '''
      anonymous function to set tick size of z ticks

      Syntax
      ----------
      [] = p._set_ztick_size(ax, size)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      ax : object
   	   matplotlib axes handle

      Optional Parameters
      --------------------
      size : float, *optional* [default = 4]
   	   font size (points)
      '''
      for tick in ax.zaxis.get_major_ticks():
          tick.label.set_fontsize(size)
      return ax

   # =========================================================
   def _rmticks(self, ax):
      '''
      anonymous function to remove major and minor ticks

      Syntax
      ----------
      [] = p._rmticks(ax)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      ax : object
   	   matplotlib axes handle
      '''
      for tic in ax.xaxis.get_major_ticks():
         tic.tick1On = tic.tick2On = False
      for tic in ax.yaxis.get_major_ticks():
         tic.tick1On = tic.tick2On = False
      return ax

   # =========================================================
   def _make_axislabels(self, ax, size=6):
      '''
      anonymous function to label a generic 3d plot

      Syntax
      ----------
      [] = p._make_axislabels(ax, size)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      ax : object
   	   matplotlib axes handle

      Optional Parameters
      --------------------
      size : float, *optional* [default = 6]
   	   font size (points)
      '''
      ax.set_xlabel('Horizontal distance (length)', fontsize=size)
      ax.set_ylabel('Lateral distance (length)', fontsize=size)
      ax.set_zlabel('Vertical distance (length)',fontsize=size) 
      return ax

   # =========================================================
   def _make_xy_axislabels(self, ax, size=6):
      '''
      anonymous function to label x and y axes

      Syntax
      ----------
      [] = p._make_xy_axislabels(ax, size)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      ax : object
   	   matplotlib axes handle

      Optional Parameters
      --------------------
      size : float, *optional* [default = 6]
   	   font size (points)
      '''
      ax.set_xlabel('Horizontal distance (length)', fontsize=size)
      ax.set_ylabel('Lateral distance (length)', fontsize=size)
      return ax
      
   # =========================================================
   def _divide_colorbar(self, ax, im, label, size=6):
      '''
      anonymous function to make up a colorbar scaled with axes size

      Syntax
      ----------
      [] = p._divide_colorbar(ax, im, label, size)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      ax : object
   	   matplotlib axes handle

      im : object
   	   matplotlib image handle

      label : str
   	   colorbar ylabel string

      Optional Parameters
      --------------------
      size : float, *optional* [default = 6]
   	   font size (points)
      '''
      from mpl_toolkits.axes_grid1 import make_axes_locatable
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=-0.1)
      cbar = self._plt.colorbar(im, cax=cax)
      cbar.ax.set_ylabel(label, fontsize=size) 

      locs, labels = self._plt.yticks()
      self._plt.yticks(locs[::3])
      locs, labels = self._plt.yticks()
      self._plt.setp(labels, fontsize=size)  
      
   #==================================================        
   def _get_3d_fig(self, mag=5):
      '''
      anonymous function to create a 3d mayavi figure

      Syntax
      ----------
      [] = p._get_3d_fig(mag)

      Parameters
      ------------
      p : instance
   	   pysesa.plot instance returned by pysesa::plot
           e.g. p = pysesa.plot()

      Optional Parameters
      --------------------
      mag : float, *optional* [default = 5]
   	   the magnification is the scaling between the pixels on the screen, and the pixels in the file saved. 
           If you do not specify it, it will be calculated so that the file is saved with the specified size. 
           If you specify a magnification, Mayavi will use the given size as a screen size, 
           and the file size will be ‘magnification * size’.
      '''
      from mayavi import mlab
      fig = mlab.figure() 
      mlab.clf()
      fig.scene.set(magnification=mag)
      fig.scene.set(background=(0.5,0.5,0.5))
      return fig


