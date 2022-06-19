# koppenpasta
A collection of scripts for use with ExoPlaSim, with the intention of allowing the user to import arbitrary topography and then produce climate maps of the output

Includes:

koppenpasta.py
  Script to take data from NetCDF files produced by ExoPlaSim and produce maps of Koppen climate zones (or Holdridge Life Zones).
  Requres NetCDF4, PIL, numpy, and scipy
  
  As of 1.2, can be imported and run as the koppenpasta.MakeMap() function. Easiest way to do this is to create a config and then run it the form:
  
  koppenpasta.MakeMap('infile.nc', 'outfile', 'config.ini')
  
  where infile.nc is a NetCDF file, outfile the name of the output imag, and config.ini is the config. But if you desire it can take all the following optional arguments, most corresponding to the configuration options shown when run in konsole (y/n options should be input as 1/0 here):
  
  in_files: string showing the path to a NetCDF file or folder containing NetCDF files, all of which will be used, or a list of such strings.
  
  output_name: name for the output image
  
  cfg_loadname: name for a config .ini to be used (will overwrite all other options)
  
  land_type: Land Climate Zones (0-5)
  
  sea_type: Sea Climate Zones (0-5)
  
  color_type: Color List (0-3)
  
  col_list_path: path to custom color list
  
  blend: Land/Sea Blend (0-1)
  
  bin_num: Bin Size (should be at least 1)
  
  interp: Interpolation Factor (0 for none)
  
  dum_ice: Dummy Sea Ice (1/0)
  
  use_topo: Adjust Temp by Topography (1/0)
  
  topo_path: path to topography map to be used for above
  
  maxel: Maximum elevation of topography map (m)
  
  minel: Minimum elevation of topography map (m)
  
  gravity: Surface gravity for use with topography map (m/s^2)
  
  blend_topo: Blend by Topography Map (1/0)
  
  sealev: Sea Level (m)
  
  sum_def: Summer Definition (0-1)
  
  temp_def: Temperate/Continental Boundary (0-1)
  
  arid_def: Hot/Cold Arid Boundary (0-1)
  
  med_def: Mediterranean Definition (0-1)
  
  wet_def: Wet-Summer Definition (0-1)
  
  add_def: Additional Wet Season Requirements (0-3)
  
  prio_def: Med/Wet-Summer Priority (0-1)
  
  ice_def: Arid/Polar Priority (0-2)
  
  sea_def: Sea Ice Definition (0-1)
  

defaultcolor.ini
  config file to allow user to select their own rgb values for the maps output by the above script

Image2sra_2.0.3.py
  Script to take greyscale heightmap images and convert them to scaled .sra files that can be read by ExoPlaSim as topography input. Created by Alex (Ostimeus), distributed with permission.
  Ost has made their own repository for this script bundled into an exoplasim configuration tool https://github.com/MegalexMaster/ExoPlaSim-InCon
  I'll leave this script in place here but it may not include later bugfixes (and the tool is nifty anyway).
  
eps_avg.py
  Script to average together multiple NetCDF outputs from ExoPlaSim together, find annual averages, shift longitudes to show better in Panoply, and offset results by a given number of months. Requires NCO (http://nco.sourceforge.net/); get "nco" and dependencies from yast and run "pip install nco" in konsole (both; first installs binaries, second installs python API)
  
  Can also be imported and run as eps_avg.main(), with the arguments:
  
  in_files: string showing the path to a NetCDF file or folder containing NetCDF files, all of which will be used, or a list of such strings.
  
  outname: name of output file
  
  ann_avg: Annual Averages (1/0)
  
  offset: Month Offset (0 for none)
  
  rotate: Rotate Longitudes (1/0)
  
  convert: Convert precipitation units (0-4)

See the full tutorial for their use at https://worldbuildingpasta.blogspot.com/2021/11/an-apple-pie-from-scratch-part-vi.html
