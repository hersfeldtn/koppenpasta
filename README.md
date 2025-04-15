# koppenpasta
A collection of scripts for use with ExoPlaSim, with the intention of allowing the user to import arbitrary topography and then produce climate maps of the output
Currently supported climate classification systems are:
- Koppen-Geiger
- Koppen-Trewartha
- Holdridge Life Zones
- Thornthwaite-Feddemma
- Prentice et al. 1992 BIOME1 model
- Hersfeldt/Pasta Bioclimate Zones
- Whittaker Biomes
- Woodward Vegetation Types
- IPCC Climate Zones
- World Climate Regions
- Two-Parameter Koppen-Alike
- Koppen-Geiger Unproxied

Includes:

koppenpasta.py
  Script to take data from NetCDF files produced by ExoPlaSim and produce maps of climate zones (or Holdridge Life Zones).
  Requres NetCDF4, PIL, numpy, and scipy
  Can be run directly and used for command-line setup
  Or can be loaded and run in python with koppenpasta.Make_map(files,opts)
    where files is the name of a file, a list of files, or a directory which will be searched for .nc files
    and opts is a dictionary of options or name of a config file, or can be left out to run with defaults

kpasta_options.cfg
  Example options, consult the internal documentation for use of each option.
  If in same directory as koppenpasta, the script will automatically load options to override its internal defaults
  Can be renamed and altered, then selected for use during script setup or by Make_Map

defaultcolor.ini
  config file to allow user to select their own rgb values for the maps output by the above script

im2sra_n.py
  Script to take greyscale heightmap images and convert them to scaled .sra files that can be read by ExoPlaSim as topography input.
  Based on Image2sra below (and retains some of Alex's code) but updated for more flexible input options.

Image2sra_2.0.3.py
  Original script for converting heightmaps to .sra files. Created by Alex (Ostimeus), distributed with permission.
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
