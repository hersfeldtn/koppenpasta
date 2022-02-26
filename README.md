# koppenpasta
A collection of scripts for use with ExoPlaSim, with the intention of allowing the user to import arbitrary topography and then produce climate maps of the output

Includes:

koppenpasta.py
  Script to take data from NetCDF files produced by ExoPlaSim and produce maps of Koppen climate zones (or Holdridge Life Zones).
  Requres NetCDF4, PIL, numpy, and scipy

defaultcolor.ini
  config file to allow user to select their own rgb values for the maps output by the above script

Image2sra_2.0.3.py
  Script to take greyscale heightmap images and convert them to scaled .sra files that can be read by ExoPlaSim as topography input. Created by Alex (Ostimeus), distributed with permission.
  Ost has made their own repository for this script bundled into an exoplasim configuration tool https://github.com/MegalexMaster/ExoPlaSim-InCon
  I'll leave this script in place here but it may not include later bugfixes (and the tool is nifty anyway).

See the full tutorial for their use at worldbuildingpasta.blogspot.com
