import os
import numpy as np
from PIL import Image

seas_zero = True    #whether to ensure all areas with <0.5 land are kept flat at sea level

path = os.path.join(os.path.dirname(__file__), '')

def Prompt_bool(prompt):
    ans = input(prompt)
    if ('y') in ans or ('Y') in ans or ('1') in ans:
        return True
    else:
        return False

def Pick_from_list(prompt, opt_list):
    ans = input(prompt)
    while True:
        try:
            choice = opt_list[int(ans)]
            break
        except:
            ans = input(' Invalid input, try again: ')
    return choice
	
def File_search(fpath):
	if not os.path.exists(fpath):
		fpath = path+fpath
	if os.path.exists(fpath):
		try:
			im = Image.open(fpath)
			print(" File found")
			return im
		except:
			print(" Error: Could not open "+str(fpath)+" as image")
			return None
	else:
		print(" No file found at "+str(fpath))
		return None



def get_Input():
	print('''
im2sra_n file converter
Converts greyscale heightmaps to .sra topo input files for ExoPlaSim GCM
Written 2025 by Nikolai Hersfeldt
Based on Image2sra_2.0.3 written 2021 by OstimeusAlex
''')
	while True:
		map_im = File_search(str(input('Input heightmap file: ')))
		if map_im is not None:
			break
	
	fill_seas = False
	sealev = 0
	mask_im = None

	ocean_type = Pick_from_list('''
Ocean type
0: Flat oceans; global oceans at constant sea level, read from heightmap
1: Separate mask; upload separate land/sea mask image
2: No oceans; mark all areas as land
Select ocean type: ''',
		('flat','mask','none'))
	print('''
Elevation
Elevation extremes should apply to highest (whitest) and lowest (blackest) pixel in input image
All elevation in meters
''')
	maxel = float(input('Maximum elevation (m): '))
	minel = float(input('Minimum elevation (m): '))

	if ocean_type == 'flat':
		sealev = float(input('Sea level (m): '))
		if minel < sealev:
			fill_seas = Prompt_bool('''
Fill in all areas below sea level to sea level?
ExoPlaSim topography should be measured from ocean surface;
	including bathymetry (depth of ocean floor below water) will give incorrect results
If you have land areas below sea level, reccomend using separate mask option
Fill in oceans? (y/n): ''')
			
	elif ocean_type == 'mask':
		print('''
Land/Sea mask image
Should be filled in white for all land areas
and black for all sea areas''')
		while True:
			mask_im = File_search(str(input('Input mask file: ')))
			if mask_im is not None:
				break
	
	else:
		sealev = minel-100	#set below minimum elevation for no seas

	gravity = float(input('''
Surface gravity
Used to scale elevation to geopotential
Should be acceleration in meters/seconds^2
Earth gravity is 9.81
Gravity (m/s^2): '''))
	
	res = int(input('''
Latitude resolution
Height of the grid in the ExoPlaSim model
Grid longitudinal width will be twice this value
Appropriate values for standard resolutions:
	T21:	  32
	T42:	  64
	T63:	  96
	T85:	 128
	T106:	 160
	T127:	 192
	T170:	 256
Resolution: '''))
	
	debug = Prompt_bool('''
Debug images
Will produce following images:
	debug_mask_big.png      Land/Sea mask at input resolution
	debug_mask_small.png    Land/Sea mask at model resolution
	debug_topo_small.png    Topography heightmap at model resolution
Produce debug images? (y/n): ''')
	
	outname = str(input('''
Output name: '''))
	
	return map_im, mask_im, maxel, minel, sealev, gravity, res, fill_seas, debug, outname

#Reads topography from greyscale image and scales to geopotential
# returns elevation array and mask
def Read_topo(map_im, mask_im, maxel, minel, sealev, gravity, res, fill_seas):

	print('Reading topography file...')
	hmap = map_im.convert('F')    # Convert to float for precision on downscaling
	hmapext = hmap.getextrema() # Find minimum and maximum before resizing for accurate elevation scaling

	if mask_im is None:
		if minel < sealev:
			hmap_bin_ar = np.ones_like(np.asarray(hmap)) * 255	#all land
		else:
			print(' Extracting land/sea mask...')
			thresh = (sealev - minel) * (hmapext[1] - hmapext[0]) / (maxel - minel) # Find greyscale value corresponding to sea level
			hmap_bin_ar = np.where(np.asarray(hmap) > thresh, 255, 0) #convert to array, binarize, then convert back to image
		try:
			hmap_bin = Image.fromarray(hmap_bin_ar)
		except:
			hmap_bin = Image.fromarray(hmap_bin_ar.astype(np.uint8))    #I dunno why the first version fails for some people, this might help
	else:
		print('Reading land/sea mask file...')
		hmap_bin = mask_im.convert('F')
	print(' Downscaling maps...')
	hmap_bin_sm = hmap_bin.resize((res*2, res), Image.Resampling.BILINEAR)
	topo_bin = np.asarray(hmap_bin_sm)
	topo_mask = topo_bin/255
	hmap = hmap.resize((res*2, res), Image.Resampling.BILINEAR)
	elev = np.asarray(hmap)
	elev = elev * (maxel - minel) / (hmapext[1] - hmapext[0]) * gravity #scale to geopotential
	if fill_seas:
		print(' Filling seas to sea level...')
		elev = np.where(elev < sealev*gravity, sealev*gravity, elev)
	if seas_zero and (fill_seas or (mask_im is None and sealev == minel)):  #fill in all sea cells as exactly sea level
		elev = np.where(topo_mask < 0.5, sealev, elev)
	return elev, topo_mask, hmap_bin, hmap_bin_sm, hmap

#Ost's code from Image2sra
def writeSRA(path,name,kcode,fmap,NLAT,NLON): #Format array and header into sra file, as well as saving it!
    """Write a lat-lon field to a formatted .sra file
    """
    label=os.path.join(path,name+'_surf_%04d.sra'%kcode)
    header=[kcode,0,11111111,0,NLON,NLAT,0,0]
    sheader = ''
    for h in header:
        sheader+=" %9d"%h
    lines=[]
    i=0
    while i<NLAT*NLON/8:
        l=''
        for n in fmap[i,:]:
            l+=' %9.3f'%n
        lines.append(l)
        i+=1
    text=sheader+'\n'+'\n'.join(lines)
    f=open(label,'w')
    f.write(text)
    f.close()


def Make_files(map_im, mask_im=None, maxel=1000, minel=0, sealev=0, gravity=9.81, res=64, fill_seas=False, debug=False, outname='out'):
	if isinstance(map_im, str):
		try:
			map_im = Image.open(map_im)
		except:
			map_im = Image.open(path+map_im)
                
	if isinstance(mask_im, str):
		try:
			map_im = Image.open(mask_im)
		except:
			map_im = Image.open(path+mask_im)
	        
	elev, mask, mask_big, mask_small, height_small = Read_topo(map_im, mask_im, maxel, minel, sealev, gravity, res, fill_seas)

	if debug:
		print('Saving debug images...')
		for im, n in zip((mask_big, mask_small, height_small), ('mask_big', 'mask_small', 'topo_small')):
			outim = im.convert('L')
			outim.save(path+'debug_'+n+'.png')
	
	print('Converting to sra files...')

	#Ost's code from image2sra
	list_image = elev.flatten() #Flattens 2d array into list
	b_w = mask.flatten()
	sra_129 = np.array(list_image).reshape(-1, 8) #Rearranges list into 2d array that matches sra format
	sra_172 = np.array(b_w).reshape(-1, 8)
	fl_height = float(res) #Apparently it doesn't like integers
	fl_width = float(res*2)
	print("Conversion successful!") #Nice!
	writeSRA(path,outname,172,sra_172,fl_height,fl_width)
	print('Saved land/sea mask data to '+path+outname+'_surf_172.sra')
	writeSRA(path,outname,129,sra_129,fl_height,fl_width)
	print('Saved topography data to '+path+outname+'_surf_129.sra')
	print('Conversion complete!')


if __name__ == "__main__":
	map_im, mask_im, maxel, minel, sealev, gravity, res, fill_seas, debug, outname = get_Input()
	Make_files(map_im, mask_im, maxel, minel, sealev, gravity, res, fill_seas, debug, outname)
	a = input('Press enter to close ')

