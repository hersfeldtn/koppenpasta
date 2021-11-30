import sys
import os
import numpy as np
import time
import math
from matplotlib.image import imread
from PIL import Image, ImageOps

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

def chunk(arr, nrows, ncols): #Split an array into chunks
    '''
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array should look like n subblocks with
    each subblock preserving the "physical" layout of arr.
    '''
    h, w = arr.shape
    assert h % nrows == 0, "{} rows is not evenly divisble by {}".format(h, nrows)
    assert w % ncols == 0, "{} cols is not evenly divisble by {}".format(w, ncols)
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))

def color_ocean(inarray, outlist, w, h): #turn 0/1 array to black/white (ocean/land) array
	for y in range(w):
		for x in range(h):
			col = inarray[x,y]
			if col > 0:
				rgb = (int(col*255),int(col*255),int(col*255))
			else:
				rgb = (0,0,0)
			outlist.append(rgb)
def color_land(inarray, outlist, w, h): #turn array to white (land) array
	for y in range(w):
		for x in range(h):
			rgb = (255,255,255)
			outlist.append(rgb)

#opening info
print('''
Image to .sra files
Script written 2021 by OstimeusAlex
For use with .sra input files for ExoPlaSim GCM
''')

#set path to current location
path = os.path.join(os.path.dirname(__file__), '')

#ask for input file
while True:
    infile = path+input('''
Filename

Input image filename: ''')
    if os.path.exists(infile):
        break
    print('No file found at '+str(infile))
print('File found')

#Prompts for configuration
#The yes/no prompts are intentionally flexible; any input containing 'y' will be read as 'yes', anything else will be read as 'no'
res = input('''
Planetary Oceans

Does your planet contain oceans? (y/n): ''')
if res in ('y'):
	oceans = 1 #Yes, this planet has oceans
	floor_value = int(input('''
Ocean Threshold

Black: 0, White: 255
Example: Entering '128' will mark values 128 and below as ocean.

Enter beginning ocean value: ''')) #at what value should land start?
	peak_value = float(input("Enter highest elevation (m above sea level): ")) #how heigh do the mountains get?
else:
	oceans = 0 #No, this planet is dry
	trench_value = abs(float(input("Enter deepest elevation (m): "))) #how deep does the land go?
	peak_value = float(input("Enter highest elevation (m): ")) #how heigh do the mountains get?
	range_value = trench_value+peak_value #what is the range between the depths and peaks?

#used to get geopotential
gravity = float(input('''
Planetary Gravity

Enter gravity of the planet (m/s^2)

Gravity: '''))

#self-explanatory
image_op = input('''
Image debug

Will produce images showing land masks of both the original image and the resized image, which may allow you to see if anything problematic has occured.

Would you like this script to produce debug images? (y/n): ''')
if image_op in ('y'):
	debug_img = 1
else:
	debug_img = 0

red_img = imread(infile)[:,:,0] #opens red channel
green_img = imread(infile)[:,:,1] #opens green channel
blue_img = imread(infile)[:,:,2] #opens blue channel
grey_img = (red_img+green_img+blue_img)/3 #averages color channels and scales from 0-1
img_width = len(grey_img[0]) #finds image width
img_height = len(grey_img) #finds image height
img_rgb = []

#converts 0-1 to elevation
if oceans == 1:
	greyscale_img = grey_img*255 #converts 0-1 to 0-255
	max_img = 255-(floor_value) #max value lowered by floor value
	rescaled_img = (greyscale_img-(floor_value)) #array lowered by floor value
	rescaled_img[rescaled_img <= 0] = 0 #any negative value becomes 0 (ocean)
	rescaled_img = (rescaled_img/max_img)*peak_value*gravity #array converted to 0-1, before multiplied by max height and gravity (geopotential)
	min = float(np.min(rescaled_img)) #finds the min height (should be 0)
	max = float(np.max(rescaled_img)) #finds the max height (should be near absolute max height)
	color_ocean(rescaled_img, img_rgb, img_width, img_height)
if oceans == 0:
	rescaled_img = ((grey_img*range_value)-trench_value)*gravity #converts 0-1 to trench-peak times gravity (geopotential)
	min = float(np.min(rescaled_img)) #finds the min height (should be near deepest depth)
	max = float(np.max(rescaled_img)) #finds the max height (should be near absolute max height)
	color_land(rescaled_img, img_rgb, img_width, img_height)

if debug_img == 1:
	img = Image.new('RGB',(img_height,img_width)) #print original image
	img.putdata(img_rgb)
	img = img.transpose(Image.ROTATE_90)
	img = ImageOps.flip(img)
	img.save(path+"OceanMaskOriginal.png")

#desired resolution
height = int(input('''
Latitude Resolution

Input latitude resolution
Example: '64' will produce a 128*64 output

Resolution: '''))
width = int(height * 2) #such that the result is always equirectangular

#print("Image width: "+str(img_width))
#print("Image height: "+str(img_height))
#print("Final width: "+str(width))
#print("Final height: "+str(height))
scale_width = img_width/width #finds the ratio between starting and final width
scale_height = img_height/height #finds the ratio between starting and final height
print("Width ratio: "+str(scale_width))
print("Height ratio: "+str(scale_height))

#tests to see whether original image is a multiple of desired resolution, before splitting the array into chunks if it is
if scale_width == math.floor(scale_width):
	if scale_height == math.floor(scale_height):
		width_ratio = int(scale_width) #width ratio between original and resized image
		height_ratio = int(scale_height) #height ratio between original and resized image
		chunked = chunk(rescaled_img, height_ratio, width_ratio) #splits the geopotential array into a 3d array
		length_ratio = len(chunked) #finds the length of the 3d array
		print("All ratios are good!")
	else:
		print("Height ratio incompatible, please try to have the image resolution as a multiple of the output resolution.")
		time.sleep(5)
		exit()
else:
	print("Width ratio incompatible, please try to have the image resolution as a multiple of the output resolution.")
	time.sleep(5)
	exit()

img_rgb = []
#different methods of rescaling depending on whether or not oceans are involved
if oceans == 1:
	b_w = chunked.copy()[:,:,:] #new empty 3d array, avoids overriding original array!
	for x in range(length_ratio):
		for y in range(width_ratio):
			for z in range(height_ratio):
				checker = chunked[x,y,z]
				if checker > 0: #converts 3d geopotential array to 3d land mask array
					b_w[x,y,z] = 1
				else:
					b_w[x,y,z] = 0
	b_w = b_w.mean(axis=(1,2)) #averages land mask array into 1d list
	dechunked = chunked.mean(axis=(1,2)) #simpler approach, since there's no ocean there's no need to check for it, so just average all the land to 1d list
	dechunked = np.reshape(dechunked, (height, width)) #reshapes list to equirectangular 2d array
	b_w = np.reshape(b_w, (height, width))
	color_ocean(b_w, img_rgb, width, height) #converts to land mask
if oceans == 0:
	dechunked = chunked.mean(axis=(1,2)) #simpler approach, since there's no ocean there's no need to check for it, so just average all the land to 1d list
	dechunked = np.reshape(dechunked, (height, width)) #reshapes list to equirectangular 2d array
	color_land(dechunked, img_rgb, width, height) #converts to land mask

if debug_img == 1:
	img = Image.new('RGB',(height,width)) #print resized image
	img.putdata(img_rgb)
	img = img.transpose(Image.ROTATE_90)
	img = ImageOps.flip(img)
	img.save(path+"OceanMaskSmall.png")

list_image = dechunked.flatten() #Flattens 2d array into list
b_w = b_w.flatten()
sra_129 = np.array(list_image).reshape(-1, 8) #Rearranges list into 2d array that matches sra format
sra_172 = np.array(b_w).reshape(-1, 8)
if oceans == 0:
	sra_172[:,:] = 1 #Everything becomes land
fl_height = float(height) #Apparently it doesn't like integers
fl_width = float(width)
print("Conversion successful!") #Nice!

name = input("Enter name of save file: ") #Might as well give it a custom name!
writeSRA(path,name,172,sra_172,fl_height,fl_width)
writeSRA(path,name,129,sra_129,fl_height,fl_width)
print("Saving successful!") #Nice!