import os
from nco import Nco
from netCDF4 import Dataset as ds

#Main Function
def main(in_files, outname, ann_avg, offset, rotate):
    nco = Nco()
    print(" Averaging input files...")
    nco.nces(input=in_files, output = outname+".nc")
    if offset != 0:
        print(" Offsetting months...")
        length = len(ds(in_files[0])['time'][:])
        print(offset+length)
        print(offset+2*length-1)
        nco.ncrcat(input=[outname+".nc",outname+".nc",outname+".nc"], output = outname+".nc", options=["-d time,"+str(offset+length)+","+str(offset+2*length-1)])
        
    if rotate == 1:
        print(" Rotating longitudes...")
        nco.ncap2(input=outname+".nc", output = outname+".nc", options=["-O -s 'lon=lon-180'"])

    if ann_avg == 1:
        print(" Producing annual averages...")
        nco.ncra(input=outname+".nc", output = outname+"_ann_av.nc")
    
    print("Done!")
    return

if __name__ == "__main__":
    print('''
ExoPlaSim Output Averager
Averages together data from multiple EPS .nc outputs
Made 2022 by Nikolai Hersfeldt

''')

    #set path to current location
    path = os.path.join(os.path.dirname(__file__), '')
    in_files = []

    #Configuration prompts
    infile = ""
    while infile != "stop":
        infile = path+input('Input NetCDF filename or folder of files ("stop" for no more inputs): ')
        if infile in "stop":
            print("")
        elif os.path.exists(infile):
            if os.path.isdir(infile):
                found = False
                for f in os.listdir(infile):
                    if f.endswith(".nc"):
                        in_files.append(infile+"/"+f)
                        found = True
                        print("  Found "+str(f))
                if found:
                    print(" Found all files")
                    break
                else:
                    print(" No files found in "+str(infile))
            else:
                in_files.append(infile)
                print(' File found')
        else:
            print(' No file found at '+str(infile))


    offset = int(input('''
Month Offset
Offsets output forward an integer number of months compared to input files
e.g. if a file starts in November, an offset of 2 will shift the start to January
Offset ("0" for no offset): '''))
    rotate = input('''
Rotate Longitudes
Convert longitudes from [0,360) range to [-180,180) range
with 0 longitude where 180 longitude formerly was
Rotate (y/n): ''')
    if rotate in ('y') or ann_avg in ('1'):
        rotate = 1
    ann_avg = input('''
Annual Averages
Outputs an additional .nc file which has every data type averaged to a single value representing the average across the whole year
Ouput annual averages (y/n): ''')
    if ann_avg in ('y') or ann_avg in ('1'):
        ann_avg = 1
    outname = path+input('''
Output Name: ''')

    main(in_files, outname, ann_avg, offset, rotate)

    
    
