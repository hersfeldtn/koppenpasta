import os
from nco import Nco
from nco.custom import Atted
from netCDF4 import Dataset as ds

#List of data types for conversion:
precip = ['evap',
         'mrro',
         'pr',
         'prc',
         'prl',
         'prsn',
         'sndc',
         'snm',
         'wfn']
veg = ['veggpp',
       'veghresp',
       'veglitter',
       'vegnogrow',
       'vegnpp',
       'vegppl',
       'vegppw']

#Main Function
def main(in_files="", outname='out', ann_avg=0, offset=0, rotate=0, convert=0):
    nco = Nco()
    if isinstance(in_files, str):
        in_files = [in_files]
    for d in in_files:
        if os.path.isdir(d):
            for f in os.listdir(d):
                if f.endswith(".nc"):
                    in_files.append(d+"/"+f)
            in_files.remove(d)
    print(" Averaging input files...")
    nco.nces(input=in_files, output = outname+".nc")
    if convert > 0:
        print(" Converting units...")
        a = 60 * 60 * 24
        a_unit = 'mm day-1'
        b_unit = 'kg C m-2 day-1'
        if convert > 1:
            a *= 30
            a_unit = 'mm month-1'
            b_unit = 'kg C m-2 month-1'
        if convert > 2:
            a *= 12.175
            a_unit = 'mm year-1'
            b_unit = 'kg C m-2 year-1'
        b = a
        if convert < 4:
            a *= 1000
        else:
            a_unit = 'm year-1'
        opts = ["-O"]
        opts += ["-s '{0}={0}*{1}'".format(p,str(a)) for p in precip]
        opts += ["-s '{0}={0}*{1}'".format(v,str(b)) for v in veg]
        nco.ncap2(input=outname+".nc", output = outname+".nc", options=opts)
        com = "/usr/bin/ncatted -O" + "".join([" -a units,'{0}',o,c,'{1}'".format(p,a_unit) for p in precip]) + "".join([" -a units,'{0}',o,c,'{1}'".format(v,b_unit) for v in veg]) + " --output={0} {0}".format(outname+".nc")
        os.system(com)
        #for p in precip:
            #nco.ncap2(input=outname+".nc", output = outname+".nc", options=["-O -s '{0}={0}*{1}'".format(p,str(a))])
            #nco.ncatted(input=outname+".nc", output = outname+".nc", options=["-O units,\"{0}\",\"{1}\"".format(p,a_unit)])
            #os.system("/usr/bin/ncatted -O -a units,'{0}',o,c,'{1}', --output={2} {2}".format(p,a_unit,outname+".nc"))
        #opt = ["-O"]
        #opt += [Atted('o', 'units', p, a_unit) for p in precip]
        #nco.ncatted(input=outname+".nc", output = outname+".nc", options=opt)
        #nco.ncatted(input=outname+".nc", output = outname+".nc", options=['units',p,a_unit])
        #for v in veg:
            #nco.ncap2(input=outname+".nc", output = outname+".nc", options=["-O -s '{0}={0}*{1}'".format(v,str(b))])
            #nco.ncatted(input=outname+".nc", output = outname+".nc", options=["-O units,'{0}','{1}''".format(v,b_unit)])
    if offset != 0:
        print(" Offsetting months...")
        length = len(ds(in_files[0])['time'][:])
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
    while infile != path+"stop":
        infile = path+input('Input NetCDF filename or folder of files ("stop" for no more inputs): ')
        if infile in path+"stop":
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
    if rotate in ('y') or rotate in ('1'):
        rotate = 1
    convert = int(input('''
Convert Units
Converts precipitation units from m/s and vegetation production units from kg C / m^2 s to:
0: Do not convert
1: Convert to mm/day and kg C / m^2 day
2: Convert to mm/month and kg C / m^2 month
3: Convert to mm/year and kg C / m^2 year
4: Convert to m/year and kg C / m^2 year
Note: may take a while
Conversion: '''))
    ann_avg = input('''
Annual Averages
Outputs an additional .nc file which has every data type averaged to a single value representing the average across the whole year
Ouput annual averages (y/n): ''')
    if ann_avg in ('y') or ann_avg in ('1'):
        ann_avg = 1
    outname = path+input('''
Output Name: ''')

    main(in_files, outname, ann_avg, offset, rotate, convert)
