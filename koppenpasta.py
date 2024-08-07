import netCDF4 as nc
import statistics as stat
import math
import numpy as np
from scipy.interpolate import RectSphereBivariateSpline as spint
import configparser
from PIL import Image
import os

#Assign integers to climate zones, because it's easier than working with arrays of strings 
#Note: by default these are stored in the climate zone arrays as 8-bit signed integers, so they must remain within the range (-127,127)
# though it wouldn't be too hard to just alter the data type for that array if necessary
Af=1
Am=2
As=3
Aw=4
BWh=5
BWk=6
BSh=7
BSk=8
Csa=9
Csb=10
Csc=11
Cwa=12
Cwb=13
Cwc=14
Cfa=15
Cfb=16
Cfc=17
Dsa=18
Dsb=19
Dsc=20
Dsd=21
Dwa=22
Dwb=23
Dwc=24
Dwd=25
Dfa=26
Dfb=27
Dfc=28
Dfd=29
ET=30
EF=31

A=32
B=33
C=34
D=35
E=36

TropRainforest=37
TropMonsoon=38
TropSavanna=39
HotDesert=40
HotSteppe=41
ColdDesert=42
ColdSteppe=43
Med=44
Subtropical=45
Oceanic=46
Continental=47
Subarctic=48
Tundra=49
IceCap=50

TropRainforestnos=51
TropSavannanos=52
HotDesertnos=53
HotSteppenos=54
ColdDesertnos=55
ColdSteppenos=56
Oceanicnos=57
Tundranos=58
IceCapnos=59

SeaTrop=60
SeaTemp=61
SeaSeasonalIce=62
SeaPermIce=63
SeaFlat=64

SeaTropnos=65
SeaTempnos=66
SeaPermIcenos=67

HPolDesert=68
HDryTundra=69
HMoistTundra=70
HWetTundra=71
HRainTundra=72
HBorDesert=73
HDryScrub=74
HBorMoistForest=75
HBorWetForest=76
HBorRainForest=77
HCoolDesert=78
HCoolDesertScrub=79
HSteppe=80
HCoolMoistForest=81
HCoolWetForest=82
HCoolRainForest=83
HWarmDesert=84
HWarmDesertScrub=85
HThornSteppe=86
HWarmDryForest=87
HWarmMoistForest=88
HWarmWetForest=89
HWarmRainForest=90
HTropDesert=91
HTropDesertScrub=92
HThornWood=93
HVDryForest=94
HTropDryForest=95
HTropMoistForest=96
HTropWetForest=97
HTropRainForest=98

Sea_Zones=1

#Create list of colors corresponding to each climate type, depending on the color list
def make_colmap(color_type, col_list_path):
    colmap = np.empty((255,3),dtype=np.uint8)
    colmap[0] = [0,0,0] #if no recognized climate zone is attached to this position, mark it black for debug purposes
    if color_type == 3:
        #Reads custom color list file
        cfg = configparser.ConfigParser()
        cfg.read(col_list_path)
        fulkop = cfg['Full Koppen']
        kopgroup = cfg['Koppen Groups']
        redkop = cfg['Reduced Koppen']
        noskop = cfg['Seasonless Koppen']
        seazon = cfg['Sea Zones']
        redsea = cfg['Reduced Sea']
        hold = cfg['Holdridge']
        nozon = cfg['Debug']
        colist = np.empty((255),dtype=str)
        colist[0] = nozon.get('None')
        colist[Af] = fulkop.get('Af')
        colist[Am] = fulkop.get('Am')
        colist[Aw] = fulkop.get('Aw')
        colist[As] = fulkop.get('As')
        colist[BWh] = fulkop.get('BWh')
        colist[BWk] = fulkop.get('BWk')
        colist[BSh] = fulkop.get('BSh')
        colist[BSk] = fulkop.get('BSk')
        colist[Csa] = fulkop.get('Csa')
        colist[Csb] = fulkop.get('Csb')
        colist[Csc] = fulkop.get('Csc')
        colist[Cwa] = fulkop.get('Cwa')
        colist[Cwb] = fulkop.get('Cwb')
        colist[Cwc] = fulkop.get('Cwc')
        colist[Cfa] = fulkop.get('Cfa')
        colist[Cfb] = fulkop.get('Cfb')
        colist[Cfc] = fulkop.get('Cfc')
        colist[Dsa] = fulkop.get('Dsa')
        colist[Dsb] = fulkop.get('Dsb')
        colist[Dsc] = fulkop.get('Dsc')
        colist[Dsd] = fulkop.get('Dsd')
        colist[Dwa] = fulkop.get('Dwa')
        colist[Dwb] = fulkop.get('Dwb')
        colist[Dwc] = fulkop.get('Dwc')
        colist[Dwd] = fulkop.get('Dwd')
        colist[Dfa] = fulkop.get('Dfa')
        colist[Dfb] = fulkop.get('Dfb')
        colist[Dfc] = fulkop.get('Dfc')
        colist[Dfd] = fulkop.get('Dfd')
        colist[ET] = fulkop.get('ET')
        colist[EF] = fulkop.get('EF')
        colist[A] = kopgroup.get('A')
        colist[B] = kopgroup.get('B')
        colist[C] = kopgroup.get('C')
        colist[D] = kopgroup.get('D')
        colist[E] = kopgroup.get('E')
        colist[TropRainforest] = redkop.get('TropRainforest')
        colist[TropMonsoon] = redkop.get('TropMonsoon')
        colist[TropSavanna] = redkop.get('TropSavanna')
        colist[HotDesert] = redkop.get('HotDesert')
        colist[HotSteppe] = redkop.get('HotSteppe')
        colist[ColdDesert] = redkop.get('ColdDesert')
        colist[ColdSteppe] = redkop.get('ColdSteppe')
        colist[Med] = redkop.get('Med')
        colist[Subtropical] = redkop.get('Subtropical')
        colist[Oceanic] = redkop.get('Oceanic')
        colist[Continental] = redkop.get('Continental')
        colist[Subarctic] = redkop.get('Subarctic')
        colist[Tundra] = redkop.get('Tundra')
        colist[IceCap] = redkop.get('IceCap')
        colist[TropRainforestnos] = noskop.get('TropRainforestnos')
        colist[TropSavannanos] = noskop.get('TropSavannanos')
        colist[HotDesertnos] = noskop.get('HotDesertnos')
        colist[HotSteppenos] = noskop.get('HotSteppenos')
        colist[ColdDesertnos] = noskop.get('ColdDesertnos')
        colist[ColdSteppenos] = noskop.get('ColdSteppenos')
        colist[Oceanicnos] = noskop.get('Oceanicnos')
        colist[Tundranos] = noskop.get('Tundranos')
        colist[IceCapnos] = noskop.get('IceCapnos')
        colist[SeaTrop] = seazon.get('SeaTrop')
        colist[SeaTemp] = seazon.get('SeaTemp')
        colist[SeaSeasonalIce] = seazon.get('SeaSeasonalIce')
        colist[SeaPermIce] = seazon.get('SeaPermIce')
        colist[SeaTropnos] = redsea.get('SeaTropnos')
        colist[SeaTempnos] = redsea.get('SeaTempnos')
        colist[SeaPermIcenos] = redsea.get('SeaPermIcenos')
        colist[SeaFlat] = redsea.get('SeaFlat')
        colist[HPolDesert] = hold.get('HPolDesert')
        colist[HDryTundra] = hold.get('HDryTundra')
        colist[HMoistTundra] = hold.get('HMoistTundra')
        colist[HWetTundra] = hold.get('HWetTundra')
        colist[HRainTundra] = hold.get('HRainTundra')
        colist[HBorDesert] = hold.get('HBorDesert')
        colist[HDryScrub] = hold.get('HDryScrub')
        colist[HBorMoistForest] = hold.get('HBorMoistForest')
        colist[HBorWetForest] = hold.get('HBorWetForest')
        colist[HBorRainForest] = hold.get('HBorRainForest')
        colist[HCoolDesert] = hold.get('HCoolDesert')
        colist[HCoolDesertScrub] = hold.get('HCoolDesertScrub')
        colist[HSteppe] = hold.get('HSteppe')
        colist[HCoolMoistForest] = hold.get('HCoolMoistForest')
        colist[HCoolWetForest] = hold.get('HCoolWetForest')
        colist[HCoolRainForest] = hold.get('HCoolRainForest')
        colist[HWarmDesert] = hold.get('HWarmDesert')
        colist[HWarmDesertScrub] = hold.get('HWarmDesertScrub')
        colist[HThornSteppe] = hold.get('HThornSteppe')
        colist[HWarmDryForest] = hold.get('HWarmDryForest')
        colist[HWarmMoistForest] = hold.get('HWarmMoistForest')
        colist[HWarmWetForest] = hold.get('HWarmWetForest')
        colist[HWarmRainForest] = hold.get('HWarmRainForest')
        colist[HTropDesert] = hold.get('HTropDesert')
        colist[HTropDesertScrub] = hold.get('HTropDesertScrub')
        colist[HThornWood] = hold.get('HThornWood')
        colist[HVDryForest] = hold.get('HVDryForest')
        colist[HTropDryForest] = hold.get('HTropDryForest')
        colist[HTropMoistForest] = hold.get('HTropMoistForest')
        colist[HTropWetForest] = hold.get('HTropWetForest')
        colist[HTropRainForest] = hold.get('HTropRainForest')
        for i in range(len(colist)):
            co = colist[i]
            if ',' in co:
                co = co.split(',')
                co = [int(j) for j in co]
                co = tuple(co)
                colmap[i] = np.asarray[co]
    else:
        colmap[HPolDesert] = [255,255,255]
        colmap[HDryTundra] = [128,128,128]
        colmap[HMoistTundra] = [96,128,128]
        colmap[HWetTundra] = [64,128,144]
        colmap[HRainTundra] = [32,128,192]
        colmap[HBorDesert] = [160,160,128]
        colmap[HDryScrub] = [128,160,128]
        colmap[HBorMoistForest] = [96,160,128]
        colmap[HBorWetForest] = [64,160,144]
        colmap[HBorRainForest] = [32,160,192]
        colmap[HCoolDesert] = [192,192,128]
        colmap[HCoolDesertScrub] = [160,192,128]
        colmap[HSteppe] = [128,192,128]
        colmap[HCoolMoistForest] = [96,192,128]
        colmap[HCoolWetForest] = [64,192,144]
        colmap[HCoolRainForest] = [32,192,192]
        colmap[HWarmDesert] = [224,224,128]
        colmap[HWarmDesertScrub] = [192,224,128]
        colmap[HThornSteppe] = [160,224,128]
        colmap[HWarmDryForest] = [128,224,128]
        colmap[HWarmMoistForest] = [96,224,128]
        colmap[HWarmWetForest] = [64,224,144]
        colmap[HWarmRainForest] = [32,224,192]
        colmap[HTropDesert] = [255,255,128]
        colmap[HTropDesertScrub] = [224,255,128]
        colmap[HThornWood] = [192,255,128]
        colmap[HVDryForest] = [160,255,128]
        colmap[HTropDryForest] = [128,255,128]
        colmap[HTropMoistForest] = [96,255,128]
        colmap[HTropWetForest] = [64,255,144]
        colmap[HTropRainForest] = [32,255,160]
        colmap[SeaTrop] = [9,120,171]
        colmap[SeaTropnos] = [9,120,171]
        colmap[SeaTemp] = [113,171,216]
        colmap[SeaTempnos] = [113,171,216]
        colmap[SeaFlat] = [113,171,216]
        colmap[SeaSeasonalIce] = [185,227,255]
        colmap[SeaPermIce] = [226,248,255]
        colmap[SeaPermIcenos] = [226,248,255]
        if color_type == 0:
            colmap[Af] = [0,0,254]
            colmap[TropRainforest] = [0,0,254]
            colmap[TropRainforestnos] = [0,0,254]
            colmap[A] = [0,0,254]
            colmap[Am] = [0,119,255]
            colmap[TropMonsoon] = [0,119,255]
            colmap[Aw] = [70,169,250]
            colmap[TropSavanna] = [70,169,250]
            colmap[TropSavannanos] = [70,169,250]
            colmap[As] = [127,201,255]
            colmap[BWh] = [254,0,0]
            colmap[HotDesert] = [254,0,0]
            colmap[HotDesertnos] = [254,0,0]
            colmap[B] = [254,0,0]
            colmap[BWk] = [254,150,149]
            colmap[ColdDesert] = [254,150,149]
            colmap[ColdDesertnos] = [254,150,149]
            colmap[BSh] = [245,163,1]
            colmap[HotSteppe] = [245,163,1]
            colmap[HotSteppenos] = [245,163,1]
            colmap[BSk] = [255,219,99]
            colmap[ColdSteppe] = [255,219,99]
            colmap[ColdSteppenos] = [255,219,99]
            colmap[Csa] = [255,255,0]
            colmap[Med] = [255,255,0]
            colmap[Csb] = [198,199,0]
            colmap[Csc] = [150,150,0]
            colmap[Cwa] = [150,255,150]
            colmap[Subtropical] = [150,255,150]
            colmap[Cwb] = [99,199,100]
            colmap[Cwc] = [50,150,51]
            colmap[Cfa] = [198,255,78]
            colmap[Cfb] = [102,255,51]
            colmap[Cfc] = [51,199,1]
            colmap[Oceanic] = [51,199,1]
            colmap[Oceanicnos] = [51,199,1]
            colmap[C] = [51,199,1]
            colmap[Dsa] = [255,0,254]
            colmap[Dsb] = [198,0,199]
            colmap[Dsc] = [150,50,149]
            colmap[Dsd] = [150,100,149]
            colmap[Dwa] = [171,177,255]
            colmap[Dwb] = [90,119,219]
            colmap[Dwc] = [76,81,181]
            colmap[Dwd] = [50,0,135]
            colmap[Dfa] = [0,255,255]
            colmap[Continental] = [0,255,255]
            colmap[D] = [0,255,255]
            colmap[Dfb] = [56,199,255]
            colmap[Dfc] = [0,126,125]
            colmap[Subarctic] = [0,126,125]
            colmap[Dfd] = [0,69,94]
            colmap[ET] = [178,178,178]
            colmap[Tundra] = [178,178,178]
            colmap[Tundranos] = [178,178,178]
            colmap[EF] = [104,104,104]
            colmap[IceCap] = [104,104,104]
            colmap[IceCapnos] = [104,104,104]
            colmap[E] = [104,104,104] 
        elif color_type == 1:
            colmap[Af] = [147,1,1]
            colmap[TropRainforest] = [147,1,1]
            colmap[TropRainforestnos] = [147,1,1]
            colmap[A] = [147,1,1]
            colmap[Am] = [254,0,0]
            colmap[TropMonsoon] = [254,0,0]
            colmap[Aw] = [255,207,207]
            colmap[TropSavanna] = [255,207,207]
            colmap[TropSavannanos] = [255,207,207]
            colmap[As] = [254,154,154]
            colmap[BWh] = [255,207,0]
            colmap[HotDesert] = [255,207,0]
            colmap[HotDesertnos] = [255,207,0]
            colmap[B] = [255,207,0]
            colmap[BWk] = [255,254,101]
            colmap[ColdDesert] = [255,254,101]
            colmap[ColdDesertnos] = [255,254,101]
            colmap[BSh] = [207,143,20]
            colmap[HotSteppe] = [207,143,20]
            colmap[HotSteppenos] = [207,143,20]
            colmap[BSk] = [206,170,84]
            colmap[ColdSteppe] = [206,170,84]
            colmap[ColdSteppenos] = [206,170,84]
            colmap[Csa] = [0,254,0]
            colmap[Med] = [0,254,0]
            colmap[Csb] = [149,255,0]
            colmap[Csc] = [203,255,0]
            colmap[Cwa] = [180,101,0]
            colmap[Subtropical] = [180,101,0]
            colmap[Cwb] = [150,102,4]
            colmap[Cwc] = [94,64,0]
            colmap[Cfa] = [0,48,0]
            colmap[Cfb] = [1,80,1]
            colmap[Oceanic] = [1,80,1]
            colmap[Oceanicnos] = [1,80,1]
            colmap[C] = [1,80,1]
            colmap[Cfc] = [0,120,0]
            colmap[Dsa] = [254,108,253]
            colmap[Dsb] = [254,182,255]
            colmap[Dsc] = [230,202,253]
            colmap[Dsd] = [202,204,203]
            colmap[Dwa] = [204,182,255]
            colmap[Dwb] = [153,124,178]
            colmap[Dwc] = [138,89,178]
            colmap[Dwd] = [109,36,179]
            colmap[Dfa] = [48,0,48]
            colmap[Continental] = [48,0,48]
            colmap[D] = [48,0,48]
            colmap[Dfb] = [101,1,100]
            colmap[Dfc] = [203,0,203]
            colmap[Subarctic] = [203,0,203]
            colmap[Dfd] = [199,21,135]
            colmap[ET] = [101,255,255]
            colmap[Tundra] = [101,255,255]
            colmap[Tundranos] = [101,255,255]
            colmap[EF] = [99,150,255]
            colmap[IceCap] = [99,150,255]
            colmap[IceCapnos] = [99,150,255]
            colmap[E] = [99,150,255]
        elif color_type == 2:
            colmap[Af] = [42,65,15]
            colmap[TropRainforest] = [42,65,15]
            colmap[TropRainforestnos] = [42,65,15]
            colmap[A] = [42,65,15]
            colmap[Am] = [53,74,19]
            colmap[TropMonsoon] = [53,74,19]
            colmap[Aw] = [73,87,32]
            colmap[As] = [73,87,32]
            colmap[TropSavanna] = [73,87,32]
            colmap[TropSavannanos] = [73,87,32]
            colmap[BWh] = [213,183,133]
            colmap[HotDesert] = [213,183,133]
            colmap[HotDesertnos] = [213,183,133]
            colmap[B] = [213,183,133]
            colmap[BWk] = [178,153,112]
            colmap[ColdDesert] = [178,153,112]
            colmap[ColdDesertnos] = [178,153,112]
            colmap[BSh] = [123,112,66]
            colmap[HotSteppe] = [123,112,66]
            colmap[HotSteppenos] = [123,112,66]
            colmap[BSk] = [128,117,74]
            colmap[ColdSteppe] = [128,117,74]
            colmap[ColdSteppenos] = [128,117,74]
            colmap[Csa] = [112,104,58]
            colmap[Med] = [112,104,58]
            colmap[Csb] = [66,75,31]
            colmap[Csc] = [66,75,31]
            colmap[Cwa] = [78,88,36]
            colmap[Subtropical] = [78,88,36]
            colmap[Cwb] = [79,81,38]
            colmap[Cwc] = [135,114,68]
            colmap[Cfa] = [65,80,27]
            colmap[Cfb] = [62,77,27]
            colmap[Oceanic] = [62,77,27]
            colmap[Oceanicnos] = [62,77,27]
            colmap[C] = [62,77,27]
            colmap[Cfc] = [68,78,48]
            colmap[Dsa] = [148,131,85]
            colmap[Dsb] = [96,92,50]
            colmap[Dsc] = [66,70,31]
            colmap[Dsd] = [57,66,23]
            colmap[Dwa] = [74,89,34]
            colmap[Dwb] = [67,83,32]
            colmap[Dwc] = [57,72,23]
            colmap[Dwd] = [65,71,28]
            colmap[Dfa] = [66,85,29]
            colmap[Continental] = [66,85,29]
            colmap[D] = [66,85,29]
            colmap[Dfb] = [55,75,21]
            colmap[Dfc] = [52,64,20]
            colmap[Subarctic] = [52,64,20]
            colmap[Dfd] = [62,71,25]
            colmap[ET] = [122,119,92]
            colmap[Tundra] = [122,119,92]
            colmap[Tundranos] = [122,119,92]
            colmap[EF] = [204,217,232]
            colmap[IceCap] = [204,217,232]
            colmap[IceCapnos] = [204,217,232]
            colmap[E] = [204,217,232]
            colmap[SeaTrop] = [20,30,66]
            colmap[SeaTropnos] = [20,30,66]
            colmap[SeaTemp] = [20,30,66]
            colmap[SeaTempnos] = [20,30,66]
            colmap[SeaFlat] = [20,30,66]
            colmap[SeaSeasonalIce] = [20,30,66]
            colmap[SeaPermIce] = [190,208,226]
            colmap[SeaPermIcenos] = [190,208,226]
    return colmap

#defines a function to convert climate zone array to a list of rgb values
def color(inarray, colmap, latl, lonl, lsm, blend, land_type, sea_type, color_type, col_list_path, inarraysea=Sea_Zones):
    if blend == 1 or sea_type == 5: #if blending is off or there is no sea, just reads from land climate array
        climmap = inarray
    else:
        lsm = np.swapaxes(lsm,1,2)
        climmap = np.where(lsm[0,:] >= 0.5, inarray, inarraysea) #lsm should be binary, but just in case it isn't in the future I set the land/sea boundary at 0.5
    climmap = np.where(climmap > 0.1, climmap, 0) #Should fill any empty spaces with 0
    outmap = np.empty((lonl,latl,3),dtype=np.uint8)
    for y in range(latl):
        for x in range(lonl):
            clim = climmap[x,y]
            outmap[x,y] = colmap[clim]
    outmap = np.swapaxes(outmap,0,1)
    return outmap

def MakeMap(in_files="", output_name="output", cfg_loadname="",
         land_type=0,
         sea_type=0,
         blend=0,
         color_type=0,
         col_list_path="",
         bin_num=1,
         interp=0,
         dum_ice=0,
         use_topo=0,
         topo_path="",
         maxel=0,
         minel=0,
         gravity=0,
         blend_topo=0,
         sealev=0,
         sum_def=0,
         temp_def=0,
         arid_def=0,
         med_def=0,
         wet_def=0,
         add_def=0,
         prio_def=0,
         ice_def=0,
         sea_def=0
         ):

    if isinstance(in_files, str):
        in_files = [in_files]
    for d in in_files:
        if os.path.isdir(d):
            for f in os.listdir(d):
                if f.endswith(".nc"):
                    in_files.append(d+"/"+f)
            in_files.remove(d)
    in_num = len(in_files)
    if in_num > 1:
        print('\nInterpreting '+in_files[0]+' et al')
    else:
        print('\nInterpreting '+in_files[0])
    if cfg_loadname!="":
        print('Loading Config '+cfg_loadname)
        cfg = configparser.ConfigParser()
        cfg.read(cfg_loadname)
        typ = cfg['Zone Types']
        land_type = typ.getint('Land Climate Zones')
        sea_type = typ.getint('Sea Climate Zones')
        blend = typ.getint('Land/Sea Blend')
        color_type = typ.getint('Color List')
        col_list_path = typ.get('Custom Color List')
        adv = cfg['Advanced Functions']
        bin_num = adv.getint('Bin Size')
        interp = adv.getint('Interpolation Factor')
        dum_ice = adv.getint('Dummy Sea Ice')
        use_topo = adv.getint('Adjust Temp by Topography')
        topo_path = adv.get('Topography Map')
        maxel = adv.getfloat('Max Elevation')
        minel = adv.getfloat('Min Elevation')
        gravity = adv.getfloat('Surface Gravity')
        blend_topo = adv.getint('Blend by Topography Map')
        sealev = adv.getfloat('Sea Level')
        defs = cfg['Climate Zone Definitions']
        sum_def = defs.getint('Summer Definition')
        temp_def = defs.getint('Temperate/Continental Boundary')
        arid_def = defs.getint('Hot/Cold Arid Boundary')
        med_def = defs.getint('Mediterranean Definition')
        wet_def = defs.getint('Wet-Summer Definition')
        add_def = defs.getint('Additional Wet Season Requirements')
        prio_def = defs.getint('Med/Wet-Summer Priority')
        ice_def = defs.getint('Arid/Polar Priority')
        sea_def = defs.getint('Sea Ice Definition')

    
    if sea_type == 0 or sea_type == 2:   #checks if we need ts arrays for sea surface temperature
        use_ts = 1
    elif sea_def == 1 and (sea_type == 1 or sea_type == 3):
        use_ts = 1
    else:
        use_ts = 0

    

    #retrieve data from NetCDF file
    print('Retrieving Data...')
    ds = nc.Dataset(in_files[0])
    time = ds['time'][:]    #month
    lat = ds['lat'][:]  #latitude
    lon = ds['lon'][:]  #longitude
    timefull = int(len(time))   #number of months in input file(s)
    timel = math.ceil(timefull/bin_num)   #number of months after binning
    if timel%2 == 0:
        halfl = int(timel/2)
        timodd = 0
    else:   #if number of months is odd, round down for number of months used for definitions depending on half-year averages
        halfl = int((timel-1)/2)
        timodd = 1
    thirdl = int(round(timel/3))
    lonl = len(lon)
    latl = len(lat)
    

    if in_num < 2:
        tas_ar = ds['tas'][:]    #2-meter air temperature, Kelvin
        pr_ar = ds['pr'][:]    #precipitation, meters/second
        lsm_ar = ds['lsm'][:]  #binary land/sea mask
        sic_ar = ds['sic'][:]  #fractional sea ice cover
        if use_topo == 1:
            grnz_ar = ds['grnz'][:]   #surface geopotential, meters^2/seconds^2
        if sum_def == 1:
            czen_ar = ds['czen'][:]   #cosine of solar zenith angle
        if use_ts == 1:
            ts_ar = ds['ts'][:]    #surface temperature, Kelvin
    else:   #routine to average together data across multiple input files
        tas_raw = np.empty((in_num,timefull,latl,lonl))  #create 4d arrays containing data from all files
        pr_raw = np.empty((in_num,timefull,latl,lonl))
        lsm_raw = np.empty((in_num,timefull,latl,lonl))
        sic_raw = np.empty((in_num,timefull,latl,lonl))
        if use_topo == 1:
            grnz_raw = np.empty((in_num,timefull,latl,lonl))
        if sum_def == 1:
            czen_raw = np.empty((in_num,timefull,latl,lonl))
        if use_ts == 1:
            ts_raw = np.empty((in_num,timefull,latl,lonl))
        for n in range(in_num):
            ds = nc.Dataset(in_files[n])
            tas_raw[n,:,:,:] = ds['tas'][:]   #load data into those arrays, file by file
            pr_raw[n,:,:,:] = ds['pr'][:]
            lsm_raw[n,:,:,:] = ds['lsm'][:]
            sic_raw[n,:,:,:] = ds['sic'][:]
            if use_topo == 1:
                grnz_raw[n,:,:,:] = ds['grnz'][:]
            if sum_def == 1:
                czen_raw[n,:,:,:] = ds['czen'][:]
            if use_ts == 1:
                ts_raw[n,:,:,:] = ds['ts'][:]
        tas_ar = np.mean(tas_raw, axis=0)     #average the data from the 4d arrays into 3d arrays
        pr_ar = np.mean(pr_raw, axis=0)
        lsm_ar = np.mean(lsm_raw, axis=0)
        sic_ar = np.mean(sic_raw, axis=0)
        if use_topo == 1:
            grnz_ar = np.mean(grnz_raw, axis=0)
        if sum_def == 1:
            czen_ar = np.mean(czen_raw, axis=0)
        if use_ts == 1:
            ts_ar = np.mean(ts_raw, axis=0)
    
    
    if use_topo == 1:
        grnz_bin = np.mean(grnz_ar, axis=0)    #geopotential is averaged into a 2d array, as it doesn't vary with time

    #Routine to bin adjacent months together
    if bin_num > 1:
        print('Binning Months...')
        tas_bin = np.empty((timel,latl,lonl))    #create new arrays with fewer months
        pr_bin = np.empty((timel,latl,lonl))
        lsm_bin = np.empty((timel,latl,lonl))
        sic_bin = np.empty((timel,latl,lonl))
        if sum_def == 1:
            czen_bin = np.empty((timel,latl,lonl))
        if use_ts == 1:
            ts_bin = np.empty((timel,latl,lonl))
        for t in range(timel):
            t_start = (t*bin_num)     #determines the month range to put in each bin
            t_end = t_start + bin_num
            if t_end > timefull:   #if there aren't enough months to fill last bin, loop around to start of year
                t_end -= timefull
                for x in range(lonl):
                    for y in range(latl):
                        tas_bin[t,y,x] = (sum(tas_ar[t_start:,y,x]) + sum(tas_ar[:t_end,y,x])) / bin_num
                        pr_bin[t,y,x] = (sum(pr_ar[t_start:,y,x]) + sum(pr_ar[:t_end,y,x])) / bin_num
                        lsm_bin[t,y,x] = (sum(lsm_ar[t_start:,y,x]) + sum(lsm_ar[:t_end,y,x])) / bin_num
                        sic_bin[t,y,x] = (sum(sic_ar[t_start:,y,x]) + sum(sic_ar[:t_end,y,x])) / bin_num
                        if sum_def == 1:
                            czen_bin[t,y,x] = (sum(czen_ar[t_start:,y,x]) + sum(czen_ar[:t_end,y,x])) / bin_num
                        if use_ts == 1:
                            ts_bin[t,y,x] = (sum(ts_ar[t_start:,y,x]) + sum(ts_ar[:t_end,y,x])) / bin_num
            else:
                for x in range(lonl):
                    for y in range(latl):
                        tas_bin[t,y,x] = stat.mean(tas_ar[t_start:t_end,y,x])     #fill each bin with averages from above-specified month range
                        pr_bin[t,y,x] = stat.mean(pr_ar[t_start:t_end,y,x])
                        lsm_bin[t,y,x] = stat.mean(lsm_ar[t_start:t_end,y,x])
                        sic_bin[t,y,x] = stat.mean(sic_ar[t_start:t_end,y,x])
                        if sum_def == 1:
                            czen_bin[t,y,x] = stat.mean(czen_ar[t_start:t_end,y,x])
                        if use_ts == 1:
                            ts_bin[t,y,x] = stat.mean(ts_ar[t_start:t_end,y,x])
    else:   #if months aren't binned, just rename the data arrays for use in the next section
        tas_bin = tas_ar
        pr_bin = pr_ar
        lsm_bin = lsm_ar
        sic_bin = sic_ar
        if sum_def == 1:
            czen_bin = czen_ar
        if use_ts == 1:
            ts_bin = ts_ar

    #Interpolation routine. 
    if interp > 0:
        print('Interpolating Data to Higher Resolution...')
        if dum_ice == 1:
            print(' Adding "dummy" sea ice...')
            for x in range(lonl):
                for y in range(latl):
                    if lsm_bin[0,y,x] < 0.5:
                        keep = 0
                    else:
                        if y > 0 and lsm_bin[0,y-1,x] < 0.5:
                            checka = 1
                            keep = 1
                        else:
                            checka = 0
                            keep = 0
                        if y < latl-1 and lsm_bin[0,y+1,x] < 0.5:
                            checkb = 1
                            keep = 1
                        else:
                            checkb = 0
                            keep = max(keep,0)
                        if lsm_bin[0,y,x-1] < 0.5:
                            checkc = 1
                            keep = 1
                        else:
                            checkc = 0
                            keep = max(keep,0)
                        if (x < lonl-1 and lsm_bin[0,y,x+1] < 0.5) or (x == lonl-1 and lsm_bin[0,y,0] < 0.5):
                            checkd = 1
                            keep = 1
                        else:
                            checkd = 0
                            keep = max(keep,0)
                    if keep == 0:
                        pass
                    elif keep == 1:
                        for t in range(timel):
                            if checka == 1:
                                sic_bin[t,y,x] = sic_bin[t,y-1,x]
                            if checkb == 1:
                                sic_bin[t,y,x] = max(sic_bin[t,y,x],sic_bin[t,y+1,x])
                            if checkc == 1:
                                sic_bin[t,y,x] = max(sic_bin[t,y,x],sic_bin[t,y,x-1])
                            if checkd == 1:
                                if x == lonl-1:
                                    sic_bin[t,y,x] = max(sic_bin[t,y,x],sic_bin[t,y,0])
                                else:
                                    sic_bin[t,y,x] = max(sic_bin[t,y,x],sic_bin[t,y,x+1])
                                    
        if use_topo == 1:
            print(' Estimating atmospheric lapse rate...')
            lapse_first = np.empty((timel,latl,lonl))
            lapse_bin = np.empty((timel,latl,lonl))
            for t in range(timel):
                lapses = []
                for y in range(latl):
                    for x in range(lonl):   #first, find appropriate neighbor cells
                        n = 0
                        nsum = 0
                        if x == 0:
                            l = lonl-1
                        else:
                            l = x-1
                        if x == lonl-1:
                            r = 0
                        else:
                            r = x+1
                        if y == 0:
                            neigh = [(l,y), (l,y+1), (x,y+1), (r,y+1), (r,y)]
                        elif y == latl-1: 
                            neigh = [(l,y), (l,y-1), (x,y-1), (r,y-1), (r,y)]
                        else:
                            neigh= [(l,y), (l,y+1), (x,y+1), (r,y+1), (r,y), (l,y-1), (x,y-1), (r,y-1)]
                        for ne in neigh:
                            if abs(grnz_bin[ne[1],ne[0]] - grnz_bin[y,x]) > 981:   #checks that there's at least 100 meter equivalent difference in elevation to make the comparison
                                n += 1
                                nsum += (tas_bin[t,ne[1],ne[0]] - tas_bin[t,y,x]) / (grnz_bin[ne[1],ne[0]] - grnz_bin[y,x])
                        if n == 0:
                            lapse_first[t,y,x] = 10000  #obviously unphysical placeholder if there isn't enough data
                        else:
                            lapse_first[t,y,x] = nsum / n   #otherwise put in average
                            lapses.append(lapse_first[t,y,x])
                tot_avg = stat.mean(lapses)
                for y in range(latl):
                    for x in range(lonl):
                        n = 1
                        if lapse_first[t,y,x] < 1000:
                            nsum = lapse_first[t,y,x]
                        else:
                            nsum = tot_avg
                        if x == 0:
                            l = lonl-1
                        else:
                            l = x-1
                        if x == lonl-1:
                            r = 0
                        else:
                            r = x+1
                        if y == 0:
                            neigh = [(l,y), (l,y+1), (x,y+1), (r,y+1), (r,y)]
                        elif y == latl-1: 
                            neigh = [(l,y), (l,y-1), (x,y-1), (r,y-1), (r,y)]
                        else:
                            neigh= [(l,y), (l,y+1), (x,y+1), (r,y+1), (r,y), (l,y-1), (x,y-1), (r,y-1)]
                        for ne in neigh:
                            if lapse_first[t,ne[1],ne[0]] < 1000:
                                n += 1
                                nsum += lapse_first[t,ne[1],ne[0]]
                        lapse_bin[t,y,x] = nsum / n
                        
                        
                        
            
        print(' Reticulating splines...')
        tas_interp = np.empty((timel,latl*interp,lonl*interp))     #create our upscaled data arrays
        pr_interp = np.empty((timel,latl*interp,lonl*interp))
        lsm_interp = np.empty((timel,latl*interp,lonl*interp))
        sic_interp = np.empty((timel,latl*interp,lonl*interp))
        if use_topo == 1:
            grnz_interp = np.empty((latl*interp,lonl*interp))
            lapse_interp = np.empty((timel,latl*interp,lonl*interp))
        if sum_def == 1:
            czen_interp = np.empty((timel,latl*interp,lonl*interp))
        if use_ts == 1:
            ts_interp = np.empty((timel,latl*interp,lonl*interp))
        lat_in = np.empty((latl))
        lon_in = np.empty((lonl))
        for i in range(latl):
            lat_in[i] = (90-lat[i])*np.pi/180
        for i in range(lonl):
            lon_in[i] = lon[i]*np.pi/180
        lat_out = np.linspace(180/(latl*interp),180,latl*interp,False)*np.pi/180
        lon_out = np.linspace(0,360,lonl*interp,False)*np.pi/180
        lat_out,lon_out = np.meshgrid(lat_out,lon_out)
        
        if use_topo == 1:
            grnz_spline = spint(lat_in,lon_in,grnz_bin)
            grnz_out = grnz_spline.ev(lat_out.ravel(),lon_out.ravel()).reshape((lonl*interp,latl*interp)).T
            grnz_interp[:,:] = grnz_out[:,:]
        
        for t in range(timel):
            tas_in = tas_bin[t,:,:]
            pr_in = pr_bin[t,:,:]
            sic_in = sic_bin[t,:,:]
            if use_topo == 1:
                lapse_in = lapse_bin[t,:,:]
            if sum_def == 1:
                czen_in = sic_bin[t,:,:]
            if use_ts == 1:
                ts_in = tas_bin[t,:,:]
            tas_spline = spint(lat_in,lon_in,tas_in)
            tas_out = tas_spline.ev(lat_out.ravel(),lon_out.ravel()).reshape((lonl*interp,latl*interp)).T
            tas_interp[t,:,:] = tas_out[:,:]
            pr_spline = spint(lat_in,lon_in,pr_in)
            pr_out = pr_spline.ev(lat_out.ravel(),lon_out.ravel()).reshape((lonl*interp,latl*interp)).T
            pr_interp[t,:,:] = pr_out[:,:]
            sic_spline = spint(lat_in,lon_in,sic_in)
            sic_out = sic_spline.ev(lat_out.ravel(),lon_out.ravel()).reshape((lonl*interp,latl*interp)).T
            sic_interp[t,:,:] = sic_out[:,:]
            if use_topo == 1:
                lapse_spline = spint(lat_in,lon_in,lapse_in)
                lapse_out = lapse_spline.ev(lat_out.ravel(),lon_out.ravel()).reshape((lonl*interp,latl*interp)).T
                lapse_interp[t,:,:] = lapse_out[:,:]
            if sum_def == 1:
                czen_spline = spint(lat_in,lon_in,czen_in)
                czen_out = czen_spline.ev(lat_out.ravel(),lon_out.ravel()).reshape((lonl*interp,latl*interp)).T
                czen_interp[t,:,:] = czen_out[:,:]
            if use_ts == 1:
                ts_spline = spint(lat_in,lon_in,ts_in)
                ts_out = ts_spline.ev(lat_out.ravel(),lon_out.ravel()).reshape((lonl*interp,latl*interp)).T
                ts_interp[t,:,:] = ts_out[:,:]
            for x in range(lonl*interp):
                for y in range(latl*interp):
                    y_l = math.floor(y/interp)
                    x_l = math.floor(x/interp)
                    lsm_interp[t,y,x] = lsm_bin[t,y_l,x_l]
        
        latl *= interp
        lonl *= interp
        tas = tas_interp      #at the end of all interpolation, rename the data arrays to the final name used in interpretation
        pr = pr_interp
        sic = sic_interp
        lsm = lsm_interp
        if sum_def == 1:
            czen = czen_interp
        if use_ts == 1:
            ts = ts_interp
        if use_topo == 1:
            grnz = grnz_interp
            lapse = lapse_interp
            
            print(' Uploading Higher-Resolution Topogaphy...')
            hmap = Image.open(topo_path)
            if blend_topo == 1:     #Use topo map to create new lsm
                hmap_bin = hmap.convert('L')    #Convert to 8-bit to avoid floating-point issues
                hmap_binext = hmap_bin.getextrema()
                thresh = (sealev - minel) * (hmap_binext[1] - hmap_binext[0]) / (maxel - minel) #Find greyscale value corresponding to sea level
                hmap_bin = hmap_bin.point( lambda p: 255 if p > thresh else 0 ) #Binarize to land/sea mask
                hmap_bin = hmap_bin.resize((lonl,latl))
                hmap_bin = hmap_bin.point( lambda p: 1 if p > 127 else 0 )  #Re-binarize after downscaling
                topo_lsm = np.asarray(hmap_bin)
                lsm = np.ones(lsm.shape) * topo_lsm[np.newaxis,:,:]    #Copy onto array with same shape as lsm
            hmap = hmap.convert('F')    #Convert to float for precision on downscaling
            hmapext = hmap.getextrema()
            hmap = hmap.resize((lonl,latl))
            elev = np.asarray(hmap)
            elev = elev * (maxel - minel) / (hmapext[1] - hmapext[0]) * gravity #Scale to elevation range
            elev = elev[np.newaxis,:,:] #Add time axis for appropriate broadcasting
            grnz = grnz[np.newaxis,:,:]
            print(' Adjusting Temperature by Topography...')
            tas = tas + lapse * (elev - grnz)
            if use_ts == 1:
                ts = ts - lapse * (grnz - sealev * gravity)
                
                
    else:   #if no interpolation is used, just rename the data arrays for use in next section
        tas = tas_bin
        pr = pr_bin
        sic = sic_bin
        lsm = lsm_bin
        if sum_def == 1:
            czen = czen_bin
        if use_ts == 1:
            ts = ts_bin


    print('Interpreting Data...')
    #create our output arrays
    Koppen_Full = np.empty((lonl,latl),dtype=np.int8)
    Sea_Zones = np.empty((lonl,latl),dtype=np.int8)

    #Extract appropriate averages and extremes with numpy operations
    Avg_Temp_ar = np.mean(tas, axis=0)-273.15 #Converts to Celsius
    Total_Precip_ar = np.mean(pr, axis=0)*2592000000*12 #Converts to mm for whole year
    Max_Temp_ar = np.amax(tas, axis=0)-273.15
    Min_Temp_ar = np.amin(tas, axis=0)-273.15
    Min_Precip_ar = np.amin(pr, axis=0)*2592000000  #Converts to mm for month
    if sea_def != 1:
        Avg_Ice_ar = np.mean(sic, axis=0)
        Max_Ice_ar = np.amax(sic, axis=0)
        Min_Ice_ar = np.amin(sic, axis=0)
    if use_ts == 1:
        if sea_type == 2 or sea_type == 3:
            Avg_Seatemp_ar = np.mean(ts, axis=0)-273.15
        else:
            Min_Seatemp_ar = np.amin(ts, axis=0)-273.15
            if sea_def == 1:
                Max_Seatemp_ar = np.amax(ts, axis=0)-273.15
    if land_type == 5:       #routine specific to Holdridge Life Zone
        Avg_Biot_ar = np.mean(np.clip(tas-273.15, 0, 30), axis=0) #constrains temperature range to [0,30] C
    elif land_type < 3 or land_type == 4:
        if timodd != 1:
            if sum_def == 1:
                long_ar = np.concatenate((czen,czen[:halfl,:,:]),axis=0)    #create 1.5-length array of zenith angle or temp
            else:
                long_ar = np.concatenate((tas,tas[:halfl,:,:]), axis=0)
            sum_ar = np.sum(np.stack([long_ar[m:m+halfl] for m in range(halfl)], axis=0), axis=0)   #stack half-year slices together and then sum them together to produce array of total values over following half-year for each month
            pr_long = np.concatenate((pr,pr),axis=0)    #create double-length array of precipitation
            pr_stack = np.stack([pr_long[m:m+timel] for m in range(timel)], axis=0)     #stack full-year slices together starting on each month of the year
            sum_max = np.argmax(sum_ar, axis=0)
            precips_ar = np.empty((timel,latl,lonl))
            for y in range(latl):       #Bit of an inefficient approach, but haven't quite figured out how to do this as a numpy operation
                for x in range(lonl):
                    precips_ar[:,y,x] = pr_stack[sum_max[y,x],:,y,x]    #choose appropriate slice corresponding to the maximum half-year total of indicator value            
            Summer_Precip_ar = np.mean(precips_ar[:halfl,:,:], axis=0)*2592000000*6
            Max_Sum_Precip_ar = np.amax(precips_ar[:halfl,:,:], axis=0)*2592000000
            Min_Sum_Precip_ar = np.amin(precips_ar[:halfl,:,:], axis=0)*2592000000
            Max_Win_Precip_ar = np.amax(precips_ar[halfl:,:,:], axis=0)*2592000000
            Min_Win_Precip_ar = np.amin(precips_ar[halfl:,:,:], axis=0)*2592000000
    
        else:   #For the moment, odd numbers of months requires reverting to the old method of iterating over arrays to find appropriate seasonal values
            Summer_Precip_ar = np.empty((lonl,latl))
            Max_Sum_Precip_ar = np.empty((lonl,latl))
            Min_Sum_Precip_ar = np.empty((lonl,latl))
            Max_Win_Precip_ar = np.empty((lonl,latl))
            Min_Win_Precip_ar = np.empty((lonl,latl))
            for y in range(latl):
                for x in range(lonl):
                    if sum_def == 1:
                        zen_ar = czen[:,y,x]
                        zen_li = zen_ar.tolist()
                        zen_li.extend(zen_li[:])
                    else:
                        Temps_ar = tas[:,y,x]
                        Temps_li = Temps_ar.tolist()
                        Temps_li.extend(Temps_li[:])
                    Precips_ar = pr[:,y,x]
                    Precips_li = Precips_ar.tolist()
                    Precips_li.extend(Precips_li[:])
                    #Determines start of summer
                    for m in range(timel):
                        if sum_def == 1:
                            avg = stat.mean(zen_li[m:m+halfl])
                            if zen_li[m-1] > zen_li[m+halfl]:
                                Sum_add = 0
                                avg = (avg * halfl + zen_li[m-1])/(halfl+0.5)
                            else:
                                Sum_add = 1
                                avg = (avg * halfl + zen_li[m+halfl])/(halfl+0.5)
                        else:
                            avg = stat.mean(Temps_li[m:m+halfl])    #for each month, average the temps for the following half a year
                            if Temps_li[m-1] > Temps_li[m+halfl]:
                                Sum_add = 0
                                avg = (avg * halfl + Temps_li[m-1])/(halfl+0.5)
                            else:
                                Sum_add = 1
                                avg = (avg * halfl + Temps_li[m+halfl])/(halfl+0.5)
                        if m == 0:      #start by assuming the first month is the start of summer
                            max_avg = avg
                            Sum_start = 0
                            Sum_length = 0
                        elif avg > max_avg:     #whenever the above process gives a higher average temp than the previous highest, make that the new start of summer
                            Sum_start = m
                            max_avg = avg
                    
                    #Extracts relevant seasonal precipitation data
                    Summer_Precip_ar[y,x] = stat.mean(Precips_li[Sum_start:Sum_start+halfl])*2592000000*6
                    Max_Sum_Precip_ar[y,x] = max(Precips_li[Sum_start:Sum_start+halfl])*2592000000
                    Min_Sum_Precip_ar[y,x] = min(Precips_li[Sum_start:Sum_start+halfl])*2592000000
                    if Sum_add == 1:
                        Max_Win_Precip_ar[y,x] = max(Precips_li[Sum_start+halfl+1:Sum_start+timel])*2592000000
                        Min_Win_Precip_ar[y,x] = min(Precips_li[Sum_start+halfl+1:Sum_start+timel])*2592000000
                    else:
                        Max_Win_Precip_ar[y,x] = max(Precips_li[Sum_start+halfl:Sum_start+timel])*2592000000
                        Min_Win_Precip_ar[y,x] = min(Precips_li[Sum_start+halfl:Sum_start+timel])*2592000000
                    if Sum_add == 0:
                        Summer_Precip_ar[y,x] = Summer_Precip_ar[y,x] + Precips_li[Sum_start-1]*2592000000*6/timel
                    else:
                        Summer_Precip_ar[y,x] = Summer_Precip_ar[y,x] + Precips_li[Sum_start+halfl]*2592000000*6/timel
            
        Summer_Length_ar = np.sum(np.where(tas>283.15,1,0),axis=0) #add together months above 10 C to determine "summer length" used for some zones

    #Determines climate zones from data
    #Identifying both land and sea zones for every cell is perhaps inefficient for blended outputs, but I presume that anyone doing a very high-resolution, interpolated output is probably going to want unblended outputs anyway
    print('Determining Climate Zones...')
    for y in range(latl):
        for x in range(lonl):
            clim = 0    #makes debug easier
            
            #find data for specific location
            Avg_Temp = Avg_Temp_ar[y,x]
            Max_Temp = Max_Temp_ar[y,x]
            Min_Temp = Min_Temp_ar[y,x]
            Total_Precip = Total_Precip_ar[y,x]
            Min_Precip = Min_Precip_ar[y,x]
            if land_type == 5:
                Avg_Biot = Avg_Biot_ar[y,x]
            elif land_type < 3 or land_type == 4:
                Summer_Precip = Summer_Precip_ar[y,x]
                Max_Sum_Precip = Max_Sum_Precip_ar[y,x]
                Min_Sum_Precip = Min_Sum_Precip_ar[y,x]
                Max_Win_Precip = Max_Win_Precip_ar[y,x]
                Min_Win_Precip = Min_Win_Precip_ar[y,x]
                Summer_Length = Summer_Length_ar[y,x]
            if sea_def != 1:    #only pull up sea ice data if we're actually marking sea ice
                Max_Ice = Max_Ice_ar[y,x]
                Min_Ice = Min_Ice_ar[y,x]
                Avg_Ice = Avg_Ice_ar[y,x]
            if use_ts == 1:
                if sea_type == 2 or sea_type == 3:
                    Avg_Seatemp = Avg_Seatemp_ar[y,x]
                else:
                    Min_Seatemp = Min_Seatemp_ar[y,x]
                    if sea_def == 1:
                        Max_Seatemp = Max_Seatemp_ar[y,x]
            else:
                if sea_type == 2 or sea_type == 3:
                    Avg_Seatemp = Avg_Temp
                else:
                    Min_Seatemp = Min_Temp
                    if sea_def == 1:
                        Max_Seatemp = Max_Temp
            #land zones
            if blend == 1 or lsm[0,y,x] >= 0.5:   #Only check for land zones if necessary
                
                #Groups only
                
                if land_type == 4:
                
                    #determine threshold for arid zones
                    if Summer_Precip > Total_Precip*0.7:
                        adjust = 280
                    elif Summer_Precip > Total_Precip*0.3:
                        adjust = 140
                    else:
                        adjust = 0
                    Arid_threshold = Avg_Temp*20 + adjust
                    
                    #Find Koppen Group
                    if Total_Precip < Arid_threshold and (Max_Temp > 10 or (Max_Temp > 0 and ice_def != 2) or (ice_def != 1 and ice_def != 2)):  #the and/or functions are accounting for different definitions set at the start
                        clim = B
                    elif Max_Temp < 10:
                        clim = E
                    elif Min_Temp > 18:
                        clim = A
                    elif Min_Temp > 0 or (temp_def == 1 and Avg_Temp > -3):
                        clim = C
                    else:
                        clim = D
                
                #seasonless zones
                
                elif land_type == 3:
                    if Total_Precip < Avg_Temp*20+140 and (Avg_Temp > 10 or (Avg_Temp > 0 and ice_def != 2) or (ice_def != 1 and ice_def != 2)):
                        if Total_Precip < Avg_Temp*10+70:
                            if Avg_Temp > 18 or (arid_def != 1 and Avg_Temp > 0) or (arid_def != 1 and temp_def == 1 and Avg_Temp > -3):
                                clim = HotDesertnos
                            else:
                                clim = ColdDesertnos
                        else:
                            if Avg_Temp > 18 or (arid_def != 1 and Avg_Temp > 0) or (arid_def != 1 and temp_def == 1 and Avg_Temp > -3):
                                clim = HotSteppenos
                            else:
                                clim = ColdSteppenos
                    elif Avg_Temp < 10:
                        if Avg_Temp < 0:
                            clim = IceCapnos
                        else:
                            clim = Tundranos
                    elif Avg_Temp > 18:
                        if Total_Precip/12 > 60:
                            clim = TropRainforestnos
                        else:
                            clim = TropSavannanos
                    else:
                        clim = Oceanicnos
                        
                #Full Koppen set
                
                elif land_type < 3:
                    #determine threshold for arid zones
                    if Summer_Precip > Total_Precip*0.7:
                        adjust = 280
                    elif Summer_Precip > Total_Precip*0.3:
                        adjust = 140
                    else:
                        adjust = 0
                    Arid_threshold = Avg_Temp*20 + adjust
                    
                    #Finds Koppen zone
                    if Total_Precip < Arid_threshold and (Max_Temp > 10 or (Max_Temp > 0 and ice_def != 2) or (ice_def != 1 and ice_def != 2)):
                        if Total_Precip < Arid_threshold/2:
                            if (arid_def != 1 and Min_Temp > 0) or (arid_def == 1 and Avg_Temp > 18) or (arid_def != 1 and temp_def == 1 and Min_Temp > -3):
                                clim = BWh
                            else:
                                clim = BWk
                        else:
                            if (arid_def != 1 and Min_Temp > 0) or (arid_def == 1 and Avg_Temp > 18) or (arid_def != 1 and temp_def == 1 and Min_Temp > -3):
                                clim = BSh
                            else:
                                clim = BSk
                    elif Max_Temp < 10:
                        if Max_Temp < 0:
                            clim = EF
                        else:
                            clim = ET
                    elif Min_Temp > 18:
                        if Min_Precip > 60:
                            clim = Af
                        elif Min_Precip > 100-Total_Precip/25:
                            clim = Am
                        else:
                            if Summer_Precip < Total_Precip/2:
                                clim = As
                            else:
                                clim = Aw
                    elif Min_Temp > 0 or (temp_def == 1 and Avg_Temp > -3):
                        if (Min_Sum_Precip < 30 or (med_def == 1 and Min_Sum_Precip < 40)) and Min_Win_Precip > Min_Sum_Precip and Max_Win_Precip > Min_Sum_Precip*3 and ((add_def == 1 or add_def == 3) or ((add_def == 0 or add_def == 2) and Summer_Precip < Total_Precip*0.5)) and (prio_def == 0 or (prio_def == 1 and ((wet_def != 1 and Max_Sum_Precip <= Min_Win_Precip*10) or (wet_def == 1 and Summer_Precip <= Total_Precip*0.7)))):
                            if Summer_Length < thirdl:
                                clim = Csc
                            elif Max_Temp > 22:
                                clim = Csa
                            else:
                                clim = Csb
                        elif ((wet_def != 1 and Max_Sum_Precip > Min_Win_Precip*10) or (wet_def == 1 and Summer_Precip > Total_Precip*0.7)) and ((add_def == 0 or add_def == 3) or ((add_def == 1 or add_deff == 2) and Summer_Precip > Total_Precip*0.5)):
                            if Summer_Length < thirdl:
                                clim = Cwc
                            elif Max_Temp > 22:
                                clim = Cwa
                            else:
                                clim = Cwb
                        else:
                            if Summer_Length < thirdl:
                                clim = Cfc
                            elif Max_Temp > 22:
                                clim = Cfa
                            else:
                                clim = Cfb
                    else:
                        if (Min_Sum_Precip < 30 or (med_def == 1 and Min_Sum_Precip < 40)) and Min_Win_Precip > Min_Sum_Precip and Max_Win_Precip > Min_Sum_Precip*3 and ((add_def == 1 or add_def == 3) or ((add_def == 0 or add_def == 2) and Summer_Precip < Total_Precip*0.5)) and (prio_def == 0 or (prio_def == 1 and ((wet_def != 1 and Max_Sum_Precip <= Min_Win_Precip*10) or (wet_def == 1 and Summer_Precip <= Total_Precip*0.7)))):
                            if Summer_Length < thirdl:
                                if Min_Temp < -38:
                                    clim = Dsd
                                else:
                                    clim = Dsc
                            elif Max_Temp > 22:
                                clim = Dsa
                            else:
                                clim = Dsb
                        elif ((wet_def != 1 and Max_Sum_Precip > Min_Win_Precip*10) or (wet_def == 1 and Summer_Precip > Total_Precip*0.7)) and ((add_def == 0 or add_def == 3) or ((add_def == 1 or add_deff == 2) and Summer_Precip > Total_Precip*0.5)):
                            if Summer_Length < thirdl:
                                if Min_Temp < -38:
                                    clim = Dwd
                                else:
                                    clim = Dwc
                            elif Max_Temp > 22:
                                clim = Dwa
                            else:
                                clim = Dwb
                        else:
                            if Summer_Length < thirdl:
                                if Min_Temp < -38:
                                    clim = Dfd
                                else:
                                    clim = Dfc
                            elif Max_Temp > 22:
                                clim = Dfa
                            else:
                                clim = Dfb
                
                #Reduced sets; determined from the above rather than having a new routine for every option
                
                if land_type == 1:
                    if clim == As:
                        clim = Aw
                
                if land_type == 2:
                    if clim == Af:
                        clim = TropRainforest
                    elif clim == Am:
                        clim = TropMonsoon
                    elif clim == As or clim == Aw:
                        clim = TropSavanna
                    elif clim == BWh:
                        clim = HotDesert
                    elif clim == BWk:
                        clim = ColdDesert
                    elif clim == BSh:
                        clim = HotSteppe
                    elif clim == BSk:
                        clim = ColdSteppe
                    elif clim == Csa or clim == Csb or clim == Csc:
                        clim = Med
                    elif clim == Cwa or clim == Cfa:
                        clim = Subtropical
                    elif clim == Cwb or clim == Cwc or clim == Cfb or clim == Cfc:
                        clim = Oceanic
                    elif clim == Dsa or clim == Dsb or clim == Dwa or clim == Dwb or clim == Dfa or clim == Dfb:
                        clim = Continental
                    elif clim == Dsc or clim == Dsd or clim == Dwc or clim == Dwd or clim == Dfc or clim == Dfd:
                        clim = Subarctic
                    elif clim == ET:
                        clim = Tundra
                    elif clim == EF:
                        clim = IceCap
                
                #Holdridge Life Zones
                
                if land_type == 5:
                    if Avg_Biot > 24:
                        if Total_Precip > 8000:
                            clim = HTropRainForest
                        elif Total_Precip > 4000:
                            clim = HTropWetForest
                        elif Total_Precip > 2000:
                            clim = HTropMoistForest
                        elif Total_Precip > 1000:
                            clim = HTropDryForest
                        elif Total_Precip > 500:
                            clim = HVDryForest
                        elif Total_Precip > 250:
                            clim = HThornWood
                        elif Total_Precip > 125:
                            clim = HTropDesertScrub
                        else:
                            clim = HTropDesert
                    elif Avg_Biot > 12:
                        if Total_Precip > 4000:
                            clim = HWarmRainForest
                        elif Total_Precip > 2000:
                            clim = HWarmWetForest
                        elif Total_Precip > 1000:
                            clim = HWarmMoistForest
                        elif Total_Precip > 500:
                            clim = HWarmDryForest
                        elif Total_Precip > 250:
                            clim = HThornSteppe
                        elif Total_Precip > 125:
                            clim = HWarmDesertScrub
                        else:
                            clim = HWarmDesert
                    elif Avg_Biot > 6:
                        if Total_Precip > 2000:
                            clim = HCoolRainForest
                        elif Total_Precip > 1000:
                            clim = HCoolWetForest
                        elif Total_Precip > 500:
                            clim = HCoolMoistForest
                        elif Total_Precip > 250:
                            clim = HSteppe
                        elif Total_Precip > 125:
                            clim = HCoolDesertScrub
                        else:
                            clim = HCoolDesert
                    elif Avg_Biot > 3:
                        if Total_Precip > 1000:
                            clim = HBorRainForest
                        elif Total_Precip > 500:
                            clim = HBorWetForest
                        elif Total_Precip > 250:
                            clim = HBorMoistForest
                        elif Total_Precip > 125:
                            clim = HDryScrub
                        else:
                            clim = HBorDesert
                    elif Avg_Biot > 1.5:
                        if Total_Precip > 500:
                            clim = HRainTundra
                        elif Total_Precip > 250:
                            clim = HWetTundra
                        elif Total_Precip > 125:
                            clim = HMoistTundra
                        else:
                            clim = HDryTundra
                    else:
                        clim = HPolDesert
                
                #Write to array
                Koppen_Full[x,y] = clim
            
            #sea zones
            if blend == 1 or lsm[0,y,x] < 0.5:
            
                #no seas
                
                if sea_type == 5:
                    pass
                else:
                    
                    #flat seas
                    
                    if sea_type == 4:
                        seaclim = SeaFlat
                    
                    #seasonless seas
                    
                    elif sea_type == 2 or sea_type == 3:
                        if Avg_Seatemp > 18:
                            seaclim = SeaTropnos
                        elif (sea_def != 1 and Avg_Ice > 0.5) or (sea_def == 1 and Avg_Seatemp < -2):
                                seaclim = SeaPermIcenos
                        else:
                            seaclim = SeaTempnos
                    
                    #seasonal seas
                    
                    else:
                        if Min_Seatemp > 18:
                            seaclim = SeaTrop
                        elif (sea_def != 1 and Min_Ice > 0.8) or (sea_def == 1 and Max_Seatemp < -2):
                            seaclim = SeaPermIce
                        elif (sea_def != 1 and Max_Ice > 0.2) or (sea_def == 1 and Min_Seatemp < -2):
                            seaclim = SeaSeasonalIce
                        else:
                            seaclim = SeaTemp
                    
                    #reduced sets
                    
                    if sea_type == 1:
                        if seaclim == SeaTrop:
                            seaclim = SeaTemp
                    
                    if sea_type == 3:
                        if seaclim == SeaTropnos:
                            seaclim = SeaTempnos
                    
                    #Write to array
                    Sea_Zones[x,y] = seaclim


    #Applies the coloring function to our climate zone arrays and draws images from them
    print('Drawing Climate Maps...')
    colmap = make_colmap(color_type, col_list_path)
    if blend == 1:
        Kop_rgb = color(Koppen_Full, colmap, latl, lonl, lsm, blend, land_type, sea_type, color_type, col_list_path)
        landmap = Image.fromarray(Kop_rgb)
        savename = output_name + '_land.png'
        landmap.save(savename)
        if sea_type != 5:
            Sea_rgb = color(Sea_Zones, colmap, latl, lonl, lsm, blend, land_type, sea_type, color_type, col_list_path)
            seamap = Image.fromarray(Sea_rgb)
            savename = output_name + '_sea.png'
            seamap.save(savename)
    else:
        blend_rgb = color(Koppen_Full, colmap, latl, lonl, lsm, blend, land_type, sea_type, color_type, col_list_path, inarraysea=Sea_Zones)
        blendmap = Image.fromarray(blend_rgb)
        savename = output_name + '.png'
        blendmap.save(savename)

    #Announce operation complete and point to output files
    print("Operation Complete")
    if blend == 1:
        if sea_type == 5:
            print ('Map Written to ' + str(output_name) + '_land.png')
        else:
            print ('Maps Written to ' + str(output_name) + '_land.png')
            print (' and ' + str(output_name) + '_sea.png')
    else:
        print('Map Written to ' + str(output_name) + '.png')
    return

if __name__ == "__main__":
    #opening info
    print('''
NetCDF to Koppen Climate Zones
Script written 2021 by Nikolai Hersfeldt
    of worldbuildingpasta.blogspot.com
For use with NetCDF output files from ExoPlaSim GCM

''')

    #set path to current location
    path = os.path.join(os.path.dirname(__file__), '')
    in_files = []
    in_num = 0

    #ask for input file
    while True:
        infile = path+input('Input NetCDF filename or folder of files: ')
        if os.path.exists(infile):
            if os.path.isdir(infile):
                in_num = 0
                found = False
                for f in os.listdir(infile):
                    if f.endswith(".nc"):
                        in_files.append(infile+"/"+f)
                        in_num += 1
                        found = True
                        print(" Found "+str(f))
                if found:
                    print("Found all files")
                    break
                else:
                    print("No files found in "+str(infile))
            else:
                in_files.append(infile)
                in_num = 0
                print('File found')
                break
        print('No file found at '+str(infile))
        


    #Prompts for configuration
    #The yes/no prompts are intentionally flexible; any input containing 'y' or '1' will be read as 'yes', anything else will be read as 'no'
    res = input('Advanced setup? (y/n): ')
    if res in ('y') or res in ('1'):
        res = 1
        extra = input('''
Additional Input Files
If added, values across all years will be averaged together.
CAUTION: All files must have same resolution and number of months.
Add additional inputs? (y/n): ''')
        if extra in ('y') or extra in ('1'):
            if len(in_files) < 1:
                in_files.append(infile)
                in_num += 1
            while True:
                nextin = path+input('Add input file ("stop" for no more input files): ')
                if nextin in ('stop'):
                    break
                if os.path.exists(nextin):
                    in_files.append(nextin)
                    in_num += 1
                else:
                    print('No file found at '+str(nextin))
            print('All files found')

    #Establish all inputs with default values first
    output_name="output"
    cfg_loadname=""
    land_type=0
    sea_type=0
    blend=0
    color_type=0
    col_list_path="none"
    bin_num=1
    interp=0
    dum_ice=0
    use_topo=0
    topo_path="none"
    maxel=0
    minel=0
    gravity=0
    blend_topo=0
    sealev=0
    sum_def=0
    temp_def=0
    arid_def=0
    med_def=0
    wet_def=0
    add_def=0
    prio_def=0
    ice_def=0
    sea_def=0

    #Main loop, so we can use the same files for multiple outputs:
    
    while True:
        if res == 1:
            cfg_op = input('''
Load config file? (y/n): ''')
            if cfg_op in ('y') or cfg_op in ('1'):
                cfg_op = 1
                while True:
                    cfg_loadname = path+input('Config filename: ')
                    if os.path.exists(cfg_loadname):
                        break
                    else:
                        print('No file found at '+str(cfg_loadname))
                print('Config found')
                cfg = configparser.ConfigParser()
                cfg.read(cfg_loadname)
                typ = cfg['Zone Types']
                land_type = typ.getint('Land Climate Zones')
                sea_type = typ.getint('Sea Climate Zones')
                blend = typ.getint('Land/Sea Blend')
                color_type = typ.getint('Color List')
                col_list_path = typ.get('Custom Color List')
                if color_type == 3:
                    while True:
                        if os.path.exists(col_list_path):
                            break
                        print('No Color List found at '+str(col_list_path))
                        col_list = input('New Color List (input "0" to choose from default options): ')
                        if int(col_list) == 0:
                            color_type = int(input('''0: Default (blue rainforest); used by Wikipedia, Artifexian, etc.
1: Alternate (red rainforest); used in some research papers, usually with older definitions (uses same sea zone colors, as those are my invention)
2: "True" color; fills each zone with their average color on Earth, based on satellite imagery (all sea zones except permanent ice given same color)
Set Color List: '''))
                        else:
                            col_list_path = path+col_list
                            break
                adv = cfg['Advanced Functions']
                bin_num = adv.getint('Bin Size')
                interp = adv.getint('Interpolation Factor')
                dum_ice = adv.getint('Dummy Sea Ice')
                use_topo = adv.getint('Adjust Temp by Topography')
                topo_path = adv.get('Topography Map')
                maxel = adv.getfloat('Max Elevation')
                minel = adv.getfloat('Min Elevation')
                gravity = adv.getfloat('Surface Gravity')
                blend_topo = adv.getint('Blend by Topography Map')
                sealev = adv.getfloat('Sea Level')
                if use_topo == 1:
                    if topo_path == 'prompt':
                        topo_map = input('New Topography Map for temperature adjustmet (input "0" to not use topography): ')
                        if str(topo_map) == "0":
                            use_topo = 0
                            break
                        else:
                            topo_path = path+topo_map
                            use_topo = 2
                    while True:
                        if os.path.exists(topo_path):
                            break
                        use_topo = 2
                        print('No Topography Map found at '+str(topo_path))
                        topo_map = input('New Topography Map for temperature adjustmet (input "0" to not use topography): ')
                        if str(topo_map) == "0":
                            use_topo = 0
                            break
                        else:
                            topo_path = path+topo_map
                    if use_topo == 2:
                        use_topo = 1
                        print('Topography Map found')
                        maxel = float(input('Highest Map Elevation (m): '))
                        minel = float(input('Lowest Map Elevation (m): '))
                        gravity = float(input('Surface Gravity (m/s^2): '))
                        sealev = float(input('Sea Level (m): '))
                defs = cfg['Climate Zone Definitions']
                sum_def = defs.getint('Summer Definition')
                temp_def = defs.getint('Temperate/Continental Boundary')
                arid_def = defs.getint('Hot/Cold Arid Boundary')
                med_def = defs.getint('Mediterranean Definition')
                wet_def = defs.getint('Wet-Summer Definition')
                add_def = defs.getint('Additional Wet Season Requirements')
                prio_def = defs.getint('Med/Wet-Summer Priority')
                ice_def = defs.getint('Arid/Polar Priority')
                sea_def = defs.getint('Sea Ice Definition')
                print('Config parameters loaded')
        #If a config is loaded, the rest of configuration is skipped. Otherwise:
            else:
                land_type = int(input('''
Land Climate Zones
0: Full Koppen set of 31 climate zones (default)
1: As above, but exclude wet-summer savanna (Aw) / dry-summer savanna (As) distinction (common in formal maps)
2: Reduced set of 14 climate zones (common in worldbuilding circles)
3: Seasonless: Ignore seasonal changes and use annual averages only, producing 9 zones
    (appropriate for very short years or tidal-locked worlds without libration)
4: 5 Koppen Groups (A,B,C,D,E) only
5: Holdridge Life Zones approximation (references only biotemperature and precipitation)
Set Land Climate Zone Type: '''))
                sea_type = int(input('''
Sea Climate Zones
0: Set of 4 zones of my own invention; tropical, temperate, seasonal ice, permanent ice (default)
1: As above but exclude tropical/temperate distinction
2: Seasonless: Mark tropical, temperate, and permanent ice based on annual averages
3: As above but exclude tropical/temperate distinction
4: No sea climate zones, just flat blue
5: Exclude seas entirely, produce no output (appropriate for all-land planets)
Set Sea Climate Type: '''))
                if land_type == 5:
                    color_type = int(input('''
Color List
0: Default; based on the key on Wikipedia's Holdridge Life Zones page.
1: Import custom color list (see defaultcolor.ini for template)
Set Color List: '''))
                    if color_type > 0:
                        color_type = 3
                else:
                    color_type = int(input('''
Color List
0: Default (blue rainforest); used by Wikipedia, Artifexian, etc.
1: Alternate (red rainforest); used in some research papers, usually with older definitions (uses same sea zone colors, as those are my invention)
2: "True" color; fills each zone with their average color on Earth, based on satellite imagery (all sea zones except permanent ice given same color)
3: Import custom color list (see defaultcolor.ini for template)
Set Color List: '''))
                if color_type == 3:
                    while True:
                        col_list = input('Custom Color List filename: ')
                        col_list_path = path + col_list
                        if os.path.exists(col_list_path):
                            break
                        else:
                            print('No file found at '+str(col_list_path))
                    print('Custom Color List found')
                blend = int(input('''
Land/Sea Blend
0: Use land/sea mask in the NetCDF file to blend land and sea climates into a single map:
    tiles marked land show land climate zones, tiles marked sea show sea climate zones (default)
1: Do not blend: Produce distinct land and sea climate maps (appropriate for final output, so you can
    combine them with a high-res land/sea mask yourself)
Set Blend setting: '''))
                bin_num = int(input('''
Bin Size
Bins adjacent "months" in the NetCDF file together.
e.g. if you have a NetCDF file with 36 months, a value of "3" will average together data from every 3 months
    to produce 12 months, which will then be used for determining climate zones.
Ideally there should be a whole number of bins per year.
For an input of "1" or "0" the script will not bin months together.
Set Bin Size: '''))
                #bin_num = 0     #not ready for public use, hasn't worked properly in testing yet.
                if bin_num < 2:
                    bin_num = 1
                interp = int(input('''
Interpolation Rescaling Factor
Multiplies the latitude and longitude of the output and interpolates the climate data up to that resolution
e.g. applying a rescaling factor of 4 to a model with 32x64 resolution will get an output with 128x256 resolution
Interpolation essentially just averages the values of adjacent cells, so isn't as accurate as actual modelling
    and in particular will not reflect the impact of topographical features smaller than those input to ExoPlaSim;
    but it does reflect the underlying data better than trying to interpolate "by eye" from the output maps.
The rescaling factor need not be an integer, but the output must have an integer resolution
    (e.g. you could apply a factor of 2.5 to a 32x64 model to get 80x160 resolution)
CAUTION: High rescaling factors can significantly extend the script's runtime.
Input "0" for no interpolation.
Set Interpolation Rescaling Factor: '''))
                if interp > 0:
                    dum_ice = input('''
"Dummy" Sea Ice Cover
Adds "dummy" sea ice to land whenever there is sea ice in at least 1 adjacent cell, before interpolation
This helps ensure sea ice zones remain contiguous with land when interpolated
Add "dummy" sea ice? (y/n): ''')
                    if dum_ice in ('y') or dum_ice in ('1'):
                        dum_ice = 1
                    else:
                        dum_ice = 0
                    use_topo = input('''
Adjust Temperature to Topography
Allows you to upload a higher-resolution heightmap, and then after interpolation, temperature is adjusted
    based on deviation between that heightmap and an interpolated heightmap from the model output.
Precipitation is not adjusted, so this is still not as accurate as higher-resolution modelling,
    but it is an improvement and gives the map some "texture".
Upload Topography? (y/n): ''')
                    if use_topo in ('y') or use_topo in ('1'):
                        use_topo = 1
                    else:
                        use_topo = 0


        if use_topo == 1:
            if cfg_op == 1:
                pass
            else:
                ds = nc.Dataset(in_files[0])
                lat = ds['lat'][:]
                maph = len(lat) * interp
                mapw = maph * 2
                print('''Input Topography Map
CAUTION: Must be greyscale (white high, black low) and seas should be marked as flat surfaces
Interpolated map resolution:''')
                print(str(maph) + 'x' + str(mapw))
                while True:
                    topo_map = input('Topography Map filename: ')
                    topo_path = path + topo_map
                    if os.path.exists(topo_path):
                        break
                    else:
                        print('No file found at '+str(topo_path))
                print('Topography Map found')
                maxel = float(input('Highest Map Elevation (m): '))
                minel = float(input('Lowest Map Elevation (m): '))
                sealev = float(input('Sea Level (m): '))
                gravity = float(input('Surface Gravity (m/s^2): '))
                if blend != 1:
                    blend_topo = input('''Blend by Topography Map
Blends Land and Sea Maps at Interpolated Resolution using Topography Map
    with specified sea level
Blend by Map? (y/n): ''')
                    if blend_topo in ('y') or blend_topo in ('1'):
                        blend_topo = 1
                    else:
                        blend_topo = 0

        if res == 1:
            if cfg_op == 1:
                pass
            else:
                if land_type != 5:
                    con_def = input('''
Configure Climate Zone Definitions? (y/n): ''')
                else:
                    con_def = 'n'
                if con_def in ('y') or con_def in ('1'):
                    sum_def = int(input('''
Summer Definition
0: Summer is the half of the year with the highest average temperature (default)
1: Summer is the half of the year with the highest average zenith angle (heigh of the sun in the sky;
    more consistent across each hemisphere, but doesn't account well for eccentricity)
Set Definition: '''))
                    temp_def = int(input('''
Temperate/Continental (and Hot/Cold Arid) Boundary
0: Temperate (C) and hot arid zones where [minimum temperate > 0 C], continental (D) and cold arid zones elsewhere (default)
1: C and hot arid zones where [minimum temperature > -3 C], D and cold arid zones elsewhere
    (less common, approximates maximum extent of frozen ground in winter)
Set Definition: '''))
                    arid_def = int(input('''
Hot/Cold Arid Boundary
0: Hot arid (BWh, BSh) zones where [minimum temperature > 0 C] (or -3 C if 1 set in last step),
    cold arid (BWk, BSk) elsewhere (default)
1: Hot arid where [annual average temperature > 18 C], cold arid elsewhere (used in older maps)
Set Definition: '''))
                    if land_type < 3:
                        med_def = int(input('''
Mediterranean Definition
0: Mediterranean (Csa, Csb, Csc, Dsa, Dsb, Dsc, Dsd) zones require driest summer month to have <30 mm of precipitation
    and <1/3 the precipitation of the wettest winter month (default)
1: Mediterranean zones require driest summer month to have <40 mm of precipitation
    and <1/3 the precipitation of the wettest winter month (less common)
Set Definition: '''))
                        wet_def = int(input('''
Wet-Summer Definition
0: Wet-summer (Cwa, Cwb, Cwc, Dwa, Dwb, Dwc Dwd) zones require >10 times more precipitation
    in wettest summer month than driest winter month (default)
1: Wet-summer zones require >70 percent of their precipitation to be in the summer half of the year (less common)
Set Definition: '''))
                        add_def = int(input('''
Additional Wet Season Requirements
0: Mediterranean (Cs, Ds) zones must additionally have more total rain in winter than summer (default)
1: Wet-summer (Cw, Dw) zones must additionally have more total rain in summer than winter
2: Both Med and Wet-summer zones must have more total rain in their wet seasons
3: No additional requirements for Med or Wet-summer zones
Set Definition: '''))
                        if add_def == 2:
                            prio_def = 0
                        else:
                            prio_def = int(input('''
Mediterranean / Wet-Summer Priority
0: Where a region meets the requirements for both, it will be marked as Med (Cs, Ds) rather than Wet-summer (Cw, Dw) (default)
1: Where a region meets the requirements for both, it will be marked as Wet-summer (Cw, Dw) rather than Med (Cs, Ds)
Set Definition: '''))
                    ice_def = int(input('''
Arid/Polar Priority
0: All zones meeting arid definition are classed as arid (B), even if [maximum temperature < 0 C] (default)
1: All zones with [maximum temperature < 0 C] are classed as Ice Cap (EF), even if they meet arid definition
    (less common, may be cleaner for dry worlds with large ice caps)
2: All zones with [maximum temperature < 10 C] are classed as polar (E)
Set Definition: '''))
                    if sea_type < 4:
                        sea_def = int(input('''
Sea Ice Definition
0: Seasonal and permanent sea ice extent determined by checking sea ice cover mask in the NetCDF file (default)
1: Permanent sea ice zone where [maximum temperature < -2 C], seasonal where [minimum temperature < -2 C]
    (my initial definition, before I had access to sea ice modelling)
Set Definition: '''))
                cfg_save_op = input('''
Config
Saves the above configurations and definitions to a .ini file that can be loaded in later uses of this script to save time.
Note that:
    The configs do not include the input or output filenames
    Loading a config skips the rest of configuration
    The config file can only be altered later by manually opening and editting the .ini file
Save configuration to config file? (y/n): ''')
                if cfg_save_op in ('y') or cfg_save_op in ('1'):
                    topo_save = topo_path
                    if use_topo == 1:
                        topo_prompt = input('''
Topography Map Prompt
Prompts the User to input a topography map and related data every time this config is used.
Otherwise, the same topography will be used in the future (unless it can't be found)
Prompt for Topography Map? (y/n): ''')
                        if topo_prompt in ('y') or topo_prompt in ('1'):
                            topo_save = 'prompt'
                    cfg_savename = str(input('Config name: '))
                    config = configparser.ConfigParser()
                    config['Zone Types'] = {'Land Climate Zones': str(land_type),
                                            'Sea Climate Zones': str(sea_type),
                                            'Land/Sea Blend': str(blend),
                                            'Color List': str(color_type),
                                            'Custom Color List': str(col_list_path)}
                    config['Advanced Functions'] = {'Bin Size': str(bin_num),
                                                    'Interpolation Factor': str(interp),
                                                    'Dummy Sea Ice': str(dum_ice),
                                                    'Adjust Temp by Topography': str(use_topo),
                                                    'Topography Map': str(topo_save),
                                                    'Max Elevation': str(maxel),
                                                    'Min Elevation': str(minel),
                                                    'Surface Gravity': str(gravity),
                                                    'Blend by Topography Map': str(blend_topo),
                                                    'Sea Level': str(sealev)}
                    config['Climate Zone Definitions'] = {'Summer Definition': str(sum_def),
                                                        'Temperate/Continental Boundary': str(temp_def),
                                                        'Hot/Cold Arid Boundary': str(arid_def),
                                                        'Mediterranean Definition': str(med_def),
                                                        'Wet-Summer Definition': str(wet_def),
                                                        'Additional Wet Season Requirements': str(add_def),
                                                        'Med/Wet-Summer Priority': str(prio_def),
                                                        'Arid/Polar Priority': str(ice_def),
                                                        'Sea Ice Definition': str(sea_def)}
                    cfg_path = path + cfg_savename + '.ini'
                    with open(cfg_path, 'w') as cfg_file:
                        config.write(cfg_file)
                    print('Config saved to ' + str(cfg_path))
            output_name = path+input('''
Output map name: ''')
            print('''Setup Complete
''')

        MakeMap(in_files, output_name, cfg_loadname, land_type, sea_type, blend, color_type, col_list_path,
                bin_num, interp, dum_ice, use_topo, topo_path, maxel, minel, gravity, blend_topo, sealev, sum_def, temp_def,
                arid_def, med_def, wet_def, add_def, prio_def, ice_def, sea_def)
        redo = input("""
Produce another map with same inputs? (y/n): """)
        if redo in ('y') or redo in ('1'):
            res=1
        else:
            break

