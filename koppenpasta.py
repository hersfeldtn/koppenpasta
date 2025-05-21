import netCDF4 as nc
import math
import numpy as np
from scipy.interpolate import RectSphereBivariateSpline as spl_interp
from scipy.interpolate import RegularGridInterpolator as gr_interp
from scipy.ndimage import maximum_filter as max_filter
import configparser
from PIL import Image, ImageFont, ImageDraw
import os

ver_num = "2.1.1"

### GUIDANCE

# The script has been overhauled in 2.0 to be easier to modify
# To add new climate zone definitions:
# Add any new climate zones to list in "CLIMATE ZONES" section
# Add rgb color for new zones to appropriate list in "COLOR DICTIONARIES" (if in doubt, use def_color, which will always be loaded)
# Scroll down to "CLIMATE ZONE ALGORITHMS" and find the templates
# Create new _Data, _Param, and _Alg functions following these templates
# (or modify existing functions)
# If adding new functions, create a new Clim_func entry with a new key
# Set 'land_type' in the options below to desired key to make it the default
# (or try adding it to the input options if you feel ambitious)

##Contents (rough positions as of 2.0, bit outdated now):
# Options:                      76
# Climate Zones:                177
#   Color Dictionaries:         548
#   Climate Zone Names:         1282
# Main Functions:               1618
#   Input:                      1831
#       File_search:            1882
#       Get_input:              1927
#   Data Processing:            2423
#       Bin_months:             2428
#       Interp:                 2502
#       Find_Lapse:             2582
#       Read_topo:              2620
#       Calc_PET:               2675
#       Calc_GDD:               2752
#       Estimate_evap:          2796
#       Get_nc:                 2820
#       Get_pet:                2916
#       Debug_file:             2992
#       Get_params:             3128
#       Get_clims:              3206
#   Image Processing:           3269
#       Make_colmap:            3284
#       Make_key:               3400
#       Make_chart:             3434
#       Make_image:             3504
#   Command Functions:          3540
#       Save_opts:              3544
#       Make_clim:              3575
#       Make_map:               3598
#   Climate Zone Algorithms:    3624
#       Templates:              3626
#       Koppen-Geiger:          3724
#       Trewartha:              4027
#       Holdridge:              4112
#       Thornthwaite:           4316
#       Whittaker:              4496
#       Two-param KG:           4574
#       Woodward:               4655
#       IPCC:                   4758
#       WCR:                    4795
#       Common biome functions: 4849
#       Prentice:               5021
#       Pasta:                  5205
#       KG unproxied            5566
#       Sea:                    5796
#       Pasta Sea:              5933
#       Extra:                  6072


### OPTIONS

#Dictionary of default options
# priority is direct user input > config file loaded in setup > kpasta_options.cfg > script defaults

option_def = {
    
    ## Main options, Set during configuration:

    #climate type and color options:
    'land_type': 'Koppen-Geiger',       # main land climate zone type
    'land_subtype': 'full',             # land subtype
    'land_color': 'standard',           # land colors
    'sea_type': 'sea_standard',         # sea climate zone type
    'sea_subtype': 'full',              # sea climate zone subtype
    'sea_color': 'standard',            # sea colors
    'color_file': None,                 # custom color list file

    #additional processing options:
    'file_combine': 'data',             # approach for combining multiple files
    'seasonless': False,                # average all data across time before finding zones
    'blend': True,                      # blend land and sea maps based on land mask
    'bin_months': 0,                    # reduce input months by binning sequential months together
    'bin_preserve_ext': True,           # preserve extremes of maximum and minimum temperature when binning months
    'month_length': 1,                  # length of month relative to Earth for GDD calculations

    #interpolation
    'interp_scale': 0,                  # scale of interpolation
    'interp_type': 'spline',            # type of interpolation
    'dummy_ice': True,                  # add dummy ice to coastal land cells when interpolating sea ice
    'topo_map': None,                   # topography to use for temperature adjustment by topography after interpolation
    'maxel': 1000,                      # max elevation in topography map (m)
    'minel': 0,                         # min elevation in topography map (m)
    'sealev': 0,                        # sea level in topography map (m)
    'gravity': 9.81,                    # gravity used for scaling topography map to geopotential (m/s^2)
    'blend_topo': True,                 # use topography map for land/sea blend

    #output options:
    'make_key': False,                  # attempt to make a key of climates present in map
    'make_chart': False,                # make a chart of temperature and precipitation of each point, colored by zone
    'outname': 'output',                # output name

    ## Not set in configuration:

    #general
    'lapse_threshold': 981,             # threshold of geopotential difference between adjacent cells for determining empirical lapse rate
    'const_lapse_rate': None,           # constant lapse rate to use in place of empirical lapse rate (K/km, positive for increasing temp at lower elevation)
    'efficient': False,                 # use efficient version of climate algorithms, if available
    'image_scale': None,                # scale to apply to final image, using nearest-neighbor interpolation; can be number to multiply resolution or tuple of (x,y) target resoluion
    'font_size': 20,                    # font size for map key and chart
    'debug_file': False,                # produce additional .nc file containing internally used data
    'force_alt_data': False,            # force use of Alternate_Data function instead of climate-specific _Data function
    'temp_adjust': -273.15,             # factor to add for all temperature adjustments (-273.15 to convert from K to C)
    'precip_adjust': 2592000000,        # factor to multiply by for all precipitation-related adjustments (2592000000 to convert from m/s to mm/month)
    'verbose': False,                   # print additional information to command line during run

    #options used for all PET calculations (used for holdridge, prentice, pasta, and unproxied KG)
    'pet_method': 'asce-pm',            # method to calculate PET
    'pet_backup_ps': 101.3,             # backup surface pressure to use if not given (kPa)
    'pet_backup_wind': 1.0,             # backup wind speed to use if not given (m/s)
    'pet_gascon': 287.05,               # effective gas constant, in J/kg*K; can be found in DIAG file for eps outputs
    'pet_use_vegf': True,               # use vegf data to adjust PET surface data if using asce-pm method
    
    #options used for all evapotranspiration calculations (used for prentice, pasta, and unproxied KG)
    'estimate_evap': 'sea',             # where to use PET and simple soil moisture model to estimate evaporation

    #options used for all GDD calculations (used for prentice, pasta, and unproxied KG)
    'gdd_productivity_modifier': 1,     # scale to be applied to GDD to account for average photosynthetic productivity relative to Earth
    'gdd_require_contiguous': True,     # rather than count total yearly GDD, count maximum accumulated GDD in any contiguous period that remains above baseline temperature
    'gdd_indicate_inf': True,           # count GDD as infinite if temp never falls below baseline temperature (in array, returns as GDD=1 million)
    'gdd_limit_light': True,            # limit GDD based on incident radiation
    'gdd_par_ratio': 0.5,               # ratio of photosynthetically active radiation to total incident radiation

    #options for absolute temperature extremes (used for woodward, prentice, pasta, and unproxied KG)
    'temp_tunings': 'eps',              # temperature threshold tunings for different data types
    'temp_adjust_ts': True,             # adjust maxt and mint based on difference between tas and ts
    
    #sea options (for standard and pasta)
    'sea_use_ts': True,                 # use surface temperature rather than 2-meter air temperature
    'sea_ice_use_temp': False,          # use temperature <-2C to determine sea ice rather than sea ice data
    
    #koppen-geiger options                  those marked (t) also apply to trewartha
    'kg_summer_zenith': False,          # (t)define summer by solar zenith angle rather than temperature
    'kg_temperate_min': 0,              # (t)min winter temp for temperate and hot arid zones
    'kg_arid_avg':  False,              # define hot and cold arid zones by average temp of 18C rather than min temp
    'kg_med_summer_precip': 30,         # (t)max precip of driest month of summer for mediterranean zones
    'kg_wet_summer_total': False,       # (t)define wet-summer zones based on >70% of total rain in summer rather than by seasonal extremes of precip (or 10/1 ratio for trewartha)
    'kg_wet_season_req': 'med_strict',  # additional requirements for wet season defintions
    'kg_wet_summer_priority': False,    # wet-summer gets priority over med where both requirements are met, rather than reverse
    'kg_arid_polar_priority': 'arid',   # (t)priority of arid and polar regions
    'kg_trewartha_arid': False,         # use Trewartha arid threshold formula
    'kg_trewartha_seasons': False,      # use Trewartha Xs/Xw definitions
    'kg_med_as': False,                 # (t)define As zones like Cs/Ds zones, rather than just based on having more than half total rain in winter

    #holdridge options
    'h_no_pet': False,                  # index zones by precipitation and biotemperature only, with no PET input
    'h_estimate_biot_avg': False,       # estimate biotemperate from average temperature using fit to Earth values

    #pasta options                          those marked (kg) also apply to unproxied koppen-geiger and those marked (p) also aply to Prentice
    'pas_boil_pres': True,              # calculate boiling point from surface pressure, rather than using constant 100 C
    'pas_ice_def': 'ice',               # (kg)parameter for defining ice cover, for land and sea
    'pas_med_thresh': 0,                # (kg)forces alternate GrS threshold for Mediterranean zones, leave at 0 to use defaults
    'pas_gint_thresh': 1250,            # (kg)threshold of GInt for both stopping GDD accumulation and 
    'pas_simple_input': False,          # (kg)(p)alters other options at startup to accommodate only tas and pr inputs
    }


### CLIMATE ZONES

# Integer assignments for each climate zone, for ease of use in algorithms
# Can use any unique positive integer up to 65535

#Debug
debug=0

#Sea
SeaTrop=1
SeaTemp=2
SeaSeasonalIce=3
SeaPermIce=4
SeaFlat=5

#Koppen-Geiger
Af=10
Am=11
Aw=12
As=13
BWh=14
BWk=15
BSh=16
BSk=17
Csa=18
Csb=19
Csc=20
Cwa=21
Cwb=22
Cwc=23
Cfa=24
Cfb=25
Cfc=26
Dsa=27
Dsb=28
Dsc=29
Dsd=30
Dwa=31
Dwb=32
Dwc=33
Dwd=34
Dfa=35
Dfb=36
Dfc=37
Dfd=38
ET=39
EF=40

#Koppen-Geiger groups
A=50
B=51
C=52
D=53
E=54

#Extra 2-leter Koppen-Geiger
BW=60
BS=61
Cs=62
Cw=63
Cf=64
Ds=65
Dw=66
Df=67

#Koppen-Geiger reduced
TropRainforest=70
TropMonsoon=71
TropSavanna=72
HotDesert=73
HotSteppe=74
ColdDesert=75
ColdSteppe=76
Med=77
Subtropical=78
Oceanic=79
Continental=80
Subarctic=81
Tundra=82
IceCap=83

#Trewartha
TrAr = 100
TrAw = 101
TrAs = 102
TrBW = 103
TrBS = 104
TrCs = 105
TrCw = 106
TrCf = 107
TrDo = 108
TrDc = 109
TrEo = 110
TrEc = 111
TrFt = 112
TrFi = 113

#Trewartha groups
TrA = 120
TrB = 121
TrC = 122
TrD = 123
TrE = 124
TrF = 125

#Holdridge
HPolDesert=200
HDryTundra=201
HMoistTundra=202
HWetTundra=203
HRainTundra=204
HBorDesert=205
HDryScrub=206
HBorMoistForest=207
HBorWetForest=208
HBorRainForest=209
HCoolDesert=210
HCoolDesertScrub=211
HSteppe=212
HCoolMoistForest=213
HCoolWetForest=214
HCoolRainForest=215
HWarmDesert=216
HWarmDesertScrub=217
HThornSteppe=218
HWarmDryForest=219
HWarmMoistForest=220
HWarmWetForest=221
HWarmRainForest=222
HTropDesert=223
HTropDesertScrub=224
HThornWood=225
HVDryForest=226
HTropDryForest=227
HTropMoistForest=228
HTropWetForest=229
HTropRainForest=230

#Thornthwaite
THTorSat = 300
THTorWet = 301
THTorMoist = 302
THTorDry = 303
THTorSemi = 304
THTorArid = 305
THHotSat = 306
THHotWet = 307
THHotMoist = 308
THHotDry = 309
THHotSemi = 310
THHotArid = 311
THWarmSat = 312
THWarmWet = 313
THWarmMoist = 314
THWarmDry = 315
THWarmSemi = 316
THWarmArid = 317
THCoolSat = 318
THCoolWet = 319
THCoolMoist = 320
THCoolDry = 321
THCoolSemi = 322
THCoolArid = 323
THColdSat = 324
THColdWet = 325
THColdMoist = 326
THColdDry = 327
THColdSemi = 328
THColdArid = 329
THFrigSat = 330
THFrigWet = 331
THFrigMoist = 332
THFrigDry = 333
THFrigSemi = 334
THFrigArid = 335

#Thornthwaite climate variability
THPrecipLow = 340
THPrecipMed = 341
THPrecipHigh = 342
THPrecipExt = 343
THComboLow = 344
THComboMed = 345
THComboHigh = 346
THComboExt = 347
THTempLow = 348
THTempMed = 349
THTempHigh = 350
THTempExt = 351

#Prentice
PBTropRain=400
PBTropSeason=401
PBTropDry=402
PBWarmMixed=403
PBTempDecid=404
PBCoolMixed=405
PBCoolCon=406
PBTaiga=407
PBColdMixed=408
PBColdDecid=409
PBXero=410
PBWarmGrass=411
PBCoolGrass=412
PBTundra=413
PBHotDesert=414
PBSemiDesert=415
PBIce=416

#Whittaker
WTropRain = 430
WTropSeas = 431
WSubtropDes = 432
WTempRain = 433
WTempSeas = 434
WWood = 435
WTempDes = 436
WBor = 437
WTun = 438

#Woodward
WoEverTrop = 440
WoEverChill = 441
WoEverFrost = 442
WoDecidDry = 443
WoDecidWin = 444
WoConifer = 445
WoTundra = 446
WoDry = 447

#IPCC Climate Zones:
ITropWet = 460
ITropMoist = 461
ITropDry = 462
IWarmMoist = 463
IWarmDry = 464
ICoolMoist = 465
ICoolDry = 466
IBorMoist = 467
IBorDry = 468
IPolMoist = 469
IPolDry = 470

#World Climate Regions
WCRTropMoist = 480
WCRTropDry = 481
WCRTropDes = 482
WCRSubTropMoist = 483
WCRSubTropDry = 484
WCRSubTropDes = 485
WCRWarmTempMoist = 486
WCRWarmTempDry = 487
WCRWarmTempDes = 488
WCRCoolTempMoist = 489
WCRCoolTempDry = 490
WCRCoolTempDes = 491
WCRBorMoist = 492
WCRBorDry = 493
WCRBorDes = 494
WCRPolMoist = 495
WCRPolDry = 496
WCRPolDes = 497

#Pasta
TUr=500
TUrp=501
TUf=502
TUfp=503
TUs=504
TUsp=505
TUA=506
TUAp=507
TQf=508
TQfp=509
TQs=510
TQsp=511
TQA=512
TQAp=513
TF=514
TG=515
CTf=516
CTfp=517
CTs=518
CTsp=519
CDa=520
CDap=521
CDb=522
CDbp=523
CEa=524
CEap=525
CEb=526
CEbp=527
CEc=528
CEcp=529
CMa=530
CMb=531
CAMa=532
CAMb=533
CAa=534
CAap=535
CAb=536
CAbp=537
CFa=538
CFb=539
CG=540
CI=541
HTf=542
HTfp=543
HTs=544
HTsp=545
HDa=546
HDap=547
HDb=548
HDbp=549
HDc=550
HDcp=551
HMa=552
HMb=553
HMc=554
HAMa=555
HAMb=556
HAMc=557
HAa=558
HAap=559
HAb=560
HAbp=561
HAc=562
HAcp=563
HFa=564
HFb=565
HFc=566
HG=567
ETf=568
ETfp=569
ETs=570
ETsp=571
EDa=572
EDap=573
EDb=574
EDbp=575
EMa=576
EMb=577
EAMa=578
EAMb=579
EAa=580
EAap=581
EAb=582
EAbp=583
EFa=584
EFb=585
EG=586
Ada=587
Aha=588
Adc=589
Ahc=590
Adh=591
Ahh=592
Ade=593
Ahe=594

#Pasta ocean
Ofi=595
Ofd=596
Ofg=597
Og=598
Oc=599
Ot=600
Oh=601
Or=602
Oe=603

### COLOR DICTIONARIES

#Dictionaries of zone - color pairings:

#Default dictionary that's always loaded; add here to ensure color is available
def_color = {
    debug: [0,0,0]
    }

#Other dictionaries are loaded depending on climate zone type and color choice (not dependent on subtype)
sea_standard_color = {
    SeaTrop: [9,120,171],
    SeaTemp: [113,171,216],
    SeaSeasonalIce: [185,227,255],
    SeaPermIce: [226,248,255],
    SeaFlat: [113,171,216]
    }

sea_true_color = {
    SeaTrop: [20,30,66],
    SeaTemp: [20,30,66],
    SeaSeasonalIce: [20,30,66],
    SeaPermIce: [190,208,226],
    SeaFlat: [20,30,66]
    }

sea_white_color = {
    SeaTrop: [255,255,255],
    SeaTemp: [255,255,255],
    SeaSeasonalIce: [255,255,255],
    SeaPermIce: [255,255,255],
    SeaFlat: [255,255,255]
    }

koppen_standard_color = {
    Af: [0,0,254],
    Am: [0,119,255],
    Aw: [70,169,250],
    As: [127,201,255],
    BWh: [254,0,0],
    BWk: [254,150,149],
    BSh: [245,163,1],
    BSk: [255,219,99],
    Csa: [255,255,0],
    Csb: [198,199,0],
    Csc: [150,150,0],
    Cwa: [150,255,150],
    Cwb: [99,199,100],
    Cwc: [50,150,51],
    Cfa: [198,255,78],
    Cfb: [102,255,51],
    Cfc: [51,199,1],
    Dsa: [255,0,254],
    Dsb: [198,0,199],
    Dsc: [150,50,149],
    Dsd: [150,100,149],
    Dwa: [171,177,255],
    Dwb: [90,119,219],
    Dwc: [76,81,181],
    Dwd: [50,0,135],
    Dfa: [0,255,255],
    Dfb: [56,199,255],
    Dfc: [0,126,125],
    Dfd: [0,69,94],
    ET: [178,178,178],
    EF: [104,104,104]
    }

koppen_alt_red_color = {
    Af: [147,1,1],
    Am: [254,0,0],
    Aw: [255,207,207],
    As: [254,154,154],
    BWh: [255,207,0],
    BWk: [255,254,101],
    BSh: [207,143,20],
    BSk: [206,170,84],
    Csa: [0,254,0],
    Csb: [149,255,0],
    Csc: [203,255,0],
    Cwa: [180,101,0],
    Cwb: [150,102,4],
    Cwc: [94,64,0],
    Cfa: [0,48,0],
    Cfb: [1,80,1],
    Cfc: [0,120,0],
    Dsa: [254,108,253],
    Dsb: [254,182,255],
    Dsc: [230,202,253],
    Dsd: [202,204,203],
    Dwa: [204,182,255],
    Dwb: [153,124,178],
    Dwc: [138,89,178],
    Dwd: [109,36,179],
    Dfa: [48,0,48],
    Dfb: [101,1,100],
    Dfc: [203,0,203],
    Dfd: [199,21,135],
    ET: [101,255,255],
    EF: [99,150,255]
    }

koppen_true_color = {
    Af: [42,65,15],
    Am: [53,74,19],
    Aw: [73,87,32],
    As: [73,87,32],
    BWh: [213,183,133],
    BWk: [178,153,112],
    BSh: [123,112,66],
    BSk: [128,117,74],
    Csa: [112,104,58],
    Csb: [66,75,31],
    Csc: [66,75,31],
    Cwa: [78,88,36],
    Cwb: [79,81,38],
    Cwc: [135,114,68],
    Cfa: [65,80,27],
    Cfb: [62,77,27],
    Cfc: [68,78,48],
    Dsa: [148,131,85],
    Dsb: [96,92,50],
    Dsc: [66,70,31],
    Dsd: [57,66,23],
    Dwa: [74,89,34],
    Dwb: [67,83,32],
    Dwc: [57,72,23],
    Dwd: [65,71,28],
    Dfa: [66,85,29],
    Dfb: [55,75,21],
    Dfc: [52,64,20],
    Dfd: [62,71,25],
    ET: [122,119,92],
    EF: [204,217,232]
    }

#copy from main zones
koppen_group_color = {
    A: Af,
    B: BWh,
    C: Cfb,
    D: Dfa,
    E: EF
    }

koppen_twoletter_color = {
    BW: BWh,
    BS: BSh,
    Cs: Csa,
    Cw: Cwa,
    Cf: Cfb,
    Ds: Dsa,
    Dw: Dwa,
    Df: Dfa
    }

koppen_reduced_color = {
    TropRainforest: Af,
    TropMonsoon: Am,
    TropSavanna: Aw,
    HotDesert: BWh,
    ColdDesert: BWk,
    HotSteppe: BSh,
    ColdSteppe: BSk,
    Med: Csa,
    Subtropical: Cwa,
    Oceanic: Cfb,
    Continental: Dfa,
    Subarctic: Dfc,
    Tundra: ET,
    IceCap: EF
    }

trewartha_standard_color = {
    TrAr: [132,7,11],
    TrAw: [205,28,10],
    TrAs: [206,78,35],
    TrBW: [255,222,74],
    TrBS: [240,145,55],
    TrCs: [159,195,1],
    TrCw: [20,102,57],
    TrCf: [0,151,54],
    TrDo: [0,173,215],
    TrDc: [178,85,156],
    TrEo: [20,215,158],
    TrEc: [20,81,161],
    TrFt: [192,192,192],
    TrFi: [140,140,140]
    }

trewartha_kalike_color = {
    TrAr: [0,0,254],
    TrAw: [70,169,250],
    TrAs: [127,201,255],
    TrBW: [254,0,0],
    TrBS: [245,163,1],
    TrCs: [255,255,0],
    TrCw: [150,255,150],
    TrCf: [198,255,78],
    TrDo: [102,255,51],
    TrDc: [56,199,255],
    TrEo: [51,199,1],
    TrEc: [0,126,125],
    TrFt: [178,178,178],
    TrFi: [104,104,104]
    }

trewartha_group_color = {
    TrA: TrAr,
    TrB: TrBW,
    TrC: TrCf,
    TrD: TrDo,
    TrE: TrEc,
    TrF: TrFi
    }

holdridge_standard_color = {
    HPolDesert: [255,255,255],
    HDryTundra: [128,128,128],
    HMoistTundra: [96,128,128],
    HWetTundra: [64,128,144],
    HRainTundra: [32,128,192],
    HBorDesert: [160,160,128],
    HDryScrub: [128,160,128],
    HBorMoistForest: [96,160,128],
    HBorWetForest: [64,160,144],
    HBorRainForest: [32,160,192],
    HCoolDesert: [192,192,128],
    HCoolDesertScrub: [160,192,128],
    HSteppe: [128,192,128],
    HCoolMoistForest: [96,192,128],
    HCoolWetForest: [64,192,144],
    HCoolRainForest: [32,192,192],
    HWarmDesert: [224,224,128],
    HWarmDesertScrub: [192,224,128],
    HThornSteppe: [160,224,128],
    HWarmDryForest: [128,224,128],
    HWarmMoistForest: [96,224,128],
    HWarmWetForest: [64,224,144],
    HWarmRainForest: [32,224,192],
    HTropDesert: [255,255,128],
    HTropDesertScrub: [224,255,128],
    HThornWood: [192,255,128],
    HVDryForest: [160,255,128],
    HTropDryForest: [128,255,128],
    HTropMoistForest: [96,255,128],
    HTropWetForest: [64,255,144],
    HTropRainForest: [32,255,160]
    }

holdridge_vibrant_color = {
    HPolDesert: [20,82,161],
    HDryTundra: [182,254,227],
    HMoistTundra: [137,218,236],
    HWetTundra: [68,195,225],
    HRainTundra: [179,195,253],
    HBorDesert: [204,235,215],
    HDryScrub: [153,214,175],
    HBorMoistForest: [102,193,135],
    HBorWetForest: [51,173,96],
    HBorRainForest: [0,152,55],
    HCoolDesert: [236,243,204],
    HCoolDesertScrub: [220,233,163],
    HSteppe: [205,223,121],
    HCoolMoistForest: [189,215,81],
    HCoolWetForest: [174,205,42],
    HCoolRainForest: [158,194,0],
    HWarmDesert: [255,248,219],
    HWarmDesertScrub: [255,245,196],
    HThornSteppe: [255,239,172],
    HWarmDryForest: [254,235,144],
    HWarmMoistForest: [255,231,123],
    HWarmWetForest: [255,226,97],
    HWarmRainForest: [255,223,74],
    HTropDesert: [246,210,206],
    HTropDesertScrub: [238,180,173],
    HThornWood: [233,150,139],
    HVDryForest: [227,119,108],
    HTropDryForest: [219,90,76],
    HTropMoistForest: [214,59,43],
    HTropWetForest: [207,29,1],
    HTropRainForest: [186,26,0]
    }

thornthwaite_standard_color = {
    THTorSat: [67,35,108],
    THTorWet: [37,49,107],
    THTorMoist: [43,113,50],
    THTorDry: [102,112,36],
    THTorSemi: [107,69,25],
    THTorArid: [108,18,22],
    THHotSat: [116,34,144],
    THHotWet: [39,75,164],
    THHotMoist: [59,172,59],
    THHotDry: [157,169,33],
    THHotSemi: [159,104,28],
    THHotArid: [163,26,32],
    THWarmSat: [128,55,158],
    THWarmWet: [62,89,176],
    THWarmMoist: [95,193,48],
    THWarmDry: [216,228,10],
    THWarmSemi: [225,142,19],
    THWarmArid: [228,27,35],
    THCoolSat: [142,57,160],
    THCoolWet: [69,102,184],
    THCoolMoist: [110,197,43],
    THCoolDry: [237,235,11],
    THCoolSemi: [245,158,12],
    THCoolArid: [234,28,34],
    THColdSat: [180,104,186],
    THColdWet: [133,173,223],
    THColdMoist: [163,216,79],
    THColdDry: [244,242,87],
    THColdSemi: [249,205,100],
    THColdArid: [241,95,119],
    THFrigSat: [220,179,221],
    THFrigWet: [197,208,237],
    THFrigMoist: [211,236,170],
    THFrigDry: [250,248,175],
    THFrigSemi: [250,232,175],
    THFrigArid: [246,174,187]
    }

thornthwaite_variability_color = {
    THPrecipLow: [221,224,242],
    THPrecipMed: [161,173,221],
    THPrecipHigh: [69,90,173],
    THPrecipExt: [48,66,145],
    THComboLow: [220,240,207],
    THComboMed: [161,214,129],
    THComboHigh: [100,192,47],
    THComboExt: [49,131,54],
    THTempLow: [251,226,229],
    THTempMed: [243,126,148],
    THTempHigh: [234,38,44],
    THTempExt: [167,28,32]
    }

whittaker_standard_color = {
    WTropRain: [53,122,44],
    WTropSeas: [178,156,30],
    WSubtropDes: [230,190,106],
    WTempRain: [120,170,113],
    WTempSeas: [157,183,126],
    WWood: [225,125,77],
    WTempDes: [255,216,146],
    WBor: [167,201,162],
    WTun: [190,227,227]
    }

woodward_standard_color = {
    WoEverTrop: [0,38,255],
    WoEverChill: [127,201,255],
    WoEverFrost: [127,255,255],
    WoDecidDry: [182,255,0],
    WoDecidWin: [0,255,33],
    WoConifer: [255,127,237],
    WoTundra: [160,160,160],
    WoDry: [226,212,147]
    }

ipcc_standard_color = {
    ITropWet: [70,140,114],
    ITropMoist: [131, 210, 93],
    ITropDry: [246, 245, 119],
    IWarmMoist: [108,228,255],
    IWarmDry: [255,213,129],
    ICoolMoist: [201,244,137],
    ICoolDry: [190,161,209],
    IBorMoist: [158,170,218],
    IBorDry: [217,215,156],
    IPolMoist: [219,252,231],
    IPolDry: [225,225,225]
    }

wcr_standard_color = {
    WCRTropMoist: [89,59,49],
    WCRTropDry: [173,124,92],
    WCRTropDes: [208,165,133],
    WCRSubTropMoist: [240,254,175],
    WCRSubTropDry: [246,242,119],
    WCRSubTropDes: [250,250,180],
    WCRWarmTempMoist: [113,172,88],
    WCRWarmTempDry: [169,211,145],
    WCRWarmTempDes: [200,227,182],
    WCRCoolTempMoist: [66,89,73],
    WCRCoolTempDry: [173,185,147],
    WCRCoolTempDes: [207,222,201],
    WCRBorMoist: [87,111,147],
    WCRBorDry: [158,183,205],
    WCRBorDes: [210,220,230],
    WCRPolMoist: [173,181,183],
    WCRPolDry: [192,208,208],
    WCRPolDes: [228,236,238]
    }

prentice_standard_color = {
    PBTropRain: [67,93,94],
    PBTropSeason: [101,115,90],
    PBTropDry: [183,132,97],
    PBWarmMixed: [60,74,106],
    PBTempDecid: [71,126,86],
    PBCoolMixed: [82,131,145],
    PBCoolCon: [170,207,204],
    PBTaiga: [65,95,126],
    PBColdMixed: [241,130,87],
    PBColdDecid: [77,124,176],
    PBXero: [193,90,71],
    PBWarmGrass: [246,183,86],
    PBCoolGrass: [180,195,113],
    PBTundra: [159,156,203],
    PBHotDesert: [249,252,199],
    PBSemiDesert: [228,239,235],
    PBIce: [200,223,253]
    }

pasta_standard_color = {
    TUr: [0,0,255],
    TUrp: [4,0,191],
    TUf: [41,112,255],
    TUfp: [26,80,188],
    TUs: [145,180,255],
    TUsp: [95,125,196],
    TUA: [199,216,255],
    TUAp: [136,157,206],
    TQf: [55,210,192],
    TQfp: [48,141,130],
    TQs: [117,245,230],
    TQsp: [114,197,188],
    TQA: [186,253,245],
    TQAp: [174,219,213],
    TF: [83,83,147],
    TG: [30,28,109],
    CTf: [84,218,34],
    CTfp: [54,158,16],
    CTs: [167,253,129],
    CTsp: [120,192,89],
    CDa: [14,251,93],
    CDap: [0,194,65],
    CDb: [0,219,117],
    CDbp: [5,158,66],
    CEa: [172,251,214],
    CEap: [133,214,176],
    CEb: [112,240,186],
    CEbp: [54,171,120],
    CEc: [65,251,251],
    CEcp: [4,182,185],
    CMa: [180,240,51],
    CMb: [172,209,44],
    CAMa: [251,255,0],
    CAMb: [162,172,27],
    CAa: [215,194,117],
    CAap: [161,139,54],
    CAb: [197,219,118],
    CAbp: [132,171,84],
    CFa: [172,203,210],
    CFb: [180,188,192],
    CG: [153,153,153],
    CI: [94,94,94],
    HTf: [255,102,0],
    HTfp: [147,59,1],
    HTs: [253,151,83],
    HTsp: [198,89,16],
    HDa: [255,66,66],
    HDap: [223,48,48],
    HDb: [255,0,0],
    HDbp: [199,0,0],
    HDc: [181,33,48],
    HDcp: [137,11,26],
    HMa: [253,157,30],
    HMb: [217,137,18],
    HMc: [184,107,0],
    HAMa: [255,187,0],
    HAMb: [212,160,17],
    HAMc: [177,137,27],
    HAa: [245,200,163],
    HAap: [209,161,122],
    HAb: [230,164,148],
    HAbp: [202,129,109],
    HAc: [210,121,121],
    HAcp: [178,83,83],
    HFa: [154,106,106],
    HFb: [136,89,89],
    HFc: [119,60,60],
    HG: [71,31,31],
    ETf: [128,0,255],
    ETfp: [99,0,199],
    ETs: [181,115,247],
    ETsp: [143,90,196],
    EDa: [225,0,255],
    EDap: [158,0,179],
    EDb: [249,108,218],
    EDbp: [181,79,159],
    EMa: [255,26,205],
    EMb: [209,10,166],
    EAMa: [255,0,123],
    EAMb: [178,31,102],
    EAa: [193,139,159],
    EAap: [160,106,125],
    EAb: [255,184,248],
    EAbp: [200,116,193],
    EFa: [189,148,194],
    EFb: [157,118,162],
    EG: [89,52,91],
    Ada: [232,230,162],
    Aha: [255,251,204],
    Adc: [205,221,186],
    Ahc: [236,243,230],
    Adh: [234,184,184],
    Ahh: [255,224,224],
    Ade: [227,186,227],
    Ahe: [239,210,238]
    }

pastaocean_standard_color = {
    Ofi: [226,248,255],
    Ofd: [185,227,255],
    Ofg: [137,169,190],
    Og: [79,119,150],
    Oc: [113,171,216],
    Ot: [9,120,171],
    Oh: [6,79,147],
    Or: [1,45,86],
    Oe: [81,9,170]
    }

pasta_earthlike_color = {   #wip, will implement later
    TUr: [0,0,255],
    TUrp: [4,0,191],
    TUf: [41,112,255],
    TUfp: [26,80,188],
    TUs: [115,170,242],
    TUsp: [97,137,189],
    TUA: [168,194,255],
    TUAp: [140,162,212],
    TQf: [0,0,0],
    TQfp: [0,0,0],
    TQs: [0,0,0],
    TQsp: [0,0,0],
    TQA: [0,0,0],
    TQAp: [0,0,0],
    TF: [0,0,0],
    TG: [0,0,0],
    CTf: [128,0,255],
    CTfp: [99,0,199],
    CTs: [181,115,247],
    CTsp: [126,79,173],
    CDa: [14,251,93],
    CDap: [0,194,65],
    CDb: [112,240,186],
    CDbp: [54,171,120],
    CEa: [55,210,192],
    CEap: [48,141,130],
    CEb: [65,251,251],
    CEbp: [4,182,185],
    CEc: [225,0,255],
    CEcp: [158,0,179],
    CMa: [253,157,30],
    CMb: [177,137,27],
    CAMa: [251,255,0],
    CAMb: [162,172,27],
    CAa: [255,0,123],
    CAap: [178,31,102],
    CAb: [255,66,66],
    CAbp: [199,0,0],
    CFa: [186,253,245],
    CFb: [180,188,192],
    CG: [153,153,153],
    CI: [94,94,94],
    HTf: [0,0,0],
    HTfp: [0,0,0],
    HTs: [0,0,0],
    HTsp: [0,0,0],
    HDa: [0,0,0],
    HDap: [0,0,0],
    HDb: [0,0,0],
    HDbp: [0,0,0],
    HDc: [0,0,0],
    HDcp: [0,0,0],
    HMa: [0,0,0],
    HMb: [0,0,0],
    HMc: [0,0,0],
    HAMa: [0,0,0],
    HAMb: [0,0,0],
    HAMc: [0,0,0],
    HAa: [0,0,0],
    HAap: [0,0,0],
    HAb: [0,0,0],
    HAbp: [0,0,0],
    HAc: [0,0,0],
    HAcp: [0,0,0],
    HFa: [0,0,0],
    HFb: [0,0,0],
    HFc: [0,0,0],
    HG: [0,0,0],
    ETf: [0,0,0],
    ETfp: [0,0,0],
    ETs: [0,0,0],
    ETsp: [0,0,0],
    EDa: [0,0,0],
    EDap: [0,0,0],
    EDb: [0,0,0],
    EDbp: [0,0,0],
    EMa: [0,0,0],
    EMb: [0,0,0],
    EAMa: [0,0,0],
    EAMb: [0,0,0],
    EAa: [0,0,0],
    EAap: [0,0,0],
    EAb: [0,0,0],
    EAbp: [0,0,0],
    EFa: [0,0,0],
    EFb: [0,0,0],
    EG: [0,0,0],
    Ada: [232,230,162],
    Aha: [255,251,204],
    Adc: [161,184,132],
    Ahc: [214,238,191],
    Adh: [0,0,0],
    Ahh: [0,0,0],
    Ade: [0,0,0],
    Ahe: [0,0,0]
}

pasta_true_color = {
    TUr: [41,63,13],
    TUrp: [42,65,16],
    TUf: [55,74,20],
    TUfp: [59,80,24],
    TUs: [75,85,33],
    TUsp: [89,102,47],
    TUA: [107,105,53],
    TUAp: [124,116,63],
    TQf: [59,78,23],
    TQfp: [54,73,24],
    TQs: [75,80,35],
    TQsp: [67,76,30],
    TQA: [107,105,53],
    TQAp: [124,116,63],
    TF: [78,84,66],
    TG: [98,91,59],
    CTf: [59,78,23],
    CTfp: [54,73,24],
    CTs: [75,80,35],
    CTsp: [67,76,30],
    CDa: [60,78,23],
    CDap: [36,54,15],
    CDb: [55,75,21],
    CDbp: [38,62,11],
    CEa: [60,63,29],
    CEap: [38,52,18],
    CEb: [49,61,18],
    CEbp: [52,64,25],
    CEc: [62,71,24],
    CEcp: [64,74,27],
    CMa: [60,73,26],
    CMb: [51,63,22],
    CAMa: [103,97,54],
    CAMb: [118,108,68],
    CAa: [105,98,58],
    CAap: [58,68,25],
    CAb: [102,100,55],
    CAbp: [94,87,55],
    CFa: [78,84,66],
    CFb: [93,88,54],
    CG: [98,91,59],
    CI: [240,240,240],
    HTf: [55,74,20],
    HTfp: [59,80,24],
    HTs: [75,85,33],
    HTsp: [89,102,47],
    HDa: [60,78,23],
    HDap: [36,54,15],
    HDb: [55,75,21],
    HDbp: [38,62,11],
    HDc: [62,71,24],
    HDcp: [64,74,27],
    HMa: [60,73,26],
    HMb: [51,63,22],
    HMc: [51,63,22],
    HAMa: [103,97,54],
    HAMb: [118,108,68],
    HAMc: [118,108,68],
    HAa: [107,105,53],
    HAap: [124,116,63],
    HAb: [107,105,53],
    HAbp: [124,116,63],
    HAc: [107,105,53],
    HAcp: [124,116,63],
    HFa: [78,84,66],
    HFb: [93,88,54],
    HFc: [93,88,54],
    HG: [98,91,59],
    ETf: [59,78,23],
    ETfp: [54,73,24],
    ETs: [75,80,35],
    ETsp: [67,76,30],
    EDa: [60,78,23],
    EDap: [36,54,15],
    EDb: [55,75,21],
    EDbp: [38,62,11],
    EMa: [60,73,26],
    EMb: [51,63,22],
    EAMa: [103,97,54],
    EAMb: [118,108,68],
    EAa: [105,98,58],
    EAap: [58,68,25],
    EAb: [102,100,55],
    EAbp: [94,87,55],
    EFa: [78,84,66],
    EFb: [93,88,54],
    EG: [98,91,59],
    Ada: [167,137,95],
    Aha: [238,210,156],
    Adc: [177,153,110],
    Ahc: [208,181,141],
    Adh: [167,137,95],
    Ahh: [238,210,156],
    Ade: [177,153,110],
    Ahe: [208,181,141]
}

pastaocean_true_color = {
    Ofi: [240,240,240],
    Ofd: [10,10,51],
    Ofg: [10,10,51],
    Og: [10,10,51],
    Oc: [10,10,51],
    Ot: [10,10,51],
    Oh: [10,10,51],
    Or: [10,10,51],
    Oe: [10,10,51]
    }

### CLIMATE ZONE NAMES

#Zone names shown when making map key (defaults to variable name if not in list)
name_key = {
    SeaTrop: 'Tropical Sea',
    SeaTemp: 'Temperate Sea',
    SeaSeasonalIce: 'Seasonal Sea Ice',
    SeaPermIce: 'Permanent Sea Ice',
    SeaFlat: 'Sea',
    Af:  'Af: Tropical Rainforest',
    Am:  'Am: Tropical Monsoon',
    Aw:  'Aw: Tropical Savanna',
    As:  'As: Dry-summer Tropical Savanna',
    BWh: 'BWh: Hot Desert',
    BWk: 'BWk: Cold Desert',
    BSh: 'BSh: Hot Steppe',
    BSk: 'BSk: Cold Steppe',
    Csa: 'Csa: Hot-summer Mediterranean',
    Csb: 'Csb: Warm-summer Mediterranean',
    Csc: 'Csc: Cold-summer Mediterranean',
    Cwa: 'Cwa: Monsoon Humid Subtropical',
    Cwb: 'Cwb: Monsoon Temperate Oceanic',
    Cwc: 'Cwc: Monsoon Subpolar Oceanic',
    Cfa: 'Cfa: Humid Subtropical',
    Cfb: 'Cfb: Temperate Oceanic',
    Cfc: 'Cfc: Subpolar Oceanic',
    Dsa: 'Dsa: Hot-summer Mediterranean Continental',
    Dsb: 'Dsb: Warm-summer Mediterranean Continental',
    Dsc: 'Dsc: Subarctic Mediterranean Continental',
    Dsd: 'Dsd: Frigid Mediterranean Continental',
    Dwa: 'Dwa: Hot-summer Monsoon Continental',
    Dwb: 'Dwb: Warm-summer Monsoon Continental',
    Dwc: 'Dwc: Subarctic Monsoon Continental',
    Dwd: 'Dwd: Frigid Monsoon Continental',
    Dfa: 'Dfa: Hot-summer Continental',
    Dfb: 'Dfb: Warm-summer Continental',
    Dfc: 'Dfc: Subarctic Continental',
    Dfd: 'Dfd: Frigid Continental',
    ET:  'ET: Tundra',
    EF:  'EF: Ice Cap',
    A:   'A: Tropical',
    B:   'B: Arid',
    C:   'C: Temperate',
    D:   'D: Continental',
    E:   'E: Polar',
    BW:  'BW: Desert',
    BS:  'BS: Steppe',
    Cs:  'Cs: Mediterranean Temperate',
    Cw:  'Cw: Monsoon Temperate',
    Cf:  'Cf: Humid Temperate',
    Ds:  'Ds: Mediterranean Continental',
    Dw:  'Dw: Monsoon Continental',
    Df:  'Df: Humid Continental',
    TropRainforest: 'Tropical Rainforest',
    TropMonsoon: 'Tropical Monsoon',
    TropSavanna: 'Tropical Savanna',
    HotDesert: 'Hot Desert',
    ColdDesert: 'Cold Desert',
    HotSteppe: 'Hot Steppe',
    ColdSteppe: 'Cold Steppe',
    Med: 'Mediterranean',
    Subtropical: 'Humid Subtropical',
    Oceanic: 'Oceanic',
    Continental: 'Humid Continental',
    Subarctic: 'Subarctic',
    Tundra: 'Tundra',
    IceCap: 'Ice Cap',
    TrAr: 'Ar: Tropical Rainforest',
    TrAw: 'Aw: Tropical Savanna',
    TrAs: 'As: Dry-summer Tropical Savanna',
    TrBW: 'BW: Desert',
    TrBS: 'BS: Steppe',
    TrCs: 'Cs: Mediterranean Subtropical',
    TrCw: 'Cw: Monsoon Subtropical',
    TrCf: 'Cf: Humid Subtropical',
    TrDo: 'Do: Oceanic',
    TrDc: 'Dc: Continental',
    TrEo: 'Eo: Maritime Subarctic',
    TrEc: 'Ec: Continental Subarctic',
    TrFt: 'Ft: Tundra',
    TrFi: 'Fi: Ice Cap',
    TrA:  'A: Tropical',
    TrB:  'B: Arid',
    TrC:  'C: Subtropical',
    TrD:  'D: Temperate',
    TrE:  'E: Subarctic',
    TrF:  'F: Polar',
    HPolDesert: 'Polar Desert',
    HDryTundra: 'Dry Tundra',
    HMoistTundra: 'Moist Tundra',
    HWetTundra: 'Wet Tundra',
    HRainTundra: 'Rain Tundra',
    HBorDesert: 'Boreal Desert',
    HDryScrub: 'Boreal Desert Scrub',
    HBorMoistForest: 'Boreal Moist Forest',
    HBorWetForest: 'Boreal Wet Forest',
    HBorRainForest: 'Boreal Rainforest',
    HCoolDesert: 'Cool Desert',
    HCoolDesertScrub: 'Cool Desert Scrub',
    HSteppe: 'Cool Steppe',
    HCoolMoistForest: 'Cool Moist Forest',
    HCoolWetForest: 'Cool Wet Forest',
    HCoolRainForest: 'Cool Rainforest',
    HWarmDesert: 'Warm Desert',
    HWarmDesertScrub: 'Warm Desert Scrub',
    HThornSteppe: 'Warm Thorn Steppe',
    HWarmDryForest: 'Warm Dry Forest',
    HWarmMoistForest: 'Warm Moist Forest',
    HWarmWetForest: 'Warm Wet Forest',
    HWarmRainForest: 'Warm Rainforest',
    HTropDesert: 'Tropical Desert',
    HTropDesertScrub: 'Tropical Desert Scrub',
    HThornWood: 'Tropical Thorn Wood',
    HVDryForest: 'Tropical Very Dry Forest',
    HTropDryForest: 'Tropical Dry Forest',
    HTropMoistForest: 'Tropical Moist Forest',
    HTropWetForest: 'Tropical Wet Forest',
    HTropRainForest: 'Tropical Rainforest',
    THTorSat: 'Torrid Saturated',
    THTorWet: 'Torrid Wet',
    THTorMoist: 'Torrid Moist',
    THTorDry: 'Torrid Dry',
    THTorSemi: 'Torrid Semiarid',
    THTorArid: 'Torrid Arid',
    THHotSat: 'Hot Saturated',
    THHotWet: 'Hot Wet',
    THHotMoist: 'Hot Moist',
    THHotDry: 'Hot Dry',
    THHotSemi: 'Hot Semiarid',
    THHotArid: 'Hot Arid',
    THWarmSat: 'Warm Saturated',
    THWarmWet: 'Warm Wet',
    THWarmMoist: 'Warm Moist',
    THWarmDry: 'Warm Dry',
    THWarmSemi: 'Warm Semiarid',
    THWarmArid: 'Warm Arid',
    THCoolSat: 'Cool Saturated',
    THCoolWet: 'Cool Wet',
    THCoolMoist: 'Cool Moist',
    THCoolDry: 'Cool Dry',
    THCoolSemi: 'Cool Semiarid',
    THCoolArid: 'Cool Arid',
    THColdSat: 'Cold Saturated',
    THColdWet: 'Cold Wet',
    THColdMoist: 'Cold Moist',
    THColdDry: 'Cold Dry',
    THColdSemi: 'Cold Semiarid',
    THColdArid: 'Cold Arid',
    THFrigSat: 'Frigid Saturated',
    THFrigWet: 'Frigid Wet',
    THFrigMoist: 'Frigid Moist',
    THFrigDry: 'Frigid Dry',
    THFrigSemi: 'Frigid Semiarid',
    THFrigArid: 'Frigid Arid',
    THPrecipLow: 'Precipitation Low',
    THPrecipMed: 'Precipitation Medium',
    THPrecipHigh: 'Precipitation High',
    THPrecipExt: 'Precipitation Extreme',
    THComboLow: 'Combination Low',
    THComboMed: 'Combination Medium',
    THComboHigh: 'Combination High',
    THComboExt: 'Combination Extreme',
    THTempLow: 'Temperature Low',
    THTempMed: 'Temperature Medium',
    THTempHigh: 'Temperature High',
    THTempExt: 'Temperature Extreme',
    WTropRain: 'Tropical Rainforest',
    WTropSeas: 'Tropical Seasonal Forest/Savanna',
    WSubtropDes: 'Subtropical Desert',
    WTempRain: 'Temperate Rainforest',
    WTempSeas: 'Temperate Seasonal Forest',
    WWood: 'Woodland/Shrubland',
    WTempDes: 'Temperate Grassland/Cold Desert',
    WBor: 'Boreal Forest',
    WTun: 'Tundra',
    WoEverTrop: 'Tropical Broadleaf Evergreen',
    WoEverChill: 'Chill-Tolerant Evergreen Broadleaf',
    WoEverFrost: 'Frost-Tolerant Evergreen Broadleaf',
    WoDecidDry: 'Dry-Deciduous Broadleaf',
    WoDecidWin: 'Winter-Deciduous Broadleaf',
    WoConifer: 'Coniferous',
    WoTundra: 'Tundra',
    WoDry: 'Dry Grass/Shrub/Desert',
    ITropWet: 'Tropical Wet',
    ITropMoist: 'Tropical Moist',
    ITropDry: 'Tropical Dry',
    IWarmMoist: 'Warm Temperate Moist',
    IWarmDry: 'Warm Temperate Dry',
    ICoolMoist: 'Cool Temperate Moist',
    ICoolDry: 'Cool Temperate Dry',
    IBorMoist: 'Boreal Moist',
    IBorDry: 'Boreal Dry',
    IPolMoist: 'Polar Moist',
    IPolDry: 'Polar Dry',
    WCRTropMoist: 'Tropical Moist',
    WCRTropDry: 'Tropical Dry',
    WCRTropDes: 'Tropical Desert',
    WCRSubTropMoist: 'Subtropical Moist',
    WCRSubTropDry: 'Subtropical Dry',
    WCRSubTropDes: 'Subtropical Desert',
    WCRWarmTempMoist: 'Warm Temperate Moist',
    WCRWarmTempDry: 'Warm Temperate Dry',
    WCRWarmTempDes: 'Warm Temperate Desert',
    WCRCoolTempMoist: 'Cool Temperate Moist',
    WCRCoolTempDry: 'Cool Temperate Dry',
    WCRCoolTempDes: 'Cool Temperate Desert',
    WCRBorMoist: 'Boreal Moist',
    WCRBorDry: 'Boreal Dry',
    WCRBorDes: 'Boreal Desert',
    WCRPolMoist: 'Polar Moist',
    WCRPolDry: 'Polar Dry',
    WCRPolDes: 'Polar Desert',
    PBTropRain: 'Tropical Rainforest',
    PBTropSeason: 'Tropical Seasonal Forest',
    PBTropDry: 'Tropical Dry Forest/Savanna',
    PBWarmMixed: 'Evergreen/Warm Mixed Forest',
    PBTempDecid: 'Temperate Deciduous Forest',
    PBCoolMixed: 'Cool Mixed Forest',
    PBCoolCon: 'Cool Conifer Forest',
    PBTaiga: 'Taiga',
    PBColdMixed: 'Cold Mixed Forest',
    PBColdDecid: 'Cold Deciduous Forest',
    PBXero: 'Xerophytic Woods/Shrub',
    PBWarmGrass: 'Warm Grass/Shrub',
    PBCoolGrass: 'Cool Grass/Shrub',
    PBTundra: 'Tundra',
    PBHotDesert: 'Hot Desert',
    PBSemiDesert: 'Semidesert',
    PBIce: 'Ice/Polar Desert',
    TUr:  'TUr: Tropical Rainforest',
    TUrp: 'TUrp: Tropical Hyperpluvial Rainforest',
    TUf:  'TUf: Tropical Forest',
    TUfp: 'TUfp: Tropical Monsoon Forest',
    TUs:  'TUs: Tropical Moist Savanna',
    TUsp: 'TUsp: Tropical Moist Monsoon Savanna',
    TUA:  'TUA: Tropical Dry Savanna',
    TUAp: 'TUAp: Tropical Dry Monsoon Savanna',
    TQf:  'TQf: Quasitropical Forest',
    TQfp: 'TQfp: Quasitropical Monsoon Forest',
    TQs:  'TQs: Quasitropical Moist Savanna',
    TQsp: 'TQsp: Quasitropical Moist Monsoon Savanna',
    TQA:  'TQA: Quasitropical Dry Savanna',
    TQAp: 'TQAp: Quasitropical Dry Monsoon Savanna',
    TF:   'TF: Tropical Twilight',
    TG:   'TG: Tropical Dark',
    CTf:  'CTf: Subtropical Forest',
    CTfp: 'CTfp: Subtropical Monsoon Forest',
    CTs:  'CTs: Subtropical Moist Savanna',
    CTsp: 'CTsp: Subtropical Moist Monsoon Savanna',
    CDa:  'CDa: Oceanic Temperate',
    CDap: 'CDap: Oceanic Temperate Rainforest',
    CDb:  'CDb: Continental Temperate',
    CDbp: 'CDbp: Continental Temperate Rainforest',
    CEa:  'CEa: Oceanic Boreal',
    CEap: 'CEap: Oceanic Boreal Rainforest',
    CEb:  'CEb: Continental Boreal',
    CEbp: 'CEbp: Continental Boreal Rainforest',
    CEc:  'CEc: Percontinental Boreal',
    CEcp: 'CEcp: Percontinental Boreal Rainforest',
    CMa:  'CMa: Oceanic Submediterranean',
    CMb:  'CMb: Continental Submediterranean',
    CAMa: 'CAMa: Oceanic Mediterranean',
    CAMb: 'CAMb: Continental Mediterranean',
    CAa:  'CAa: Cool Dry Savanna',
    CAap: 'CAap: Cool Dry Monsoon Savanna',
    CAb:  'CAb: Cold Steppe',
    CAbp: 'CAbp: Cold Pluvial Steppe',
    CFa:  'CFa: Oceanic Tundra',
    CFb:  'CFb: Continental Tundra',
    CG:   'CG: Cold Barren',
    CI:   'CI: Ice',
    HTf:  'HTf: Supertropical Forest',
    HTfp: 'HTfp: Supertropical Monsoon Forest',
    HTs:  'HTs: Supertropical Moist Savanna',
    HTsp: 'HTsp: Supertropical Moist Monsoon Savanna',
    HDa:  'HDa: Hot Swelter',
    HDap: 'HDa: Hot Pluvial Swelter',
    HDb:  'HDb: Torrid Swelter',
    HDbp: 'HDbp: Torrid Pluvial Swelter',
    HDc:  'HDc: Boiling Swelter',
    HDcp: 'HDcp: Boiling Pluvial Swelter',
    HMa:  'HMa: Hot Subparamediterranean',
    HMb:  'HMb: Torrid Subparamediterranean',
    HMc:  'HMc: Boiling Subparamediterranean',
    HAMa: 'HAMa: Hot Paramediterranean',
    HAMb: 'HAMb: Torrid Paramediterranean',
    HAMc: 'HAMc: Boiling Paramediterranean',
    HAa:  'HAa: Hot Dry Savanna',
    HAap: 'HAap: Hot Dry Monsoon Savanna',
    HAb:  'HAb: Torrid Steppe',
    HAbp: 'HAbp: Torrid Pluvial Steppe',
    HAc:  'HAc: Boiling Steppe',
    HAcp: 'HAcp: Boiling Pluvial Steppe',
    HFa:  'HFa: Hot Parch',
    HFb:  'HFb: Torrid Parch',
    HFc:  'HDc: Boiling Parch',
    HG:   'HG: Hot Barren',
    ETf:  'ETf: Extratropical Forest',
    ETfp: 'ETfp: Extratropical Monsoon Forest',
    ETs:  'ETs: Extratropical Moist Savanna',
    ETsp: 'ETsp: Extratropical Moist Monsoon Savanna',
    EDa:  'EDa: Superseasonal Extracontinental',
    EDap: 'EDap: Superseasonal Extracontinental Rainforest',
    EDb:  'EDb: Hyperseasonal Extracontinental',
    EDbp: 'EDbp: Hyperseasonal Extracontinental Rainforest',
    EMa:  'EMa: Superseasonal Subextramediterranean',
    EMb:  'EMb: Hyperseasonal Subextramediterranean',
    EAMa: 'EAMa: Superseasonal Extramediterranean',
    EAMb: 'EAMb: Hyperseasonal Extramediterranean',
    EAa:  'EAa: Superseasonal Dry Savanna',
    EAap: 'EAap: Superseasonal Dry Monsoon Savanna',
    EAb:  'EAb: Hyperseasonal Steppe',
    EAbp: 'EAbp: Hyperseasonal Pluvial Steppe',
    EFa:  'EFa: Superseasonal Pulse',
    EFb:  'EFb: Hyperseasonal Pulse',
    EG:   'EG: Extraseasonal Barren',
    Ada:  'Ada: Warm Semidesert',
    Aha:  'Aha: Warm Desert',
    Adc:  'Adc: Cold Semidesert',
    Ahc:  'Ahc: Cold Desert',
    Adh:  'Adh: Hot Semidesert',
    Ahh:  'Ahh: Hot Desert',
    Ade:  'Ade: Extraseasonal Semidesert',
    Ahe:  'Ahe: Extraseasonal Desert',
    Ofi:  'Ofi: Permanent Sea Ice',
    Ofd:  'Ofd: Seasonal Sea Ice',
    Ofg:  'Ofg: Dark Seasonal Sea Ice',
    Og:   'Og: Dark Ocean',
    Oc:   'Oc: Cool Ocean',
    Ot:   'Ot: Tropical Ocean',
    Oh:   'Oh: Hot Ocean',
    Or:   'Or: Torrid Ocean',
    Oe:   'Oe: Extraseasonal Ocean'
    }


### MAIN FUNCTIONS

##For reference, a quick overview of the typical data flow:
# If script run directly, routine at end of script calls Get_input() for file names
#   Get_input() prompts for user input and uses File_search() to find files in provided path name(s)
# Files are then passed to Make_map() as files[] (Make_map() can also be used directly as main function, giving it a list of files or single filename)
#   Make_map() passes files[] to Make_clim()
#       Make_clim() finds the appropriate land and sea _Data(), _Param(), and _Alg() functions listed in Clim_func{}
#       files[] and functions (as land_funcs[] and sea_funcs[]) are then passed to Get_params()
#           Get_params() extracts data arrays from files[] with nc.Dataset() as dat[]
#           then runs the appropriate land _Data() function on dat[]
#               _Data() usually extracts climate data from dat[] with Get_nc(), passing it dat[] and the id key of the desired climate data
#                   Get_nc() then finds the appropriate array in dat[] and saves it as dat_ar[]
#                   -if using sequential file combination, arrays from each file are linked in sequence with np.concatenate()
#                   -or if using data-averaging file combination, arrays from each file are extracted and then averaged together.
#                   -if using interpolation, dat_ar[] is then passed to Interp() for interpolation
#                   -if binning months, dat_ar[] is then passed to Bin_months() for binning
#               each dat_ar[] is then passed back to _Data(), which combines them into all_data{} and returns it
#           Get_params() receives this as data{}
#           dat[] and data{} are then passed to the sea _Data() function,
#               which adds to data{} from dat[] as necessary and returns it
#           -if using seasonless option, all arrays in data{} are then averaged into one month
#           the appropriate land _Param() function is then run on data{}
#               _Param then processes 3d climate data (time, lat, lon) into 2d parameters (lat,lon)
#               these are all saved to all_param{} and returned
#           Get_params() receives this as params{}
#           as above, data{} and params{} are then passed to the sea _Param() function,
#               which adds to params{} from data{} and returns it
#           -if using data-averaging file combination, all data is extracted first to data{}
#            and then one data{} is passed to _Param to make params{}
#           -if using parameter-averaging file combination, data is extracted from each file individually to make separate data{}
#            and each is passed to _Param to make separate par{}, which are then averaged into one params{}
#           -if making a debug file, data{} and params{} are passed to Debug_file()
#               Debug_file() saves all constituent arrays to a netCDF file,
#                along with any land/sea masks, latitude and longitude arrays, uploaded topography, and a list of all option settings
#           Get_params() then returns the final params{}
#       Make_clim() passes params{} and functions to Get_clims()
#           Get_clims() finds an appropriate land/sea mask
#           and then iterates through the map point-by-point
#           at each point, it creates a par{} containing the value from each param{} array at that point
#           par{} is passed to the land or sea _Alg() function according to the land/sea mask
#               _Alg determines an appropriate clim id for that point and returns it
#           Get_clims() saves this to the appropriate _clims[] array.
#           an alternative "efficient" can be used to run an _Alg() function on the whole params{} at once, but is not much implemented yet
#           -if making a chart, each _clims[] array is separately passed to Make_chart()
#               which makes a chart image and returns it
#            Get_clims() saves each as a .png image file
#           -if blending land and sea, the _clims[] arrays are combined to one clims[] array using the land/sea mask
#            otherwise, they remain as separate maps
#           in either case, clims[] arrays are saved to maps{}, which is returned
#       Make_clim() receives and returns maps{}
#   Make_map() passes maps{} to Make_image()
#       Make_image() retrieves a color map with Get_colmap()
#           Get_colmap returns colmap[], which matches each clim id to an rgb color tuple
#       then for each array in maps{},
#        Make_image iterates over each point in the map,
#        reading the clim id and saving the corresponding color from colmap[] to outmap[]
#       outmap[] is converted to outim
#       -if there is any image scaling, the image is nearest-neighbor rescaled
#       all resulting outim are saved as .png image files
#       -if making a key, maps{} and colmap[] are then passed to Make_key()
#           which constructs a single key image for all clim ids in maps{} and returns it
#        Make_image() saves this as a .png file
#       it then returns the last outim, though this is not currently used in the script
#   Make_map() returns with no output
# if run directly, the script prompts for the option to run the process again, using the same files[] but separate configuration options
# in summary, the general route to follow is
# Get_input
#   > Make_map
#           > Make_clim
#                   > Get_params
#                           > _Data
#                                   > Get_nc
#                           _Data <
#                   Get_params <
#                           > _Param
#                   Get_params <
#           Make_clim <
#                   > Get_clims
#                           > _Alg
#                   Get_Clims <
#           Make_clim <
#   Make_map <
#           > Make_image
#   

##There are a couple alternatives to the standard Get_nc()
# Get_nc_if() is given both dat[] and data{}
#   it checks if the desired data is already in data{}
#   and only retrieves it from dat[] (calling reguilar Get_nc()) if not there
#   This is useful for e.g. the sea _Data functions, as it prevents repeating work
# Get_nc_adjust() is designed for allowing for interpolation
#   it's given the keys to both temperature data and elevation data
#       it reads the topography map if necessary with Read_topo()
#       finds the temperature/elevation lapse rate with Find_lapse()
#   then compares the elevation data to the topography
#   and constructs an appropriate map of temperature adjustments to apply to all temperature data


##Configuration options can be set at a number of points:
# Can be set by user input during Get_input(),
#  either by responding to prompts or directly inputting option keys and values for final prompt
#  in either case, will override all other inputs
# Get_input() also prompts for input of a custom config file
#   which is searched with Load_options(), then overwritten by user input (if not skipping rest of configuration) but overwrites all other input
# Resulting options (user + custom config) are then passed to Make_map() as in_opts
#  a dictionary of options or name of a custom config file can also be passed to Make_map() directly
#   Make_map() then passes options to Save_opts() as in_opts
#    Make_clim() and Make_image() will also call Save_opts() if given an in_opts parameter
#       Save_opts() first uses reset_default() to reset all options to the default in the script above
#       then searches for kpasta_options.cfg in same directory as script
#       if found, searches with Load_options() and adds to list of options, overriding defaults
#       if given a name of a custom config file, will search with Load_options() and add to list
#        or if given a dictionary, will load all into options
#       in either case, overriding defaults and kpasta_options
#       In all cases, Save_opts() uses add_opt() to add options to kpasta_options{}
# any other functions can then access any saved options with opt()
# Any other function can also add an option(s) to the list (formatted as a dictionary) with add_opt()
#  which later functions can then access
#  but currently no existing functions outside of Save_opts() (and one special case for Make_clim()) do so


##A few other notable data structures may be produced at various points and stored in kpasta_common{}:

# colmap[] is constructed by Get_colmap() the first time it is called
#   if making a chart, this will be in Make_chart(), within Make_clim(), within Make_map()
#   otherwise, this will be in Make_image(), within Make_map()
#   Make_key() also can call it if run directly, but shouldn't usually do so with the current script
# Get_colmap() first checks common; if no previous colmap{}, it will call Make_colmap()
#   Make_colmap() checks the options for land and sea type and color and color file name
#   if a color file has been specified in the options, it will search for the file and read it with configparser
#       for each item, it will search for the corresponding climate id in global variables,
#        then add key and rgb tuple to colmap[]
#   without a color file, it adds colors from the above color dictionaries depending on climate type and color options
#    colors for all climate subtypes are always added, as 
#   in any case, Make_colmap() uses Add_color() to update colmap[]
#       Add_color converts input dictionaries to array entries at positions in colmap[] corresponding to the clim ids
#       it also allows for dictionaries to reference previously entered ids, copying the rgb tuple from there
#   Make_colmap returns the resulting colmap[] to Get_colmap[]
# Get_colmap() saves colmap[] to kpasta_common{}, so that it can be referenced in the future rather than rebuilt for each use.

# land/sea masks are constructed a few different ways
# Typically, masks are built by get_mask()
#   this is called by some _Data() functions that need a land/sea mask
#   but otherwise will be called by Extra_Data(), constructing one from the 'lsm' parameter in ExoPlaSim outputs
#   if not using ExoPlaSim inputs, care should be taken to find an alternate method to construct a mask
# get_mask() first checks kpasta_common{} for an existing mask
#   if there is none, it will extract one from dat[] using a specified key
#    and save this to kpasta_common{} as 'mask'
#   if interpolation is being used, it will also nearest-neighbor upscale the mask to the target resolution
#    and save this to kpasta_common{} as 'mask_big'
# For interpolation with temperature adjustment, Get_nc_adjust() will use Read_topo()
#   Read_topo() finds the topography file specified in options,
#    resizes it as appropriate,
#    reads it to build a topography map and detailed land/sea mask at the target resolution,
#    and returns both.
# If blend by topography is specified, Get_nc_adjust() saves the topography-based mask to kpasta_common{} as 'mask_topo'
# When building maps, Get_clims() usually requires a mask
#   first it checks if a mask has been specified in params[]
#   then it searches kpasta_common{} for 'mask_topo'
#   then it searches for 'mask_big'
#   then it searches for 'mask'
#    such that it should always default to the largest, most detailed mask
#   if no mask is found, it treats the map as all land

# latitude and longitude arrays are used for interpolation and constructing the debug file
# Get_params() will always try to retrieve the arrays from lat and lon parameters in dat[]
#   using coords_from_file(), which then saves them to kpasta_common{}
# if this succeeds, Interp() will use these for the input scale
#   and then construct a set of larger arrays at the output resolution with get_coords(),
#    ensuring that it shares the same central longitude as the input,
#    then saves these to kpasta_common{}
# if it fails, Interp() will construct it's own coordinate arrays at both scales



## Dictionaries

global kpasta_common, kpasta_options

# empty dictionary to hold some common data arrays used across functions, like lat, lon, mask, etc.
kpasta_common = {}

# options dictionary for use across all functions
kpasta_options = option_def

# utility functions for easy and safe access to global dictionaries
def opt(option):
    return kpasta_options[option]

def common(key):
    return kpasta_common[key]

def add_opt(option):
    kpasta_options.update(option)

def add_common(key, data):
    kpasta_common[key] = data
    verb(f'    Saved {key} to common data')

def reset_default():
    global kpasta_common
    global kpasta_options
    kpasta_options = {}
    add_opt(option_def)
    kpasta_common = {}

#Current path name
path = os.path.join(os.path.dirname(__file__), '')
        

#Dictionary of climate functions, to be filled in later
Clim_func = {}

#Function for verbose output, just to save the extra line
def verb(rep):
    if opt('verbose'):
        print(rep)
    return

## Input

#Decide how an option value string should be interpreted and returns the appropriate data type
def interp_opt(opt_val):
    opt_type = 0
    if not isinstance(opt_val,str):
        return opt_val
    for c in opt_val:
        if c.isnumeric() or c=='-':     #treat as int if all numbers or minus sign
            pass
        elif c=='.':        #if also has a decimal, treat as float
            opt_type = 1
        else:               #if has other characters, treat as string
            opt_type = 2
            break
    if opt_type == 2:   #strings have a few special cases:
        if opt_val in ("none","None","NONE"):  #interpret as none type
            opt_key = None
        elif opt_val in ("true","True","TRUE"):   #interpret as boolean true
            opt_key = True
        elif opt_val in ("false","False","FALSE"):    #interpret as boolean false
            opt_key = False
        elif "(" in opt_val and "," in opt_val:        #attempt to interpret as tuple
            try:
                opt_key = tuple([int(v) for v in opt_val[1:-1].split(',')])
            except:
                opt_key = str(opt_val)      #interpret as string if tuple failed
        else:
            opt_key = str(opt_val)
    elif opt_type == 1:
        opt_key = float(opt_val)
    else:
        opt_key = int(opt_val)
    return opt_key

#load options from config file
# cfg_file: config file to load options from
# opt_lo: lower-priority options: override these with options from file
# opt_hi: higher-priority options: override file options with these
def Load_options(cfg_file, opt_lo={}, opt_hi={}):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(cfg_file)
    in_opts = opt_lo
    for section in cfg.sections():
        for k, v in cfg.items(section):
            in_opts.update({k: interp_opt(v)})
    in_opts.update(opt_hi)
    return in_opts

#Find .nc files in path, returning list of files
# if path is a .nc file, list contains that single file
# if path is folder, finds every .nc file in folder or subdirectories
def File_search(fpath):
    files = []
    if os.path.exists(fpath):
        if os.path.isdir(fpath):
            for p, n, fs in os.walk(fpath):
                for f in fs:
                    if f.endswith(".nc"):
                        n_f = os.path.join(p, f)
                        files.append(n_f)
                        print(" Found "+str(n_f))
            if files:
                print(" Found all .nc files in "+fpath)
            else:
                print(" No .nc files found in "+str(fpath))
        else:
            if fpath.endswith(".nc"):
                print(" File found")
                files.append(fpath)
            else:
                print(f" Error: {fpath} is not .nc file")
    else:
        print(" No file found at "+str(fpath))
    return files

#Wrapper for File_search() that checks for snapshots and adds warning
def Prompt_files(fpath):
    files = File_search(fpath)
    files2 = []
    snaps = []
    for f in files:
        if 'SNAP' in f:
            snaps.append(f)
        else:
            files2.append(f)
    if len(snaps) > 0:
        print('''
WARNING: input appears to include snapshot files:''')
        for s in snaps:
            print(s)
        if not Prompt_bool('''Snapshot files record data at a specific point in time,
    not as an average across whole months,
    and so are unsuitable for climate classification purposes
Include these files? (any other files will still be included) (y/n): '''):
            files = files2
    return files

#Prompts with string and returns True if input contains y, Y, or 1
def Prompt_bool(prompt):
    ans = input(prompt)
    if ('y') in ans or ('Y') in ans or ('1') in ans:
        return True
    else:
        return False

#Prompt with string and then picks from opt_list based on response and returns result
def Pick_from_list(prompt, opt_list):
    ans = input(prompt)
    while True:
        try:
            choice = opt_list[int(ans)]
            break
        except:
            ans = input(' Invalid input, try again: ')
    return choice

#Get user input by command line
# returns input files and dictionary of input options
def Get_input(first=True, fi=[]):

    in_files = fi
    in_opts = {}
    if first:       #only show header and get file input on first run
        print(f'''
Koppenpasta climate zones interpreter and mapmaker
    version {ver_num}
Script written 2021 by Nikolai Hersfeldt
    of worldbuildingpasta.blogspot.com
For use with NetCDF output files from ExoPlaSim GCM

''')

        #ask for input file
        while not in_files:
            in_files = Prompt_files(path+input('Input NetCDF filename or folder of files: '))

        #Prompts for configuration
        if not Prompt_bool('Advanced setup? (y/n): '):
            return in_files, in_opts
        
        if Prompt_bool('''
Additional Input Files
If added, values across all years will be averaged together.
All files must be .nc files with same resolution and number of months.
Add additional inputs? (y/n): '''):
            while True:
                nextin = input('Add input file ("stop" for no more input files): ')
                if nextin == ('stop') or nextin == ('STOP'):
                    break
                else:
                    in_files.append(Prompt_files(path+nextin))

    if Prompt_bool('''
Load alternative config file? (y/n): '''):
        cfg_loadname = path+input('Config filename: ')
        while True:
            if os.path.exists(cfg_loadname):
                print(' Config file found')
                try:
                    cfg_opts = Load_options(cfg_loadname)
                    if Prompt_bool(' Options loaded; skip rest of configuration? (y/n)'):
                        return in_files, cfg_opts
                    else:
                        in_opts = cfg_opts
                        break
                except:
                    print(' Error: failed to read config file')
            else:
                print(' No file found at '+str(cfg_loadname))
            new_try = input(' Try again or input "stop" to skip and continue with configuration: ')
            if new_try == ('stop') or new_try == ('STOP'):
                break
            else:
                cfg_loadname = path+new_try
            
    in_opts['land_type'] = Pick_from_list(('''
Land Climate Zone Type
0:  Koppen-Geiger climate zones (default)
1:  Koppen-Trewartha climate zones
2:  Holdridge life zones
3:  Thornthwaite-Feddemma climate types
4:  Prentice et al. 1992 BIOME1 model
5:  Pasta bioclimate zones
6:  Whittaker biomes
7:  Woodward vegetation types
8:  IPCC climate zones
9:  World climate regions (Sayre et al. 2020)
10: Two-parameter Koppen-alike
11: Un-proxied Koppen-Geiger
Set land climate zone type: '''),
        ('Koppen-Geiger',
            'Trewartha',
            'Holdridge',
            'Thornthwaite',
            'Prentice',
            'Pasta',
            'Whittaker',
            'Woodward',
            'IPCC',
            'WCR',
            'TwoParamKG',
            'KG_unproxied'
            ))
    if in_opts['land_type'] in ('Koppen-Geiger', 'KG_unproxied', 'TwoParamKG'):
        if in_opts['land_type'] == 'TwoParamKg':
            in_opts['land_subtype'] = Pick_from_list(('''
Two-Parameter Koppen-Alike Subtype
0: 12 Koppen-alike climate zones (default)
1: 5 Koppen Groups (A,B,C,D,E) only
Set climate zone subtype: '''),
            ('full',
                'groups'
            ))
        else:
            in_opts['land_subtype'] = Pick_from_list(('''
Koppen-Geiger Subtype
0: Full Koppen set of 31 climate zones (default)
1: As above, but exclude wet-summer savanna (Aw) / dry-summer savanna (As) distinction (common in formal maps)
2: Two-letter; categorize only by first 2 letters in designation, for 14 zones
3: Alternative reduced set of 14 climate zones (common in worldbuilding circles)
4: 5 Koppen Groups (A,B,C,D,E) only
Set climate zone subtype: '''),
            ('full',
                'no_As',
                'two_letter',
                'reduced',
                'groups'
            ))
        color_prompt = ('''
Koppen-Geiger Land Colors
0: Default (blue rainforest); used by Wikipedia, Artifexian, etc.
1: Alternate (red rainforest); used in some research papers, usually with older definitions
2: "True" color; fills each zone with their average color on Earth, based on satellite imagery
3: Import custom color list (see defaultcolor.ini for template)
Set land colors: ''')
        color_opts = ('standard',
            'alt_red',
            'true',
            'file')

    elif in_opts['land_type'] == 'Trewartha':
        in_opts['land_subtype'] = Pick_from_list(('''
Trewartha Subtype
0: Typical Trewartha set of 14 zones (default)
1: As above, but exclude As and Eo
2: 6 Trewartha groups (A,B,C,D,E,F) only)
Set climate zone subtype: '''),
            ('full',
                'reduced',
                'groups'
                ))
        color_prompt = ('''
Trewartha Land Colors
0: Default; based on Belda et al. 2014
1: Koppen-alike; copies default color of most similar Koppen-Geiger zone
2: Import custom color list (see defaultcolor.ini for template)
Set land colors: ''')
        color_opts = ('standard',
                        'kalike',
                        'file')

    elif in_opts['land_type'] == 'Holdridge':
        in_opts['land_subtype'] = 'full'
        color_prompt = ('''
Holdridge Land Colors
0: Default grey-yellow-brown colors, based on Wikipedia
1: More vibrant multicolor alternative, based on Audebert et al. 2024
2: Import custom color list (see defaultcolor.ini for template)
Set land colors: ''')
        color_opts = ('standard',
                        'vibrant',
                        'file')
    elif in_opts['land_type'] == 'Thornthwaite':
        in_opts['land_subtype'] = Pick_from_list(('''
Thornthwaite Subtype
0: Thermal and moisture types (default)
1: Climate variability type
Set climate zone subtype: '''),
            ('full',
                'variability'))
        color_prompt = ('''
Thornthwaite Colors
0: Default rainbow colors, based on Feddema 2005
1: Import custom color list (see defaultcolor.ini for template)
Set Land colors: ''')
        color_opts = ('standard', 'file')
    elif in_opts['land_type'] == 'Pasta':
        in_opts['land_subtype'] = Pick_from_list(('''
Pasta Subtype
0: Full set
1: As above, but exclude pluvial (Xxp) zones
2: Earthlike-only: exclude Hot (H) zone, Extraseasonal (E) zones, and their arid counterparts,
    and ignore light limitations
3: As above, but exclude pluvial (Xxp) zones as well
Set climate zone subtype: '''),
            ('full',
            'no_pluv',
            'earthlike',
            'earthlike_no_pluv'))
        if False: #in_opts ['land_subtype'] in ('earthlike', 'earthlike_no_pluv'):  #I may implement a separate earthlike-focused color set one day
            color_prompt = ('''
Pasta Colors
0: Default colors
1: Earthlike-focused: more diversity of color in Earthlike set of zones
2: "True" color; fills each zone with their average color on Earth, based on satellite imagery
3: Import custom color list (see defaultcolor.ini for template)
Set Land colors: ''')       
            color_opts = ('standard', 'earthlike', 'true', 'file')
        else:
            color_prompt = ('''
Pasta Colors
0: Default colors
1: "True" color; fills each zone with their average color on Earth, based on satellite imagery
    (or that of nearest analogue on Earth)
2: Import custom color list (see defaultcolor.ini for template)
Set Land colors: ''') 
            color_opts = ('standard', 'true', 'file')
        
    else:
        color_prompt = ('''
Land Zone Colors
0: Default
1: Import custom color list (see defaultcolor.ini for template)
Set Land colors: ''')
        color_opts = ('standard', 'file')
                        
    in_opts['land_color'] = Pick_from_list(color_prompt, color_opts)
    if in_opts['land_color'] == 'file':
        col_list = input('Custom color list filename: ')
        while True:
            col_list_path = path + col_list
            if os.path.exists(col_list_path):
                in_opts['land_color_file'] = col_list_path
                print(' Custom color list found')
                break
            else:
                print(' No file found at '+str(col_list_path))
                while True:
                    col_list = input(('Try again or input "stop" to pick other color option: '))
                    if col_list == ('stop') or col_list == ('STOP'):
                        in_opts['land_color'] = Pick_from_list(color_prompt, color_opts)
                        if in_opts['land_color'] == ('file'):
                            print('indecisive, huh?')
                        else:
                            break
                    else:
                        break
                if col_list == ('stop') or col_list == ('STOP'):
                    break
                
    sea_choice = Pick_from_list(('''
Sea Climate Zone Type
0: Set of 4 zones of my own invention; tropical, temperate, seasonal ice, permanent ice (default)
1: As above but exclude tropical/temperate distinction, showing only ice cover
2: Pasta bioclimate sea zones
3: No sea climate zones, just flat blue
4: Exclude seas entirely, produce no output (appropriate for all-land planets)
Set sea climate type:'''),
        ('full',
        'no_trop',
        'Pasta',
        'flat',
        'none'))
    if sea_choice == 'Pasta':
        in_opts['sea_type'] = 'sea_Pasta'
        in_opts['sea_subtype'] = Pick_from_list(('''
Pasta Sea Subtype
0: Full set of 10 sea zones
1: Earthlike-only: exclude dark (Ofg,Og), hot (Oh,Or), and extraseasonal (Oe) seas
2: As above, but also exclude tropical (Ot) seas, showing ice cover only
Set sea subtype: '''),
            ('full',
            'earthlike',
            'no_trop'))
        color_prompt = ('''
Pasta Sea Colors
0: Default shaded blue colors
1: "True" color; based on satellite imagery, dark blue seas and white Ofi
2: ''')
        color_opts = ('standard',
                      'true',
                      'file')
    elif sea_choice == 'none':
        in_opts['sea_type'] = ('sea_none')
        in_opts['sea_subtype'] == ('none')
    else:
        in_opts['sea_type'] = 'sea_standard'
        in_opts['sea_subtype'] = sea_choice
        color_prompt = ('''
Sea Colors
0: Default shaded blue colors
1: "True" color; based on satellite imagery, dark blue seas and blue-tinged white for permanent ice
2: White for all zones
3: ''')
        color_opts = ('standard',
            'true',
            'white',
            'file')

    if sea_choice != 'none':
        if in_opts['land_color'] == ('file'):
            color_prompt += ('Use imported custom color list (will use same list as for land)')
        else:
            color_prompt += ('Import custom color list')
        color_prompt +='''
Set sea colors: '''
        in_opts['sea_color'] = Pick_from_list(color_prompt, color_opts)
        if in_opts['sea_color'] == ('file') and in_opts('land_color') != ('file'):
            col_list = input('Custom color list filename: ')
            while True:
                col_list_path = path + col_list
                if os.path.exists(col_list_path):
                    in_opts['land_color_file'] = col_list_path
                    print(' Custom color list found')
                    break
                else:
                    print(' No file found at '+str(col_list_path))
                    while True:
                        col_list = input(('Try again or input "stop" to pick other color option: '))
                        if col_list == ('stop') or col_list == ('STOP'):
                            in_opts['sea_color'] = Pick_from_list(color_prompt, color_opts)
                            if in_opts['sea_color'] == ('file'):
                                print('indecisive, huh?')
                            else:
                                break
                        else:
                            break
                    if col_list == ('stop') or col_list == ('STOP'):
                        break

    if len(in_files) > 1:
        in_opts['file_combine'] = Pick_from_list(('''
File Combination Method
0: Average data: average climate data (e.g. monthly temp, precipitation) between files,
    then determine parameters (e.g. hottest temp, total precip) from averages, treating as one year.
1: Average parameters: determine parameters separately for each file, then average these results together;
    slower, but better if seasonal timing may vary significantly between files.
2: Sequential: treat files as portions of one long year, linking data from each in order of input
Set combination method: '''),
            ('data',
            'param',
            'seq'))
    
    if (in_opts['land_type'] in ('Prentice','Pasta','KG_unproxied')
        or (in_opts['sea_type'] == 'sea_Pasta' and in_opts['sea_subtype'] == 'full')):
        monthl = float(input('''
Month length
Length of each month (i.e. each data point in file) relative to 30-day, 720-hour Earth month
Used for counting length of growing periods for GDD or growth interruption
Input "1" or "0" for Earth-length months.
Set month length: '''))
        if monthl == 0:
            monthl = 1
        in_opts['month_length'] = monthl
    in_opts['seasonless'] = Prompt_bool('''
Seasonless Zones
Will ignore seasonal changes and use annual averages only,
    appropriate for very short years or tidal-locked worlds without libration
Use seasonless zones? (y/n): ''')
    
    in_opts['blend'] = not Prompt_bool('''
Land/Sea Maps
By default, the land/sea mask in the NetCDF file is used to blend land and sea climates into a single map:
    tiles marked land show land climate zones, tiles marked sea show sea climate zones.
But can alternatively produce distinct land and sea climate maps
    (so you can combine them with a high-res land/sea mask yourself).
Produce separate land and sea maps? (y/n): ''')

    in_opts['bin_months'] = int(input('''
Bin Size
Bins adjacent "months" in the NetCDF file together.
    e.g. if you have a NetCDF file with 36 months, a value of "3" will average together data from every 3 months
    to produce 12 months, which will then be used for determining climate zones.
Ideally there should be a whole number of bins per year
    otherwise, will loop around to start of the year for final bin to use appropriate number of months
If using sequential file combination, binning is applied afterwards to combined year
    but in all cases is applied before determining climate parameters
For an input of "1" or "0" the script will not bin months together.
Set bin size: '''))
    if in_opts['bin_months'] > 1 and in_opts['land_type'] in ('Woodward','Prentice','Pasta','KG_unproxied'):
        in_opts['bin_preserve_ext'] = not Prompt_bool('''
Binning of Temperature Extremes
By default, max temperature in each bin is maximum of the binned months
    and min temerature is the minimum
    rather than averages, as used for most parameters.
But can choose to average these instead.
Average rather than preserve temperature extremes? (y/n): ''')
    
    in_opts['interp_scale'] = float(input('''
Interpolation Rescaling Factor
Multiplies the latitude and longitude of the output and interpolates the climate data up to that resolution
    e.g. applying a rescaling factor of 4 to a model with 32x64 resolution will get an output with 128x256 resolution.
Interpolation essentially just averages the values of adjacent cells, so isn't as accurate as actual modelling
    and in particular will not reflect the impact of topographical features smaller than those input to ExoPlaSim;
    but it does reflect the underlying data better than trying to interpolate "by eye" from the output maps.
CAUTION: High rescaling factors can significantly extend the script's runtime;
    very high may cause the script to exceed available memory and crash.
Output resolution will be rounded to nearest integers
Input "0" for no interpolation.
Set interpolation rescaling factor: '''))

    if in_opts['interp_scale'] > 0:
        in_opts['interp_type'] = Pick_from_list(('''
Interpolation Method
0: spline: scipy's SmoothSphereBivariateSpline function (default)
Others are methods for scipy's RrgularGridInterpolator,
    (with input map modified to ensure proper wrapping around globe)
    consult scipy docs for more details:
1: linear: recommended alternative if spline doesn't work well
2: nearest: essentially just upscaling the map resolution without interpolation (but slower than just doing that to final map)
3: slinear
4: cubic
5: quintic
6: pchip
Interpolation method: '''),
                ('spline',
                    'linear',
                    'nearest',
                    'slinear',
                    'cubic',
                    'quintic',
                    'pchip'
                    ))
            
        if in_opts['sea_type'] != ('sea_none'):
            in_opts['dummy_ice'] = Prompt_bool('''
"Dummy" Sea Ice Cover
Adds "dummy" sea ice to land whenever there is sea ice in at least 1 adjacent cell, before interpolation;
    this helps ensure sea ice zones remain contiguous with land when interpolated
Add "dummy" sea ice? (y/n): ''')

        if Prompt_bool('''
Adjust Temperature to Topography
Allows you to upload a higher-resolution heightmap, and then after interpolation, temperature is adjusted
    based on deviation between that heightmap and an interpolated heightmap from the model output.
Precipitation is not adjusted, so this is still not as accurate as higher-resolution modelling,
    but it is an improvement and gives the map some "texture".
Upload topography? (y/n): '''):
            topo_map = input(('''
Topography Map
CAUTION: Must be greyscale (white high, black low) and seas should be marked as flat surfaces at sea level
Does not need to match output resolution, will be bilinearly interpolated to fit
Topography map filename: '''))
            while True:
                topo_path = path + topo_map
                if os.path.exists(topo_path):
                    in_opts['topo_map'] = topo_path
                    print(' Topo map found')
                    break
                else:
                    print(' No file found at '+str(topo_path))
                    topo_map = input('Try again or input "stop" to not adjust temperature by topography: ')
                    if topo_map == ('stop') or topo_map == ('STOP'):
                        in_opts['topo_map'] = None
                        break
            if in_opts['topo_map']:
                in_opts['maxel'] = float(input('Highest Map Elevation (m): '))
                in_opts['minel'] = float(input('Lowest Map Elevation (m): '))
                in_opts['sealev'] = float(input('Sea Level (m): '))
                in_opts['gravity'] = float(input('Surface Gravity (m/s^2): '))
                if in_opts['blend']:
                    in_opts['blend_topo'] = Prompt_bool('''
Blend by Topography Map
Blends Land and Sea Maps at Interpolated Resolution using Topography Map
    with specified sea level
Blend by map? (y/n): ''')

    add_ims = Pick_from_list(('''
Additional Images
Map Key: Attempts to generate key of all climate zones present in map and saves as separate image
    note: may use internal script names for climate zones
Climate Chart: Charts every point in map by annual average temperature in (x axis, degrees C)
    and total annual precipitation in (y axis, cm/Earth year), each colored by climate zone
    note: points with similar parameters may overwrite each other; will display most common for each pixel
0: make no additional images
1: make map key only
2: make climate chart only
3: make map key and climate chart
Images choice: '''),
                    (0,1,2,3))
    if add_ims == 1 or add_ims == 3:
        in_opts['make_key'] = True
    if add_ims > 1:
        in_opts['make_chart'] = True
    
    if Prompt_bool('''
Additional Configuration Options
Can input any other configuration options, including all those found in kpasta_options.cfg,
    or any others added to the script,
    by inputting key and desired value
Input more configuration options? (y/n): '''):
        while True:
            opt_key = input('''
Option key (should be written as in kpasta_options.cfg)
Or input "stop" for not more options: ''')
            if opt_key == ('stop') or opt_key == ('STOP'):
                break
            else:
                opt_val = interp_opt(input('Option value: '))
                in_opts[opt_key] = interp_opt(opt_val)

    out = input('''
Output map name: ''')

    if out in ('0','n','N'):
        while True:
            out2 = input(f''' Confirm you want to name you output "{out}" and weren't just mindlessly skipping through options
    (input 1, Y, or y) or input proper name: ''')
            if out2 != out:
                if out2 not in ('1','Y','y'):
                    out = out2
                break

    in_opts['outname'] = path+out
    print('''
Setup Complete
''')
    return in_files, in_opts


## Data processing

#Average sequential months together into bins
#   data: data array
#   nbin: number of months to bin together
#   ext: -1 for min, 0 for average, 1 for max
def Bin_months(data, nbin, ext=0):
    if not opt('bin_preserve_ext'):
        ext = 0
    d_shape = data.shape
    n_t = math.ceil(d_shape[0]/nbin)
    n_data = np.empty((n_t,d_shape[1],d_shape[2]), dtype=data.dtype)
    if ext > 0:     #find maximum in each bin
        verb('     Binning by averaging data')
        for n in range(n_t):
            if (n+1)*nbin > d_shape[0]:
                n_data[n,:,:] = np.maximum(np.amax(data[n*nbin:,:,:], 0), np.amax(data[:n*nbin+nbin-d_shape[0],:,:], 0))    #if there aren't enough months to fill last bin, loop around to start of year
            else:
                n_data[n,:,:] = np.amax(data[n*nbin:(n+1)*nbin,:,:], 0)

    elif ext < 0:   #find minimum in each bin
        verb('     Binning by finding minimum')
        for n in range(n_t):
            if (n+1)*nbin > d_shape[0]:
                n_data[n,:,:] = np.minimum(np.amin(data[n*nbin:,:,:], 0), np.amin(data[:n*nbin+nbin-d_shape[0],:,:], 0))
            else:
                n_data[n,:,:] = np.amin(data[n*nbin:(n+1)*nbin,:,:], 0)

    else:           #find average in each bin
        verb('     Binning by finding maximum')
        for n in range(n_t):
            if (n+1)*nbin > d_shape[0]:
                n_data[n,:,:] = (np.sum(data[n*nbin:,:,:], 0) + np.sum(data[:n*nbin+nbin-d_shape[0],:,:], 0))/nbin
            else:
                n_data[n,:,:] = np.mean(data[n*nbin:(n+1)*nbin,:,:], 0)
    return n_data

#Calculate appropriate target resolution from interp scale and input resolution
def make_res(in_res, scale=None):
    if not scale:
        if opt('interp_scale'):
            scale = opt('interp_scale')
        else:
            print('  No interpolation scale provided to make_res, returning input scale')
            return in_res
    if isinstance(scale, tuple):
        scale_out = (scale[1],scale[0])  #switch from x,y to lat,lon
    else:
        scale_out = (round(in_res[0] * scale), round(in_res[1] * scale))
    verb(f'   Calculated output resolution as {scale_out}')
    return scale_out

#Check if res has already been made, and if not make it, returning res either way
def get_res(in_res, scale=None):
    try:
        res = common('res')
    except:
        res = make_res(in_res, scale)
        add_common('res', res)
    return res

# return lat and lon if they have already been made, or make them
#   big for big versions of coords
#   startlon is starting longitude for left side of map
def get_coords(shape, big=False, startlon=0):
    try:
        if big:
            lat = common('lat_big')
            lon = common('lon_big')
        else:
            lat = common('lat')
            lon = common('lon')
    except:
        lat = np.linspace(math.pi/2 - math.pi/(2*shape[0]), -math.pi/2 + math.pi/(2*shape[0]), shape[0])
        lon = np.linspace(startlon + math.pi/shape[1], startlon + 2*math.pi - math.pi/shape[1], shape[1])
        if big:
            verb('   Calculated lat and lon arrays at output resolution')
            add_common('lat_big', lat)
            add_common('lon_big', lon)
        else:
            verb('   Calculated lat and lon arrays')
            add_common('lat', lat)
            add_common('lon', lon)
    return lat, lon

#get lat and lon from file if they haven't already been made
def coords_from_file(dat, latkey, lonkey, deg=True): 
    try:
        lat = common('lat')
        lon = common('lon')
    except:
        verb('    Extracting lat and lon arrays from file')
        lat = dat[latkey][:]
        lon = dat[lonkey][:]
        if deg:
            lat *= math.pi/180
            lon *= math.pi/180
        add_common('lat', lat)
        add_common('lon', lon)
    return lat,lon

#Interpolate to higher resolution
#   data: input data array
#   res: target (y,x) resolution
#   interp_type: scipy interp_type
#   coords: (lat, lon) arrays for spline method
#   mask: land/sea mask for dummy ice
#   dummy_ice: add dummy sea ice
def Interp(data, res=None, interp_type=None, coords_in=None, coords_out=None, mask=None, dummy_ice=False):
    verb('     Interpolating data')
    squeeze_at_end = False
    if data.ndim < 3:
        data = np.expand_dims(data, 0)  #add time axis to 2d array so that later functions can assume it
        squeeze_at_end = True
    d_shape = data.shape
    if not res:
        res = get_res((d_shape[1],d_shape[2]), opt('interp_scale'))
    if not interp_type:
        interp_type = opt('interp_type')
    if coords_in:
        lat = coords_in[0]
        lon = coords_in[1]
    else:
        lat,lon = get_coords((d_shape[1], d_shape[2]))
    if coords_out:
        lat_out = coords_out[0]
        lon_out = coords_out[1]
    else:
        startlon = (lon[0] + lon[-1]) / 2 - math.pi     #find input longitude midpoint and set left edge to 1 pi less
        lat_out, lon_out = get_coords(res, big=True, startlon=startlon)
    
    lat_out,lon_out = np.meshgrid(math.pi/2 - lat_out,lon_out)
    lat_out = lat_out.ravel()
    lon_out = lon_out.ravel()
    
    if dummy_ice:
        verb('     Applying dummy sea ice for interpolation')
        go = True
        if not mask:
            try:
                mask = get_mask
            except:
                print('  No land/sea mask provided to interp function; not applying dummy sea ice')
                go = False
        if go:
            data = np.where(mask, max_filter(data, (0,3,3), mode=('constant','wrap','nearest'), cval=1.0), data)   #apply dummy ice to land areas by copying max of neighboring sea ice values
    n_data = np.empty((d_shape[0],res[0],res[1]))
    if interp_type == 'spline':
        for t in range(d_shape[0]):
            lon_i = lon - (lon[0] + lon[-1]) / 2      #shift lon to center at 0
            lon_out_i = lon_out - (lon_out[0] + lon_out[-1]) / 2
            interp = spl_interp(math.pi/2 - lat, lon_i, data[t,:,:])
            n_data[t,:,:] = interp.ev(lat_out,lon_out_i).reshape((res[1],res[0])).T
    else:
        lat = np.concatenate(([2*lat[0]-lat[1]], lat, [2*lat[-1]-lat[-2]]), 0)  #if not using spherical spline, copy each end of map to other side to ensure proper wrapping
        lon = np.concatenate(([lon[-1]-2*math.pi], lon, [lon[0]+2*math.pi]), 0)
        for t in range(d_shape[0]):
            new_row1 = np.full((1,d_shape[2]), np.mean(data[t,0,:]))
            new_row2 = np.full((1,d_shape[2]), np.mean(data[t,-1,:]))
            data_t = np.concatenate((new_row1, data[t,:,:], new_row2), 0)
            data_t = np.concatenate((data_t[:,[-1]], data_t, data_t[:,[0]]), 1)
            if data_t.dtype == np.float32:
                data_t = data_t.astype(np.float64)      #workaround for scipy bug
            interp = gr_interp((math.pi/2 - lat,lon), data_t, method = interp_type)
            n_data[t,:,:] = interp((lat_out,lon_out)).reshape((res[1],res[0])).T
    if squeeze_at_end:
        n_data = np.squeeze(n_data, 0)  #remove time axis if it was added earlier
    if np.amin(data) > 0:   # prevents spline interpolation adding unphysical sign changes
        n_data = np.maximum(n_data, 0)
    elif np.amax(data) < 0:
        n_data = np.minimum(n_data, 0)

    return n_data

#Stack all 8 neighbors for each cell in a 2d array on a new axis at 0
def stack_neighbors(ar):
    up1 = np.roll(ar,1,0)
    down1 = np.roll(ar,-1,0)
    return np.stack ((
        up1,
        down1,
        np.roll(ar,1,1),
        np.roll(ar,-1,1),
        np.roll(up1,1,1),
        np.roll(up1,-1,1),
        np.roll(down1,1,1),
        np.roll(down1,-1,1)
        ), 0)

#Find lapse rate based on relative temperature and elevation between neighboring cells
# returns array of lapse rates
#   data: temperature data
#   elev: model elevation data
def Find_lapse(data, elev):
    verb('    Finding elevation differences between cells')
    elev_stack = stack_neighbors(elev)
    elev = np.broadcast_to(elev, elev_stack.shape)
    exclude = np.full_like(elev, False)     #exclude cases with cells rolled around the top and bottom of the map
    exclude[0,-1,:] = True
    exclude[1,0,:] = True
    exclude[4:6,-1,:] = True
    exclude[6:,0,:] = True
    elev_stack = np.where(exclude, elev, elev_stack)    #set excluded cases equal to main elevation so they won't be counted
    dif = elev_stack - elev
    thresh = np.where(np.absolute(dif) > opt('lapse_threshold'),1,0)    #counts cases meeting threshold
    thresh_sum = np.sum(thresh, axis=0)
    has_lapse = np.minimum(thresh_sum,1)
    has_stack = stack_neighbors(has_lapse)
    has_stack = np.where(exclude, 0, has_stack)     #excludes polar wrapping again
    has_sum = np.sum(has_stack, 0)
    all_lapse = np.zeros_like(data)
    verb('    Finding empirical lapse rates in each timestep')
    for t in range(data.shape[0]):
        dat_stack = stack_neighbors(data[t,:,:])
        dat = np.broadcast_to(data[t,:,:], dat_stack.shape)
        lapse = np.where(thresh, (dat_stack-dat)/np.where(thresh,dif,1), 0) #double-check thresh just to avoid div/0 warnings
        lapse = np.sum(lapse, 0)
        lapse = lapse/np.maximum(thresh_sum, 1)     #find average lapse rate in each cell while avoiding div/0 errors.
        lapse_av = np.sum(lapse)/np.sum(has_lapse)   #find average lapse rate for all cells with reported values
        lapse_stack = stack_neighbors(lapse)    #stack again to average resulting lapses with neighbors
        lapse = np.sum(lapse_stack, axis=0)/np.maximum(has_sum,1)
        lapse = np.where(has_lapse, lapse, lapse_av)    #fill in any cells without sufficient data with global average
        all_lapse[t,:,:] = lapse
    return all_lapse
        
#Reads topography from greyscale image and scales to geopotential
# returns elevation array and mask
#   topo_file: topo map file
#   res: target (x,y) size
#   maxel: elevation of highest pixel (m)
#   minel: elevation of lowest pixel (m)
#   sealev: greyscale value of sea level
#   gravity: surface gravitational acceleration (m/s^2)
def Read_topo(topo_map=None, res=None, maxel=None, minel=None, sealev=None, gravity=None):
    if not topo_map:
        topo_map = opt('topo_map')
    if not res:
        res = common('res')
    if not maxel:
        maxel = opt('maxel')
    if not minel:
        minel = opt('minel')
    if not sealev:
        sealev = opt('sealev')
    if not gravity:
        gravity = opt('gravity')

    res = (res[1],res[0])
    verb('     Reading map image')
    hmap = Image.open(topo_map)
    hmap = hmap.convert('F')    # Convert to float for precision on downscaling
    hmapext = hmap.getextrema() # Find minimum and maximum before resizing for accurate elevation scaling
    verb('    Extracting land/sea mask from topo map')
    thresh = (sealev - minel) * (hmapext[1] - hmapext[0]) / (maxel - minel) # Find greyscale value corresponding to sea level
    hmap_bin_ar = np.where(np.asarray(hmap) > thresh, 100, 0) #convert to array, binarize, then convert back to image
    try:
        hmap_bin = Image.fromarray(hmap_bin_ar)
    except:
        verb('     First attempt to convert mask array back to image failed, attempting to read as uint8')
        hmap_bin = Image.fromarray(hmap_bin_ar.astype(np.uint8))    #I dunno why the first version fails for some people, this might help
    hmap_bin = hmap_bin.resize(res, Image.Resampling.BILINEAR)
    topo_bin = np.asarray(hmap_bin)
    topo_mask = np.where(topo_bin > 50, True, False)  # re-binarize after downscaling
    verb('     Resampling and scaling topo elevation')
    hmap = hmap.resize(res, Image.Resampling.BILINEAR)
    elev = np.asarray(hmap)
    elev = elev * (maxel - minel) / (hmapext[1] - hmapext[0]) * gravity #scale to geopotential
    verb(f'     Topo elevation map produced with geopotential range of {np.amin(elev)} to {np.amax(elev)}')
    return elev, topo_mask

#Calculate potential evapotranspiration
# returns PET in mm/month
#   asce-pm method Based on implementation in pyet package but without xarray dependency
#       and asce guidelines https://www.mesonet.org/images/site/ASCE_Evapotranspiration_Formula.pdf 
#   hargreaves method based on here https://www.ctahr.hawaii.edu/oc/freepubs/pdf/EN-106.pdf
#   kalike method based on Koppen's 20 cm/month/C rule for aridity but tuned based on distribution of Pasta bioclimate zones on earth
#   tas: average 2-meter air temperature in C
#       if absent, average maxt and mint, but one or the other required for all methods
#   maxt: daily max temperature in C
#   mint: daily min temperature in C
#       if either absent, saturation vapor pressure calculated from tas
#   rnet: net surface radiation in W/m^2 (absorbed shortwave - emitted longwave)
#       required for asce-pm
#   rin: incident surface radiation in W/m^2
#       required for hargreaves
#   hur: average relative humidity in %
#       required for asce-pm
#   ps: surface pressure in kPa
#   elev: elevation in m
#       alternave to ps for psychrometric constant calc, assuming Earth-like atmosphere curve
#       if both absent, assume constant ps of 101.3
#   wind: average near-surface wind speed in m/s
#       if absent, assume constant of 1 m/s
#   tsoil: deep soil temperature in C or K
#       if absent, assume effectively constant (no heat flux to soil)
#   vegf: surface forest cover
#       if absent, assume constant of 0.5
def Calc_PET(method='asce_pm', tas=None, maxt=None, mint=None, rnet=None, rin=None, hur=None, ps=None, elev=None, wind=None, tsoil=None, vegf=None):
    verb('    Calculating PET')
    if tas is None:
        if not (maxt is None or mint is None):
            tas = maxt + mint / 2
        else:
            raise Exception('Calc_PET requires temperature input')
    if method == 'kalike':
        pet = tas * 7 / 30
    elif method == 'hargreaves':
            if rin is None:
                raise Exception('Calc_PET requires absorbed surface radiation input for hargreaves method')
            rinm = rin * 0.0864 #converts to MJ/m^2*day
            pet = 0.0135 * (tas + 17.8) * rinm * (238.8/(595.5 - 0.55 * tas))
    else:
        if rnet is None:
            raise Exception('Calc_PET requires absorbed surface radiation input for asce_pm or hargreaves methods')
        rnetm = rnet * 0.0864 #converts to MJ/m^2*day
        if hur is None:
            raise Exception('Calc_PET requires relative humidity input for asce_pm method')
        if ps is None:
            verb('     No surface pressure found in Calc_PET, using backup')
            ps = opt('pet_backup_ps')
            if elev is not None:
                ps *= ((293 - 0.0065*elev) / 293)**5.26     #backup options for surface pressure; probably don't work amazingly
        if wind is None:
            verb('     No wind found in Calc_PET, using backup')
            wind = opt('pet_backup_wind')    #simple placeholder just to not make this an absolute requirement
        if tsoil is None:
            verb('     No soil temperature found in Calc_PET, assuming constant')
            gsoil = 0
        else:
            gsoil = 0.07 *(np.roll(tsoil, 1, 0) - np.roll(tsoil, 1, 0))     #estimate net heat flux to soil by soil temperature change

        
        
        lambd = 2.501 - 0.002361 * tas    #latent heat of water vaporization, MJ/kg
        #asce-pm guidelines say to set lambd as constant, but I'm leaving it in
        
        #if not ts:
        #    gamma = 0.000665 * ps #backup option if for whatever reason I remove temperature dependency
        #else:
        cp = 0.001013   #specific heat of air, MJ/kg*K
        rat = 0.622 * 287.05 / opt('pet_gascon')    #ratio of molecular weights of water and dry air
        gamma = cp * ps / (rat * lambd)    #psychrometric constant, kPa/K
        

        es = 0.6108 * np.exp(17.27 * tas / (tas + 237.3)) #saturation pressure at constant t, kPa
        delt = 4098 * es / (tas + 237.3)**2  #slope of saturation curve, kPa/K

        if not (maxt is None or mint is None):
            esmax = 0.6108 * np.exp(17.27 * maxt / (maxt + 237.3))
            esmin = 0.6108 * np.exp(17.27 * mint / (mint + 237.3))
            esav = (esmax + esmin)/2
        else:
            verb('     No daily maxt and mint found in Calc_PET, ignoring diurnal temp variation')
            esav = es    #average saturation pressure, kPa
        ea = esav * hur / 100      #actual vapor pressure, kPa
        de = es - ea    #vapor pressure deficit, kPa
        if method == 'asce-pm': #may add option for more regular pm, which excludes this
            if vegf is None:    #if no vegetation input, split the difference
                verb('     No vegetation found in Calc_PET, using backup')
                vegf = 0.5
            cn = (900 + 500 * vegf)     #interpolate between tall and short reference vegetation based on forest cover
            cd = (0.34 + 0.04 * vegf)
            ga = wind * cn/(tas + 273.15) #wind factor
            denom = 1 + cd * wind
        else:
            ga = 1 + 0.537 * wind
            denom = 1
        pet = ((delt * (rnetm - gsoil) / lambd) + gamma * ga * de) / (delt + gamma * denom) #calculates as mm/day, for ease of use with asce coefficients
        
    pet = np.where(pet < 0, 0, 30*pet)  #clamp to positive values and convert to mm/month
    return pet

#Calculate growing degree-days
# returns GDD/month
#   tas: 2-meter air temperature in C
#   base: baseline temperature to count up from in C
#   plat_start: start of 'plateau'; no additional GDD is counted for degrees above this temp
#   plat_end: end of 'plateau'; GDD declines with even higher temps
#   comp: compensation temperature; zero GDD above this temp
#       with linear GDD-temp relation between plat_end and comp
def Calc_GDD(tas, base=5, plat_start=30, plat_end=10000, comp=20000):
    verb('     Calculating monthly GDD')
    GDD_max = plat_start - base
    back_slope = GDD_max / (comp - plat_end)

    GDD = tas - base    #count degrees from base
    GDD = np.where(tas > plat_start, GDD_max, GDD)  #flatten after plateau start
    GDD = np.where(tas > plat_end, GDD_max - back_slope * (tas - plat_end), GDD)    #decline after plateau end
    GDD = np.maximum(GDD, 0)    #clamp to positive

    GDD *= 30 * opt('month_length') * opt('gdd_productivity_modifier')    #convert to GDD/month and apply modifiers

    return GDD

#Calculate total gdd for year, and optionally totally gint
#   gint: options growth interruption
#   cont: count longest contiguous accumulation of GDD, rather than simple total
#   inf: where all months have GDD, set total to 1 million, treating total as effectively infinite
# full procedure works through 3 loops, each going through year twice (so accumulation carries from last month to first)
# first, GInt is accumulated forward, adding each month to the total but interrupted where GInt falls to 0
# then, the resulting totals are propogated backwards in the GInt count,
#  such that every month in each contiguous interruption period shows the total GInt accumulation for that period
# then, GDD is accumulated forward, interrupted only where GDD falls to zero during a period of sufficiently great GInt accumulation
def Calc_GDD_total(gdd, gint=None, cont=None, inf=None, th_gi=1250):
    verb('     Calculating total GDD count')
    if cont is None:
        cont = opt('gdd_require_contiguous')
    if inf is None:
        inf = opt('gdd_indicate_inf')
    if cont and len(gdd) > 1:   #check for largest contiguous accumulation of gdd rather than total; skip if seasonless
        if gint is None:
            gint_acc = np.ones_like(gdd) * 1e6  #set high so that it's always high enough to interrupt GDD
        else:
            gint_acc = np.copy(gint)
            for i in range(2):  #loop through year twice so last month loops into first
                for t in range(len(gint)):
                    gint_acc[t,:,:] = np.where(gint[t,:,:] > 0, gint[t,:,:] + gint_acc[t-1,:,:], 0)    #accumulate gint forward
            if inf:
                gint_acc[-1,:,:] = np.where(np.amin(gint_acc,0) > 0, 1e6, gint_acc[-1,:,:]) #set last month to 1 million to show effective infinity
            for i in range(2):
                for t in range(len(gint)):
                    tn = len(gint) - (t+1)
                    gint_acc[tn-1,:,:] = np.where(np.minimum(gint[tn,:,:],gint[tn-1,:,:]) > 0, gint_acc[tn,:,:], gint_acc[tn-1,:,:])    #propogate total of each gint period backwards to rest of period
            gint_tot = np.amax(gint_acc, 0)
        gdd_acc = np.copy(gdd)
        for i in range(2):
            for t in range(len(gdd)):
                gdd_n = gdd[t,:,:] + gdd_acc[t-1,:,:]
                gdd_acc[t,:,:] = np.where(gdd[t,:,:] > 0, gdd_n, np.where(gint_acc[t,:,:] > th_gi, 0, gdd_n))   #accumulate gdd forward, interrupting only in large gint periods
        gdd_tot = np.amax(gdd_acc, 0)

    else:
        gdd_tot = np.sum(gdd, 0)    #otherwise just sum total
        gdd_acc = gdd   # for below inf check
        if gint is not None:
            gint_tot = np.sum(gint, 0)
            if inf:
                gint_tot = np.where(np.amin(gint,0) > 0, 1e6, gint_tot)
    if inf:
        gdd_tot = np.where(np.amin(gdd_acc,0) > 0, 1e6, gdd_tot)  #where there is growing in all months, set gdd to 1 million to indicate effective infinity
    if gint is None:
        return gdd_tot
    else:
        return gdd_tot, gint_tot

#Estimate evaporation from precipitation and PET using simple soil water model
# returns monthly evaporation in mm/month
#   pet: total monthly potential evapotranspiration in mm/month
#   pr: total monthly precipitation in mm/month
def Estimate_evap(pet, pr):
    verb('    Estimating evapotranspiration from PET and precipitation')
    evap = np.zeros_like(pet)   #evaporation
    soilw = np.zeros_like(pet)  #soil water at end of each month
    diff = np.zeros_like(pet[0,:,:]) + 1000     #difference in soil water between first and final states
    surpdef = pet - pr  #surplus or deficit of pet over pr each month
    while np.amax(diff) > 10:   #iterate through year until there's less than a 1 cm discrepency in starting soil water anywhere in the world
        initsw = soilw[-1,:,:]
        for t in range(pet.shape[0]):
            sev = np.minimum(soilw[t-1,:,:], surpdef[t,:,:] * np.minimum(soilw[t-1,:,:],250)/250)   #soil evaporation rate limited by soil water content and saturation when below 25 cm
            evap[t,:,:] = np.where(surpdef[t,:,:] > 0, pr[t,:,:] + sev, pet[t,:,:])     #evaporation equal to pet where pr exceeds it, pr + soil evaporation otherwise
            soilw[t,:,:] = np.minimum (500, soilw[t-1,:,:] - np.where(surpdef[t,:,:] > 0, sev, surpdef[t,:,:]))     #soil water adjusted but limited to 50 cm
        diff = np.absolute(soilw[-1,:,:] - initsw)
    return evap          

#Pull data from exoplasim output netcdf file and process as necessary
# returns the processed array
#   dat: list of input netcdf data objects; output will be averaged across files
#   key: data key in file
#   coords: (lat, lon) arrays for interpolation
#   single: only extract once, from first timestep in first file
#   no_interp: don't allow interpolation 
#   adjust: adjustment map to be applied to data (for e.g. temp adjustment by topography)
#   dummy_ice: use dummy_ice option for interpolation
#   low: pick last atmospheric layer in (time, layer, lat, lon) array, i.e. near-surface layer in eps output
#   bin_ext: ext option for bin (-1 min, 0 avg, 1 max)
def Get_nc(dat, key, coords=None, res=None, single=False, no_interp=False, adjust=None, dummy_ice=False, low=False, bin_ext=0):
    verb(f'    Extracting {key}')
    if single:
        try:
            if low:
                dat_ar = dat[0][key][0,-1,:,:]
            else:
                dat_ar = dat[0][key][0,:,:]
            dat_ar = np.expand_dims(dat_ar, 0)  #keep time dimension so it's present for other functions
        except:
            dat_ar = dat[0][key][:]
    elif opt('file_combine') == 'seq':  #link data from each file along time dimension into single long array
        if low:
            dat_ar = [d[key][:,-1,:,:] for d in dat]
        else:
            dat_ar = [d[key][:] for d in dat]
        dat_ar = np.concatenate(dat_ar, 0)
    else:
        dat_ar = None
        for d in dat:
            if low:
                d_ar = d[key][:,-1,:,:]
            else:
                d_ar = d[key][:]
            if dat_ar is not None:
                dat_ar += d_ar
            else:
                dat_ar = d_ar
        dat_ar /= len(dat)  #sum values from all input files and then divide by file number to average
    if opt('bin_months') > 1 and not single:
        verb('     Binning data')
        dat_ar = Bin_months(dat_ar, opt('bin_months'), ext=bin_ext)
    if opt('interp_scale') and not no_interp:
        dat_ar = Interp(dat_ar, interp_type=opt('interp_type'), coords_in=coords, res=res, dummy_ice=dummy_ice)
    if adjust is not None:
        verb('     Applying adjustment')
        dat_ar += adjust
    return dat_ar

#If data_key already in data, return from data, otherwise get from data key
# also returns whether key was in data
# useful for sea functions, to avoid redundant work pulling data used for both land and sea
def Get_nc_if(dat, dat_key, data, data_key=None, coords=None, res=None, single=False, no_interp=False, adjust=None, dummy_ice=False, low=False, bin_ext=0):
    if not data_key:
        data_key = dat_key
    if data_key in data:
        verb(f'    {data_key} already extracted')
        return data[data_key], True
    else:
        return Get_nc(dat, dat_key, coords=coords, res=res, single=single, no_interp=no_interp, adjust=adjust, dummy_ice=dummy_ice, low=low, bin_ext=bin_ext), False

#Use temperature, file elevation, and full-resolution topography to create temperature adjustment map
# returns both the processed temp array and the adjustment map array, without adding them
#   dat: list of input netcdf data objects; output will be averaged across files
#   t_key: temperature data key in file
#   g_key: elevation data key in file
#   topo: full-res topography data array
#   t_unadjust: return temperature data without adjustment applied
def Get_nc_adjust(dat, t_key, g_key, coords=None, res=None, topo=None, t_unadjust=False, bin_ext=0):
    t_ar = Get_nc(dat, t_key, coords=coords, res=res, bin_ext=bin_ext)
    if not res and opt('interp_scale'):
        res = get_res((t_ar.shape[1],t_ar.shape[2]), scale=opt('interp_scale'))
    if not topo:
        if not res or opt('topo_map') is None:
            adjust = np.zeros_like(t_ar)    #return zero adjustment if not using topography map
            return t_ar, adjust
        else:
            try:
                topo = common('topo')   #check if topo already produced to avoid redundant work
                verb('     Topo map already uploaded; reusing')
            except:
                try:
                    print("  Uploading Higher-Resolution Topography")
                    topo, mask = Read_topo()
                    add_common('topo', topo)
                    add_common('mask_topo', mask)
                except:
                    print("  Unable to upload topography map; proceeding without temperature adjustment")
                    adjust = np.zeros_like(t_ar)
                    return t_ar, adjust
    
    try:
        ground = common('ground')
        elev_dif = common('elev_dif')
    except:
        ground = Get_nc(dat, g_key, single=True, no_interp=True)
        ground = ground[0,:,:]
        verb(f'     Comparing {g_key} and topo to find necessary adjustment')
        ground_big = Interp(ground, interp_type=opt('interp_type'), coords_in=coords, res=res)
        elev_dif = topo - ground_big
        add_common('ground', ground)
        add_common('elev_dif', elev_dif)
        if opt('debug_file'):
            add_common('ground_big', ground_big)


    clapse = opt('const_lapse_rate')
    if clapse is not None and clapse > 0:
        verb(f'    Applying constant lapse rate of {clapse} K/km')
        lapse = np.ones_like(t_ar) * -clapse / (1000*opt('gravity'))
    else:
        verb(f'    Extracting {t_key} at original scale for calculating lapse rate')
        t_ar_sm = Get_nc(dat, t_key, no_interp=True)
        verb(f'    Finding lapse rate from {t_key} and {g_key}')
        lapse = Find_lapse(t_ar_sm, ground)
        lapse = Interp(lapse, interp_type=opt('interp_type'), coords_in=coords, res=res)
    
    adjust = lapse * np.expand_dims(elev_dif, 0)

    if opt('debug_file'):
        add_common('lapse', lapse)

    if not t_unadjust:
        t_ar += adjust
    return t_ar, adjust

#standard function to extract pet and optionally evap from ExoPlaSim outputs
# dat: dat to extract parameters from
# data: previously extracted parameters
# pet_method: method for estimating pet
# no_interp: produce without interpolation
def Get_pet(dat, data, pet_method=None, no_interp=False):
    verb('    Gathering data for PET calculation')

    if not pet_method:
        pet_method = opt('pet_method')

    if 'adjust' in data and not no_interp:
        adjust = data['adjust']
    else:
        adjust = 0
    tas,w = Get_nc_if(dat, 'tas', data, adjust=adjust+opt('temp_adjust'), no_interp=no_interp)
    if pet_method == 'kalike':
        pet = Calc_PET(method = pet_method,
                       tas=tas)

    elif pet_method == 'hargreaves':
        rss,w = Get_nc_if(dat, 'rss', data, no_interp=no_interp)
        pet = Calc_PET(method = pet_method,
                       tas=tas,
                       rin=rss)
        
    else:
        rss,w = Get_nc_if(dat, 'rss', data, no_interp=no_interp)
        rls,w = Get_nc_if(dat, 'rls', data, no_interp=no_interp)
        hur,w = Get_nc_if(dat, 'hur', data, low=True, no_interp=no_interp)
        ps,w = Get_nc_if(dat, 'ps', data, no_interp=no_interp)
        spd,w = Get_nc_if(dat, 'spd', data, low=True, no_interp=no_interp)
        maxt,w = Get_nc_if(dat, 'maxt', data, 'maxt', adjust=adjust+opt('temp_adjust'), no_interp=no_interp)
        mint,v = Get_nc_if(dat, 'mint', data, 'mint', adjust=adjust+opt('temp_adjust'), no_interp=no_interp)
        if not w or not v:
            ts,w = Get_nc_if(dat, 'ts', data, adjust=adjust+opt('temp_adjust'), no_interp=no_interp)
            tdif = tas - ts#adjustment from surface to 2-meter to apply to maxt and mint
            maxt += tdif
            mint += tdif
        tsoil,w = Get_nc_if(dat, 'tso', data, no_interp=no_interp)
        tsoil = np.where(tsoil > -200-opt('temp_adjust'), tsoil+adjust, tsoil) #don't adjust tsoil over seas, which are filled in as 0
        if opt('pet_use_vegf'):
            try:
                vegf = Get_nc(dat, 'vegf', no_interp=no_interp)
            except:
                vegf = None
        else:
            vegf = None
        pet = Calc_PET(method = opt('pet_method'),
                      tas=tas,
                      maxt=maxt,
                      mint=mint,
                      rnet=rss + rls,   #total net radiation from net longwave and shortwave
                      hur=hur,
                      ps=ps/10,    #convert from hPa to kPa   
                      wind=spd,
                      tsoil=tsoil, #no need to convert units because only the difference across months matters
                      vegf=vegf)

    return pet

#get mask from common if it exists, or add it using key if it doesn't
#   convert: convert from 3d binary to 2d boolean
def get_mask(dat, key, convert=False):
    try:
        mask = common('mask')
    except:
        verb('    Retrieving land/sea mask')
        mask = Get_nc(dat, key, single=True, no_interp=True)
        if opt('sea_type') == 'sea_none':
            mask = np.ones_like(mask)
        if opt('interp_scale'):
            verb('    Interpolating mask to output resolution')
            mask_big = Interp(mask, interp_type='nearest')  #save upscaled mask as well if using interpolation
            if convert:
                mask_big = np.where(mask_big[0,:,:] > 0.5, True, False)
            add_common('mask_big', mask_big)
        if convert:
            mask = np.where(mask[0,:,:] > 0.5, True, False)
        add_common('mask', mask)
    return mask

#Produce debug netcdf file with arrays in data, param, clim, and common
# returns debug file name
def Debug_file(data, params):
    if opt('outname'):
        name = opt('outname') + '_debug.nc'
    else:
        name = 'output_debug.nc'
    debug = nc.Dataset(name, 'w')
    for k, v in data.items():   #I'm starting to wonder if there's a better way to just pull the first entry in a dictionary, but eh
        sh = v.shape
        break
    lat,lon = get_coords((sh[1], sh[2]))
    lat *= 180/math.pi
    lon = lon * 180/math.pi
    time = np.arange(0,sh[0],1)
    if opt('interp_scale'):
        names = ('lat_small','lon_small','time')
    else:
        names = ('lat','lon','time')
    for dim, n in zip((lat, lon, time), names):
        newdim = debug.createDimension(n, len(dim))
        newvar = debug.createVariable(n,'f8',n,compression='zlib', complevel=9)
        newvar[:] = dim[:]
    lat_big = None
    data.update(params)
    verb(f'    Saving {len(data)} data arrays to debug file')
    for k, v in data.items():
        verb(f'     Saving {k}')
        latdim = names[0]
        londim = names[1]
        try:
            if v.shape[-1] != len(lon):
                if lat_big is None:
                    lat_big, lon_big = get_coords((v.shape[-2], v.shape[-1]), big=True)
                    lat_big *= 180/math.pi
                    lon_big = lon_big * 180/math.pi
                    for dim, n in zip((lat_big, lon_big), ('lat','lon')):
                        newdim = debug.createDimension(n, len(dim))
                        newvar = debug.createVariable(n,'f8',n,compression='zlib', complevel=9)
                        newvar[:] = dim[:]
                latdim = 'lat'
                londim = 'lon'
            if v.ndim == 3:
                newvar = debug.createVariable(k,'f8',('time',latdim,londim),compression='zlib', complevel=9)
                newvar[:(v.shape[0]),:,:] = v[:]
            else:
                newvar = debug.createVariable(k,'f8',(latdim,londim),compression='zlib', complevel=9)
                newvar[:] = v[:]
        except:
            print(f'  Could not save {k} to debug file')
    verb(f'    Searching for extra common data arrays for debug file')
    for m in ('mask','mask_big','mask_topo','topo','ground','ground_big','elev_dif','lapse'):
        try:
            var = common(m)
            verb(f'     Saving {m}')
            if m in ('mask','ground'):
                latdim = names[0]
                londim = names[1]
            else:
                latdim = 'lat'
                londim = 'lon'
            if m == 'lapse':
                newvar = debug.createVariable(m,'f8',('time',latdim,londim),compression='zlib', complevel=9)
            else:
                newvar = debug.createVariable(m,'f8',(latdim,londim),compression='zlib', complevel=9)
            newvar[:] = var[:]
        except:
            continue
    options = debug.createGroup('options')
    verb(f'    Saving {len(kpasta_options)} options to debug file')
    for k,v in kpasta_options.items():
        if v is None:
            va = "None"
        elif isinstance(v, bool):
            if v:
                va = "True"
            else:
                va = "False"
        else:
            va = v
        setattr(options, k, va)
    debug.close()
    return name

#Alternate data function to run if using unusual input files
# As an example, this is set up for retrieving data from terraclimate data on Earth
# but you can change this or replace it with another function
# and it will be run instead of the climate-specific data function if force_alt_data is true,
# ignoring any usual file inputs
def Alternate_Data():
    dat = nc.Dataset('TerraClimate19812010_tmax.nc')
    maxt = dat['tmax'][:]
    dat.close()
    
    dat = nc.Dataset('TerraClimate19812010_tmin.nc')
    mint = dat['tmin'][:]
    dat.close()
    
    dat = nc.Dataset('TerraClimate19812010_ppt.nc')
    ppt = dat['ppt'][:]
    dat.close()
    
    dat = nc.Dataset('TerraClimate19812010_aet.nc')
    aet = dat['aet'][:]
    dat.close()
    
    dat = nc.Dataset('TerraClimate19812010_pet.nc')
    pet = dat['pet'][:]
    dat.close()

    tas = (maxt + mint) / 2

    all_data = dict(
        tas=tas,
        mint=mint,
        maxt=maxt,
        pr=ppt,
        evap=aet,
        pet=pet,
        land=np.where(mint==mint[0,0,0],0.0,1.0)
        )

    
    #mask = np.where(mint==mint[0,0,0],False,True)
    #add_common('mask', mask[0,:,:])
    fac = 8
    for k,v in all_data.items():
        new1 = np.empty((12,int(4320/fac),8640), dtype=v.dtype)
        new2 = np.empty((12,int(4320/fac),int(8640/fac)), dtype=v.dtype)
        for y in range(int(4320/fac)):
            new1[:,y,:] = np.nanmean(v[:,y*fac:(y+1)*fac,:],1)
        for x in range(int(8640/fac)):
            new2[:,:,x] = np.nanmean(new1[:,:,x*fac:(x+1)*fac],2)
        if k == 'land':
            mask = np.where(np.isnan(new2[0,:,:]), False, True)
            mask = np.where(new2[0,:,:]>0.5,True,False)
            add_common('mask', mask)
        all_data[k] = np.where(np.isnan(new2), 0, new2)
    
    return all_data
    

#Standard function to retrieve data and determine parameters
# returns dictionary of climate parameters
#   files: list of netcdf files
#   land_funcs, sea_funcs: lists of appropriate climate functions from Clim_func
def Get_params(files, land_funcs, sea_funcs):
    if not opt('force_alt_data'):
        if len(files) > 1:
            print(f" Extracting data from {files[0]} et al...")
        else:
            print(f" Extracting data from {files[0]}...")
    if opt('file_combine') == 'param' and not opt('force_alt_data'):
        params = {}
        if len(files) > 1:
            print("  Extracting and processing data to climate parameters per-file before averaging")
        if opt('seasonless'):
            print("  Averaging data across months to produce seasonless climate")
        for f in files:         # determine parameters for each year before averaging together
            verb(f'   Extracting from {f}')
            dat = nc.Dataset(f)
            try:
                coords_from_file(dat,'lat','lon') #try to ensure coords read from file for eps inputs
            except:
                pass
            verb('   Extracting data for land')
            data = land_funcs[0]([dat])
            verb('   Extracting data for sea')
            data = sea_funcs[0]([dat], data)    # sea functions take data output from land function and add to it
            verb('   Extracting any necessary extra data')
            data = Extra_Data([dat], data)  # extra data functions run in all cases
            if opt('seasonless'):
                for k, v in data.items():
                    if v.ndim > 2:
                        data[k] = np.mean(v, 0, keepdims=True)  # for seasonless zones, average along time axis but keep dimension
            verb(f'   Extracted data contains {[k for k,v in data.items()]}')
            if len(files) == 1:
                print(" Processing data to climate parameters...")
            verb('   Calculating parameters for land')
            par = land_funcs[1](data)
            verb('   Calculating parameters for sea')
            par = sea_funcs[1](data, par)
            verb('   Calculating any necessary extra parameters')
            Extra_Param(data, par)  #extra parameter functions
            if params:
                for k,p in params.items():
                    params[k] = params[k] + par[k]
            else:
                params = par
            verb(f'   Calculated parameters are {[k for k,v in params.items()]}')
        for k,p in params.items():
            params[k] = params[k] / len(files)
        
    else:
        if opt('force_alt_data'):
            print(" Using alternate data collection function...")
            if opt('file_combine') != 'data':
                print("  Note: file combination method may not be applied")
            data = Alternate_Data()
        else:
            if opt('file_combine') == 'seq' and len(files) > 1:
                print("  Linking data across files into single year")
            dats = [nc.Dataset(f) for f in files]     # average data across years, then determine parameters
            #try:
            coords_from_file(dats[0],'lat','lon') #try to ensure coords read from file for eps inputs
            #except:
            #    pass
            verb('   Extracting data for land')
            data = land_funcs[0](dats)
            verb('   Extracting data for sea')
            data = sea_funcs[0](dats, data)
            verb('   Extracting any necessary extra data')
            data = Extra_Data(dats, data)  
        if opt('seasonless'):
            print("  Averaging data across months to produce seasonless climate")
            for k, v in data.items():
                if v.ndim > 2:
                    data[k] = np.mean(v, 0, keepdims=True)
        verb(f'   Extracted data contains {[k for k,v in data.items()]}')
        if not opt('force_alt_data') and len(files) > 1:
            print(" Processing averaged data from all files to climate parameters...")
        else:
            print(" Processing data to climate parameters...")
        verb('   Calculating parameters for land')
        params = land_funcs[1](data)
        verb('   Calculating parameters for sea')
        params = sea_funcs[1](data, params)
        verb('   Calculating any necessary extra parameters')
        params = Extra_Param(data, params)
        verb(f'   Calculated parameters are {[k for k,v in params.items()]}')
    if opt('debug_file'):
        print(" Making debug file...")
        if opt('file_combine') == 'param' and len(files) > 1 and not opt('force_alt_data'):
            print("  Note: only contains data from latest file")
        debug_name = Debug_file(data, params)
        print(f"  Saved to {debug_name}")
    return params

#Standard function to determine climate zones for map
# returns dictionary of climate maps
#   params: parameters from Get_Params function
#   land_funcs, sea_funcs: lists of appropriate climate functions from Clim_func
def Get_clims(params, land_funcs, sea_funcs):
    print(" Determining climate zones...")
    try:
        mask = params['mask']   #check if mask is provided
        verb('    Using provided land/sea mask')
    except:
        try:
            mask = common('mask_topo')
            verb('   Using land/sea mask determined from topo map')
        except:
            try:
                mask = common('mask_big')   # if not, check common, looking for topo mask, then big mask
                verb('   Using upscaled land/sea mask')
            except:
                try:
                    mask = common('mask')
                    verb('   Using saved land/sea mask')
                except:
                    for k, v in params.items():      # if no mask is found, construct a universally true mask based on the first at-least 2-dimensional array in params
                        if v.ndim >= 2:
                            mask = np.full((v.shape[-2],v.shape[-1]), True)
                            print('  No land/sea mask provided to Get_clims; assuming all land')
                            break
    land_clims = np.zeros(mask.shape,dtype=np.uint16)
    sea_clims = np.zeros(mask.shape,dtype=np.uint16)
    par = {}
    if opt('blend') and opt('sea_type') != ('sea_none'):    #make masks of where to find land and sea climates
        do_land = np.where(mask, True, False)
        do_sea = np.where(do_land, False, True)
    else:
        do_land = np.full(mask.shape,True)
        do_sea = do_land
    if opt('efficient'):    # use efficient option to run algorithm on whole array at once
        verb('   Using efficient classification functions')
        land_clims = land_funcs[2](params)
        sea_clims = sea_funcs[2](params)
    else:
        verb('   Iterating through map to classify climates')
        for y in range(mask.shape[0]):     # for most cases, iterate over each cell and run function individually
            for x in range(mask.shape[1]):
                for k, v in params.items():
                    try:
                        par[k] = v[y,x]
                    except:
                        continue
                if do_land[y,x]:
                    land_clims[y,x] = land_funcs[2](par)
                if do_sea[y,x]:
                    sea_clims[y,x] = sea_funcs[2](par)
    maps = {}
    if opt('make_chart'):
        print(" Making climate chart...")
        if opt('outname'):
            outname = opt('outname')
        else:
            outname = 'output'
        for k,v in zip(('land','sea'), (land_clims,sea_clims)):
            verb(f'   Making {k} chart')
            chart_im = Make_chart(v,params)
            if chart_im is not None:
                chart_im.save(outname + '_' + k + '_chart.png')
                print("  Saved to " + outname + '_' + k + '_chart.png')
    if opt('blend'):
        verb('   Blending land and sea maps to single image')
        clims = np.where(do_land, land_clims, sea_clims)
        maps['full'] = clims
    else:
        verb('   Producing separate land and sea maps')
        maps['land'] = land_clims
        maps['sea'] = sea_clims
    return maps 

## Image processing

#Add colors in dictionary to colmap array, extending as necessary:
def Add_color(colmap, colors):
    verb(f'    Adding {len(colors)} colors to color map')
    extra = int(max(colors)) - (len(colmap) - 1)
    if extra > 0:
        colmap = np.pad(colmap, ((0,extra),(0,0)))
    for k, v in colors.items():
        if type(v) is int:
            colmap[k] = colmap[v]   #allows zones to copy from previously entered colors
        else:
            colmap[k] = v
    return colmap

#Create list of colors corresponding to each climate type, depending on the color type
def Make_colmap(land_type=None, land_color=None, sea_type=None, sea_color=None, color_file=None):
    if not land_type:
        land_type = opt('land_type')
    if not land_color:
        land_color = opt('land_color')
    if not sea_type:
        sea_type = opt('sea_type')
    if not sea_color:
        sea_color = opt('sea_color')
    if not color_file:
        color_file = opt('color_file')
        
    colmap = np.zeros((1000,3),dtype=np.uint8)
    colmap = Add_color(colmap, def_color)   # default colors always added
    
    if land_color == ('file') or sea_color == ('file'):     # add colors from list first so later additions override it
        verb('    Reading color file')
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        cfg.read(color_file)
        file_colors = {}
        for section in cfg.sections():
            for k, v in cfg.items(section):
                co = v.split(',')
                co = [int(j) for j in co]
                co = list(co)
                file_colors.update({globals()[k]: co}) # check global variables to match each key in the list to the climate zone variable
        colmap = Add_color(colmap, file_colors)

    if land_color != ('file'):
        verb('    Loading default color lists')

        if land_type == ('Koppen-Geiger') or land_type == ('KG_unproxied') or land_type == ('TwoParamKG'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, koppen_standard_color)
            elif land_color == ('alt_red'):
                colmap = Add_color(colmap, koppen_alt_red_color)
            elif land_color == ('true'):
                colmap = Add_color(colmap, koppen_true_color)
            colmap = Add_color(colmap, koppen_group_color)
            colmap = Add_color(colmap, koppen_twoletter_color)
            colmap = Add_color(colmap, koppen_reduced_color)
            
        elif land_type == ('Trewartha'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, trewartha_standard_color)
            elif land_color == ('kalike'):
                colmap = Add_color(colmap, trewartha_kalike_color)
            colmap = Add_color(colmap, trewartha_group_color)
            
        elif land_type == ('Holdridge'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, holdridge_standard_color)
                if sea_color == ('white'):
                    colmap = Add_color(colmap, {HPolDesert: [240,240,240]}) #special case so polar desert doesn't blend into white seas
            elif land_color == ('vibrant'):
                colmap = Add_color(colmap, holdridge_vibrant_color)

        elif land_type == ('Thornthwaite'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, thornthwaite_standard_color)
                colmap = Add_color(colmap, thornthwaite_variability_color)

        elif land_type == ('Whittaker'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, whittaker_standard_color)

        elif land_type == ('Woodward'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, woodward_standard_color)

        elif land_type == ('IPCC'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, ipcc_standard_color)

        elif land_type == ('WCR'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, wcr_standard_color)

        elif land_type == ('Prentice'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, prentice_standard_color)
                
        elif land_type == ('Pasta'):
            if land_color == ('standard'):
                colmap = Add_color(colmap, pasta_standard_color)
            elif land_color == ('earthlike'):
                colmap = Add_color(colmap, pasta_earthlike_color)
            elif land_color == ('true'):
                colmap = Add_color(colmap, pasta_true_color)
                
    if sea_color != ('file'):
        if sea_type == ('sea_standard'):
            if sea_color == ('standard'):
                colmap = Add_color(colmap, sea_standard_color)
            elif sea_color == ('true'):
                colmap = Add_color(colmap, sea_true_color)
            elif sea_color == ('white'):
                colmap = Add_color(colmap, sea_white_color)
        elif sea_type == ('sea_Pasta'):
            if sea_color == ('standard'):
                colmap = Add_color(colmap, pastaocean_standard_color)
            elif sea_color == ('true'):
                colmap = Add_color(colmap, pastaocean_true_color)

    return colmap

#Get colmap out of common if it is already made, otherwise make it
def Get_colmap():
    try:
        colmap = common('colmap')
    except:
        print('  Constructing color map...')
        colmap = Make_colmap()
        add_common('colmap', colmap)
    return colmap

#Attempt to produce image key of all climates in maps using colmap
def Make_key(maps, colmap=None):

    if colmap is None:
        colmap = Get_colmap()

    verb('     Finding all used climate zones')
    keys = []
    for k, v in maps.items():   #finds all climate zone types used in maps
        for x in np.nditer(v):
            if x not in keys:
                keys.append(x)
    
    verb('   Finding names for climate zones')
    zones = {}
    glob = list(globals().items())
    for k, v in glob:    #attempts to find names for each climate by searching global variables
        if type(v) is int and v in keys:
            try:
                zones[name_key[v]] = v
            except:
                zones[k] = v
    
    verb(f'   Found {len(keys)} zones, constructing key image')
    key_im = Image.new(mode="RGB", size=(1000,len(zones)*25+5), color=(0,0,0))
    im = ImageDraw.Draw(key_im)
    font_size = opt('font_size')
    font = ImageFont.load_default(size=font_size) #load default font
    h = 0
    for k, v in zones.items():
        im.rectangle([5, h+5, 25, h+25], fill=tuple(colmap[v]), outline=(100,100,100), width=1)     #draw box in zone color
        im.text((30, h+(12-int(font_size/2))), k, fill=(255,255,255), font=font)   #write zone name
        h += 25
    im_dat = np.asarray(key_im)
    col = 999
    while True:
        if np.amax(im_dat[:,col,:]) > 0 or col < 0: #scan from the right edge to find first column with color in it, to find edge of text
            break
        col -= 1
    key_im = key_im.crop((0,0,col+6,key_im.height))   #crop to 5 pixels right of longest text
    return key_im

#Make chart of points in map by average temp and total precipitation, with each point colored by zone
def Make_chart(clim, par):
    try:
        ta = par['Avg_Temp']
    except:
        print('  No appropriate average temp found in parameters; skipping chart output')
        return None
    try:
        pr = par['Total_Precip']/10   #convert to cm/year
    except:
        print('  No appropriate precipitation found in parameters; skipping chart output')
        return None
        
    maxt = math.ceil(np.amax(np.where(clim==0,ta[0,0],ta))/10)*10   #scale chart dimensions based on range of parameters
    mint = math.floor(np.amin(np.where(clim==0,ta[0,0],ta))/10)*10
    pmax = math.ceil(np.amax(np.where(clim==0,pr[0,0],pr))/100)*100
    pmin = math.floor(np.amin(np.where(clim==0,pr[0,0],pr))/100)*100

    tlen = int((maxt-mint)*10)
    plen = int(pmax-pmin)
    
    colmap = Get_colmap()

    verb('     Constructing chart image')

    font_size = opt('font_size')
    font = ImageFont.load_default(size=font_size) #load default font

    ex = int(max(0, font_size/1.5 * len(str(pmax)) - 30))
    
    chart_im = Image.new(mode="RGB", size=(ex+60+tlen, 70+plen), color=(0,0,0))
    im = ImageDraw.Draw(chart_im)



    

    im.line ([(ex+38,19),(ex+38,20+plen),(ex+40+tlen,20+plen)], fill=(255,255,255), width=2) #chart axes

    for p in range(int(1+plen/100)):
        im.line([(ex+33,19+p*100),(ex+38,19+p*100)], fill=(255,255,255), width=1)
        im.text((ex+31,18+p*100), str(pmax-p*100), fill=(255,255,255), font=font, anchor='rm')    #precip markings
    for t in range(int(1+tlen/100)):
        im.line([(ex+40+t*100,20+plen),(ex+40+t*100,26+plen)], fill=(255,255,255), width=1)   #temp markings
        im.text((ex+40+t*100,29+plen), str(mint+t*10), fill=(255,255,255), font=font, anchor='mt')
    

    tco = 40 + np.round((ta-mint)*10)   #map coordinates
    pco = 19 + plen - np.round((pr-pmin))

    verb('    Charting points on map by average climate')
    points = {}
    for x in range(clim.shape[0]):
        for y in range(clim.shape[1]):
            if clim[x,y] > 0:
                ind = str(int(tco[x,y]))+','+str(int(pco[x,y]))
                try:
                    points[ind].append(clim[x,y])   #assemble list of climates on each point
                except:
                    points[ind] = [clim[x,y]]
    verb(f'     Charted {clim.shape[0] * clim.shape[1]} map points to {len(points)} chart pixels, coloring by most common climate zone')
    for k, v in points.items():
        co = k.split(',')
        if len(v) == 1:
            col = tuple(colmap[v[0]])   #if one color on point, use that
        else:
            alts = [0]
            curmax = 0
            for c in v:
                if c not in alts:
                    count = v.count(c)  #count occurances of each climate type and assemble list of most common types
                    if count > curmax:
                        alts = [c]
                        curmax = count
                    elif count == curmax:
                        alts.append(c)
            col = tuple(colmap[alts[int(len(alts)/2)]]) #if one winner, use that climate, if multiple, pick from middle of list
        im.point((int(co[0]),int(co[1])), fill=col)
        
    
    return chart_im
        

#Convert climate arrays to color images
#   maps: dictionary containing arrays of climate zones
def Make_image(maps, outname=None, in_opts = None):
    if in_opts:
        Save_opts(in_opts)
    if not outname:
        if opt('outname'):
            outname = opt('outname')
        else:
            outname = 'map'
    colmap = Get_colmap()
    print(' Outputting climate zone map...')
    for k, v in maps.items():
        outmap = np.empty((v.shape[0],v.shape[1],3), dtype=np.uint8)    # note dimensions switched from (lat,lon) to (x,y)
        for x in range(v.shape[1]):
            for y in range(v.shape[0]):     # maybe someday I'll properly implement a numpy indexing routine
                clim = v[y,x]
                outmap[y,x,:] = colmap[clim]
        outim = Image.fromarray(outmap)
        if opt('image_scale'):
            scale = opt('image_scale')
            if type(scale) is not tuple:
                scale = (round(outim.size[0]*scale),round(outim.size[1]*scale))
            verb(f'   Scaling output image to {scale}')
            outim = outim.resize(scale, resample=Image.Resampling.NEAREST)
        if len(maps) > 1:
            savename = outname + '_' + k + '.png'
        else:
            savename = outname + '.png'
        print(f'  Saving to {savename}')
        outim.save(savename)
    if opt('make_key'):
        print(f'  Saving map key to {outname}_key.png')
        key_im = Make_key(maps, colmap)
        key_im.save(outname + '_key.png')
    return outim

 

## Command functions

#Resets options, then adds option list or reads config file if provided
# in_opts: dictionary of options or name of config file
def Save_opts(in_opts=None):
    reset_default() #reset all options to default at top of script
    if os.path.exists(path+'kpasta_options.cfg'):
        print(" Loading config options from kpasta_options.cfg (does not override previously set options)")
        add_opt(Load_options(path+'kpasta_options.cfg'))    #load options from kpasta_options, overriding defaults
    if in_opts:
        if isinstance(in_opts, str):    #loads config file if specified
            if os.path.exists(in_opts):
                in_opts = Load_options(in_opts)
            elif os.path.exists(path+in_opts):
                in_opts = Load_options(path+in_opts)
            else:
                print(f" No config file found at {in_opts}")
        add_opt(in_opts)    #input list override defaults and kpasta_options

    for o in ('color_file','topo_map'):
        f = opt(o)

        if f is not None and not os.path.exists(f):
            if os.path.exists(path+f):  #checks both raw input and input with path for files
                add_opt({o: path+f})
            else:
                print(f'WARNING: no file found at {f} for {o}')

    f = opt('outname')
    if "/" not in f and "\\" not in f:  #if outname doesn't look like a path, presume path should be added
        add_opt({'outname': path+f})
    
    try:
        evap = opt('evap_estimate_sea')     #backwards compatibility with deprecated options
        if evap:
            add_opt({'estimate_evap': 'sea'})
        else:
            add_opt({'estimate_evap': 'never'})
    except:
        pass

    if opt('pas_simple_input'):     #meta option to set all other options
        verb('   Setting all options for simple pasta input per pas_simple_input')
        add_opt({
            'pet_method': 'kalike',
            'estimate_evap': 'all',
            'gdd_limit_light': False,
            'temp_tunings': 'tavg',
            'sea_use_ts': False,
            'sea_ice_use_temp': True,
            'pas_boil_pres': False,
            'pas_ice_def': 'tavg'
            })

    return

#Produce climate zone arrays from files
# returns dictionary of arrays with climates
#  files: list of input files
#  in_opts: dictionary of options or name of config file
def Make_clim(files, in_opts=None):
    if in_opts:
        Save_opts(in_opts)
    if opt('efficient'):
        try:
            land_funcs = Clim_func[opt('land_type')+'_efficient']
            sea_funcs = Clim_func[opt('sea_type')+'_efficient']
        except:
            print(f" Efficient function for {opt('land_type')} not found; using regular function")
            land_funcs = Clim_func[opt('land_type')]
            sea_funcs = Clim_func[opt('sea_type')]
            add_opt({'efficient': False})
    else:
        land_funcs = Clim_func[opt('land_type')]
        sea_funcs = Clim_func[opt('sea_type')]
    params = Get_params(files, land_funcs, sea_funcs)
    maps = Get_clims(params, land_funcs, sea_funcs)

            
    return maps

#Main routine: finds configs, runs Make_clim, and then produces output map
# files: name of file or list containing files
# in_opts: dictionary of options or name of config file
def Make_map(files, in_opts=None):
    if isinstance(files, str):      #makes sure files is a list
        files = File_search(files)
    Save_opts(in_opts)
    maps = Make_clim(files)
    Make_image(maps)
    return



















### CLIMATE ZONE ALGORITHMS ##############################################################################################################
    
## Templates

# For each climate classification, 3 functions should be used:
#   a _Data function to get data from given input files, convert to appropriate units, and save them to a dictionary
#   a _Param function to process that data to key parameters for determining zones, and save those to another dictionary
#   a _Alg function to determine climate zones from parameters for a single cell and return the climate zone
#       (though you can use existing _Data or _Param functions if they provide the necessary data)


#_Data function:
#   dat: raw climate data, e.g. a netcdf data object
def Template_Data(dat):
    data = Get_nc(dat, 'key')   # use Get_nc to extract and process data from files
    data *= 360                 # convert units as necessary
    all_data = dict(            # add to dictionary
        data=data               # use this format when creating dictionary to store each array by its own name
        )
    if True:
        data2 = Get_nc(dat, 'key2') # use optional data as appropriate and add to dictionary
        all_data['data2'] = data2   # when adding to existing dictionary, must name keys with strings (with quotation marks)
    return all_data
# all_data can contain 3d (time, lat, lon) or 2d (lat, lon) arrays

#_Param function:
#   all_data from _Data function
def Template_Param(all_data):
    data = all_data['data']     # take data from dictionary (note need for string key)
    if True:
        data2 = all_data['data2']   # take optional data as appropriate
    param_mean = np.mean(data, 0)   # produce key parameters however required
    param_max = np.max(data, 0)
    all_param = dict(              # add all to new dictionary
        param_mean=param_mean,
        param_max=param_max
        )
    return all_param
#all_param should contain only 2d (lat, lon) arrays


#_Alg function:
#   all_param from _Param function, but provided per-cell
def Template_Alg(all_param):
    clim = 0                        # set default clim for safety and debug
    param_mean = all_param['param_mean']    # take data from dictionary
    param_max = all_param['param_max']
    if param_mean > 5:
        if param_max > 10:      # determine climate zone as appropriate
            clim = A            # using if/else trees makes for slow execution but is easier to write and modify
        else:
            clim = B
    else:
        clim = C
    return clim
#clim should be a single climate zone id

# Functions should be added to the Clim_func dictionary

Clim_func['Template'] = (Template_Data, Template_Param, Template_Alg)

#Some oprtional slightly more advanced data functions:
def Template_Data_advanced(dat):
    lat, lon = coords_from_file(dat,'lat key','lon key')    # can retrieve lat, lon coordinates from file for use with interpolation
    coords = (lat,lon)                              # but if not provided, they will be automatically searched for in the file and generated if not found

    mask = get_mask(dat, 'land/sea mask key', convert=True)  # produce a land/sea mask (use convert=True to change 0.0-1.0 scale to true/false mask cut at 0.5 on the scale)
                                                                # used for blending land and sea maps
    
    temp, adjust = Get_nc_adjust(dat, 'temp key', 'elevation key', coords)   # to allow for temperature adjustment by topography, run with temperature and model elevation data
                                                                                # (will automatically find topo map file, or return zero adjustment if not using one)
    max_t = Get_nc(dat, 'max key', bin_ext = 1)     # use bin_ext to ensure proper behavior during binning (1 for max, -1 for min, 0 or leave out for average)

    all_data = dict(
        temp=temp,
        mask=mask,
        adjust=adjust   #include temperature adjustment so it can be reused by sea algorithm
        )
    return all_data

#Template for alternate efficient algorithm function, working on whole climate array at once:
def Template_Alg_efficient(all_param):
    for k, v in all_param.items():
        clim = np.zeros_like(v) #copy first available parameter for shape of climate array
        break
    
    param_mean = all_param['param_mean']    # take data from dictionary as normal
    param_max = all_param['param_max']

    clim = np.where(param_mean > 5,
                    np.where(param_max > 10, A,      # use numpy functions to work over whole parameter arrays at once
                        B),
                            C)
    return clim #return resulting array
#clim here is a 2d array of climate zone ids for the whole map

# Create new entry for efficient algorithm
# must have same name as normal entry with _efficient added, can use same _Data and _Param functions
Clim_func['Template_efficient'] = (Template_Data, Template_Param, Template_Alg_efficient)


## Koppen-Geiger

def Koppen_Data(dat):
    tas, adjust = Get_nc_adjust(dat, 'tas', 'grnz')

    pr = Get_nc(dat, 'pr')
    
    all_data = dict(
        tas=tas + opt('temp_adjust'),    #convert from K to C
        pr=pr * opt('precip_adjust'),  #convert from m/s to mm/month
        adjust=adjust
        )

    if opt('kg_summer_zenith'):
        czen = Get_nc(dat, 'czen')
        all_data['czen'] = czen
        
    return all_data

def Koppen_Param(data):
    tas = data['tas']   #2-meter air temp in C
    pr = data['pr']     #precipitaiton in mm/month

    verb('    Finding basic parameters')
    Avg_Temp = np.mean(tas, axis=0)
    Total_Precip = np.mean(pr, axis=0)*12 #Converts to mm/year
    Max_Temp = np.amax(tas, axis=0)
    Min_Temp = np.amin(tas, axis=0)
    Min_Precip = np.amin(pr, axis=0)

    timel = len(tas)        #Should come back and optimize this whole section at some point, it's such a memory hog
    if timel > 1:
        verb('    Dividing year into summer and winter halfs')
        halfl = math.floor(timel/2)
        if opt('kg_summer_zenith'):
            verb('     Using solar zenith angle')
            czen = data['czen']     #cosine of zenith star angle
            long = np.concatenate((czen,czen[:halfl,:,:]), axis=0)
        else:
            verb('    Using monthly temperature')
            long = np.concatenate((tas,tas[:halfl,:,:]), axis=0)
        sum_ar = np.sum(np.stack([long[m:m+halfl] for m in range(halfl)], 0), axis=0)   #stack half-year slices together and then sum them together to produce array of total values over following half-year for each month
        pr_long = np.concatenate((pr,pr),axis=0)    #create double-length array of precipitation
        pr_stack = np.stack([pr_long[m:m+timel] for m in range(timel)], axis=0)     #stack full-year slices together starting on each month of the year
        sum_max = np.argmax(sum_ar, axis=0)
        verb('    Finding summer and winter precipitation')
        precips = np.empty_like(pr)
        for y in range(pr.shape[1]):       #Bit of an inefficient approach, but haven't quite bothered to figured out how to do this as a numpy operation
            for x in range(pr.shape[2]):
                precips[:,y,x] = pr_stack[sum_max[y,x],:,y,x]    #choose appropriate slice corresponding to the maximum half-year total of indicator value
        Summer_Precip = np.mean(precips[:halfl,:,:], axis=0)*6    #convert to mm/half-year
        if timel%2 != 0:
            for y in range(pr.shape[1]):
                for x in range(pr.shape[2]):
                    ind = sum_max[y,x]          #for odd total months, round down for halfl then add in precip of hottest neighboring month to average at half weight
                    if tas[-1 if ind == 0 else ind-1,y,x] > tas[0 if ind+halfl == timel else ind+halfl+1,y,x]:
                        add_odd = pr[-1 if ind == 0 else ind-1,y,x]
                    else:
                        add_odd = pr[0 if ind+halfl == timel else ind+halfl+1,y,x]
                    Summer_Precip[y,x] = (halfl*2*Summer_Precip + 6*add_odd)/(halfl*2+1)
        verb('    Finding seasonal precipitation extremes')
        Max_Sum_Precip = np.amax(precips[:halfl,:,:], axis=0)
        Min_Sum_Precip = np.amin(precips[:halfl,:,:], axis=0)
        Max_Win_Precip = np.amax(precips[halfl:,:,:], axis=0)   #winter is counted as 1-month longer for odd months, but eh it's good enough for now
        Min_Win_Precip = np.amin(precips[halfl:,:,:], axis=0)
    else:
        verb('    Using backup calculations for "seasonless" equivalents of seasonal parameters')
        Summer_Precip = Total_Precip/2      #backup options for seasonless zones
        Max_Sum_Precip = Total_Precip/12
        Min_Sum_Precip = Max_Sum_Precip
        Max_Win_Precip = Max_Sum_Precip
        Min_Win_Precip = Max_Sum_Precip
    
    verb('    Determining length of "summer" above 10 C')
    Summer_Length = np.sum(np.where(tas>10,1,0), axis=0)/len(tas) #Portion of year above 10 C

    all_param = dict(
        Avg_Temp = Avg_Temp,
        Total_Precip = Total_Precip,
        Max_Temp = Max_Temp,
        Min_Temp = Min_Temp,
        Min_Precip = Min_Precip,
        Summer_Precip = Summer_Precip,
        Max_Sum_Precip = Max_Sum_Precip,
        Min_Sum_Precip = Min_Sum_Precip,
        Max_Win_Precip = Max_Win_Precip,
        Min_Win_Precip = Min_Win_Precip,
        Summer_Length = Summer_Length
        )
    
    return all_param

def Koppen_Alg(par):
    clim=0
    Avg_Temp = par['Avg_Temp']
    Total_Precip = par['Total_Precip']
    Max_Temp = par['Max_Temp']
    Min_Temp = par['Min_Temp']
    Summer_Precip = par['Summer_Precip']

    #determine threshold for arid zone
    if opt('kg_trewartha_arid'):
        Arid_threshold = Avg_Temp*23 - 640 * (Total_Precip - Summer_Precip) / max(Total_Precip,0.001) + 410
    else:
        if Summer_Precip > Total_Precip*0.7:
            adjust = 280
        elif Summer_Precip > Total_Precip*0.3:
            adjust = 140
        else:
            adjust = 0
        Arid_threshold = Avg_Temp*20 + adjust

    arpol = opt('kg_arid_polar_priority')
    cold = opt('kg_temperate_min')
    
    #Groups
        
    if Total_Precip < Arid_threshold and (Max_Temp > 10 or arpol == 'arid' or (Max_Temp > 0 and arpol == 'EF')):
        clim = B
    elif Max_Temp < 10:
        clim = E
    elif Min_Temp > 18:
        clim = A
    elif Min_Temp > cold:
        clim = C
    else:
        clim = D

    if opt('land_subtype') == 'groups': #finish here for groups only
        return clim
                
    #Full Koppen set
        
    else:
        Min_Precip = par['Min_Precip']
        Max_Sum_Precip = par['Max_Sum_Precip']
        Min_Sum_Precip = par['Min_Sum_Precip']
        Max_Win_Precip = par['Max_Win_Precip']
        Min_Win_Precip = par['Min_Win_Precip']
        Summer_Length = par['Summer_Length']
        
        if clim == B:
            if Total_Precip < Arid_threshold/2:    #desert/steppe test
                clim = BW
            else:
                clim = BS
            if opt('land_subtype') != 'two_letter':
                Bxh = ((Min_Temp > cold and not opt('kg_arid_avg'))     #hot/cold test
                       or (Avg_Temp > 18 and opt('kg_arid_avg')))    
                if clim == BW:
                    if Bxh:
                        clim = BWh
                    else:
                        clim = BWk
                else:
                    if Bxh:
                        clim = BSh
                    else:
                        clim = BSk
        elif clim == E:
            if Max_Temp < 0:
                clim = EF
            else:
                clim = ET
        elif clim == A:
            if Min_Precip > 60:
                clim = Af
            elif Min_Precip > 100-Total_Precip/25:
                clim = Am
            else:
                Xs = False
                if opt('land_subtype') == 'full':
                    if opt('kg_med_as'):
                        if (Min_Sum_Precip < opt('kg_med_summer_precip')
                            and Max_Sum_Precip > 3 * Min_Sum_Precip
                            and (Summer_Precip < Total_Precip/2
                                or (opt('kg_wet_season_req') != 'med_strict'
                                    and opt('kg_wet_season_req') != 'all_strict'))):
                            Xs = True
                    elif Summer_Precip < Total_Precip/2:
                        Xs = True
                if Xs:
                    clim = As
                else:
                    clim = Aw
        else:
            if opt('kg_trewartha_seasons'):
                Xs = (Summer_Precip < Total_Precip/4
                      and Min_Precip < 30
                      and Total_Precip < 890)
                Xw = Summer_Precip > Total_Precip * 0.7 if opt('kg_wet_summer_total') else 10/11
            else:
                Xs = (Min_Sum_Precip < opt('kg_med_summer_precip')  #med test
                      and Max_Sum_Precip > 3 * Min_Sum_Precip
                      and (Summer_Precip < Total_Precip/2
                           or (opt('kg_wet_season_req') != 'med_strict'
                               and opt('kg_wet_season_req') != 'all_strict')))
                Xw = ((Max_Sum_Precip > 10 * Min_Win_Precip     #wet-summer test
                       and not opt('kg_wet_summer_total')
                       and (Summer_Precip > Total_Precip/2
                            or (opt('kg_wet_season_req') != 'sum_strict'
                                and opt('kg_wet_season_req') != 'all_strict')))
                      or (Summer_Precip > 0.7 * Total_Precip
                          and opt('kg_wet_summer_total')))
            if Xs and Xw:
                if opt('kg_wet_summer_priority'):   #determine med or wet-summer where they overlap
                    Xs = False
                else:
                    Xw = False
            if Summer_Length < 1/3:    #ab/cd test
                if Min_Temp > -38:  #c/d test
                    tlet = 3
                else:
                    tlet = 4
            else:
                if Max_Temp > 22:   #a/b test
                    tlet = 1
                else:
                    tlet = 2
            if opt('land_subtype') != 'reduced':    #only bother with this if not doing reduced sets
                if opt('land_subtype') == 'two_letter':
                    if clim == C:
                        if Xs:
                            clim = Cs
                        elif Xw:
                            clim = Cw
                        else:
                            clim = Cf
                    else:
                        if Xs:
                            clim = Ds
                        elif Xw:
                            clim = Dw
                        else:
                            clim = Df
                elif clim == C:       #choose within C and D by overlap of subtypes
                    if Xs:
                        if tlet == 1:
                            clim = Csa
                        elif tlet == 2:
                            clim = Csb
                        else:
                            clim = Csc
                    elif Xw:
                        if tlet == 1:
                            clim = Cwa
                        elif tlet == 2:
                            clim = Cwb
                        else:
                            clim = Cwc
                    else:
                        if tlet == 1:
                            clim = Cfa
                        elif tlet == 2:
                            clim = Cfb
                        else:
                            clim = Cfc
                else:
                    if Xs:
                        if tlet == 1:
                            clim = Dsa
                        elif tlet == 2:
                            clim = Dsb
                        elif tlet == 3:
                            clim = Dsc
                        else:
                            clim = Dsd
                    elif Xw:
                        if tlet == 1:
                            clim = Dwa
                        elif tlet == 2:
                            clim = Dwb
                        elif tlet == 3:
                            clim = Dwc
                        else:
                            clim = Dwd
                    else:
                        if tlet == 1:
                            clim = Dfa
                        elif tlet == 2:
                            clim = Dfb
                        elif tlet == 3:
                            clim = Dfc
                        else:
                            clim = Dfd
                            
        if opt('land_subtype') == 'reduced':
            if clim == Af:      #replace climates with reduced equivalents
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
            elif clim == C:      #C and D use more simple reduced breakdown
                if Xs:
                    clim = Med
                elif tlet == 1:
                    clim = Subtropical
                else:
                    clim = Oceanic
            elif clim == D:
                if tlet < 3:
                    clim = Continental
                else:
                    clim = Subarctic
            elif clim == ET:
                clim = Tundra
            elif clim == EF:
                clim = IceCap
                
    return clim                    

Clim_func['Koppen-Geiger'] = (Koppen_Data, Koppen_Param, Koppen_Alg)

## Trewartha

def Trewartha_Alg(par):
    clim=0
    Avg_Temp = par['Avg_Temp']
    Total_Precip = par['Total_Precip']
    Max_Temp = par['Max_Temp']
    Min_Temp = par['Min_Temp']
    Summer_Precip = par['Summer_Precip']
    Summer_Length = par['Summer_Length']

    arpol = opt('kg_arid_polar_priority')

    #determine threshold for arid zones
    Arid_threshold = Avg_Temp*23 - 640 * (Total_Precip - Summer_Precip) / max(Total_Precip,0.001) + 410
    
    #Groups
        
    if Total_Precip < Arid_threshold and (Max_Temp > 10 or arpol == 'arid' or (Max_Temp > 0 and arpol == 'ET')):
        clim = TrB
    elif Max_Temp < 10:
        clim = TrF
    elif Min_Temp > 18:
        clim = TrA
    elif Summer_Length < 1/3:
        clim = TrE
    elif Summer_Length < 2/3:
        clim = TrD
    else:
        clim = TrC

    if opt('land_subtype') == 'groups': #finish here for groups only
        return clim
                
    #Full set
        
    else:
        Min_Precip = par['Min_Precip']
        
        if clim == TrB:
            if Total_Precip < Arid_threshold/2:
                clim = TrBW
            else:
                clim = TrBS
                
        elif clim == TrF:
            if Max_Temp < 0:
                clim = TrFi
            else:
                clim = TrFt
                
        elif clim == TrA:
            if Min_Precip > 60:
                clim = TrAr
            else:
                Xs = False
                if opt('land_subtype') == 'full':
                    if opt('kg_med_as'):
                        if (Summer_Precip < Total_Precip/4
                            and Min_Precip < opt('kg_med_summer_precip')
                            and Total_Precip < 890):
                            Xs = True
                    elif Summer_Precip < Total_Precip/2:
                        Xs = True
                if Xs:
                    clim = TrAs
                else:
                    clim = TrAw
        elif clim == TrC:
            if (Summer_Precip < Total_Precip/4
                and Min_Precip < opt('kg_med_summer_precip')
                and Total_Precip < 890):
                clim = TrCs
            elif Summer_Precip > Total_Precip * (0.7 if opt('kg_wet_summer_total') else 10/11):
                clim = TrCw
            else:
                clim = TrCf
        else:
            Xo = Min_Temp > opt('kg_temperate_min')
            if clim == TrD:
                if Xo:
                    clim = TrDo
                else:
                    clim = TrDc
            elif opt('land_subtype') == 'full':
                if Xo:
                    clim = TrEo
                else:
                    clim = TrEc
                
    return clim                    

Clim_func['Trewartha'] = (Koppen_Data, Koppen_Param, Trewartha_Alg) #Shares data and parameters with koppen

## Holdridge Life Zones

def Holdridge_Data(dat):
    tas, adjust = Get_nc_adjust(dat, 'tas', 'grnz')

    pr = Get_nc(dat, 'pr')
    
    all_data = dict(
        tas=tas + opt('temp_adjust'),    #convert from K to C
        pr=pr * opt('precip_adjust'),  #convert from m/s to mm/month
        adjust=adjust
        )
    
    if not opt('h_no_pet'):
        pet = Get_pet(dat, all_data)
        all_data['pet'] = pet
        
    return all_data

def Holdridge_Param(data):
    tas = data['tas']   #2-meter air temperature in C
    pr = data['pr']     #precipitation in mm/month

    verb('    Calculating total precipitation')
    Total_Precip = np.mean(pr, 0)*12

    if opt('h_estimate_biot_avg'):  #rough attempt to estimate average biotemperature based only on a single annual average temperature
        verb('    Estimating biotemperature from average temperature')
        Ta = tas - 15
        Tb = tas - 25
        biot = np.where(tas < -20, 0,
                    np.where(tas < 15, tas + 0.013*Ta**2 - 0.12*Ta,
                        np.where(tas < 25, tas,
                            np.where(tas < 35, tas - 0.046*Tb**2 - 0.02*Tb,
                                30))))
    else:
        verb('    Calculating biotemperature')
        biot = np.where(tas < 0, 0,
                    np.where(tas > 30, 30,
                             tas))
    Avg_Biot = np.mean(biot, 0)
    
    all_param = dict(
        Total_Precip = Total_Precip,
        Avg_Biot = Avg_Biot
        )

    if not opt('h_no_pet'):
        verb('    Calculating PETR')
        pet = data['pet']   #potential evapotranspiration in mm/month
        PETR = 12*np.mean(pet, 0) / np.maximum(Total_Precip, 0.001)    #avoid div0 error
        all_param['PETR'] = PETR
    
    return all_param

def Holdridge_Alg(par):
    clim = 0
    Avg_Biot = par['Avg_Biot']
    Total_Precip = par['Total_Precip']
    
    if opt('h_no_pet'):
        
        if Avg_Biot > 24:       #simplifiied square-grid indexing by biotemp and precip
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
            
    else:
        PETR = par['PETR']
        pr_ind = math.log(max(Total_Precip,0.001)/62.5, 2)
        pe_ind = math.log(max(PETR,0.001), 2) + 3
        ind = pr_ind - pe_ind # indicator of how far left or right the intersection of precip and PETR is on the ternary diagram
        if Avg_Biot > 24:       #staggered-grid indexing by biotemp and precip/PETR indicator
            if ind > 6:
                clim = HTropRainForest
            elif ind > 4:
                clim = HTropWetForest
            elif ind > 2:
                clim = HTropMoistForest
            elif ind > 0:
                clim = HTropDryForest
            elif ind > -2:
                clim = HVDryForest
            elif ind > -4:
                clim = HThornWood
            elif ind > -6:
                clim = HTropDesertScrub
            else:
                clim = HTropDesert
        elif Avg_Biot > 12:
            if ind > 5:
                clim = HWarmRainForest
            elif ind > 3:
                clim = HWarmWetForest
            elif ind > 1:
                clim = HWarmMoistForest
            elif ind > -1:
                clim = HWarmDryForest
            elif ind > -3:
                clim = HThornSteppe
            elif ind > -5:
                clim = HWarmDesertScrub
            else:
                clim = HWarmDesert
        elif Avg_Biot > 6:
            if ind > 4:
                clim = HCoolRainForest
            elif ind > 2:
                clim = HCoolWetForest
            elif ind > 0:
                clim = HCoolMoistForest
            elif ind > -2:
                clim = HSteppe
            elif ind > -4:
                clim = HCoolDesertScrub
            else:
                clim = HCoolDesert
        elif Avg_Biot > 3:
            if ind > 3:
                clim = HBorRainForest
            elif ind > 1:
                clim = HBorWetForest
            elif ind > -1:
                clim = HBorMoistForest
            elif ind > -3:
                clim = HDryScrub
            else:
                clim = HBorDesert
        elif Avg_Biot > 1.5:
            if ind > 2:
                clim = HRainTundra
            elif ind > 0:
                clim = HWetTundra
            elif ind > -2:
                clim = HMoistTundra
            else:
                clim = HDryTundra
        else:
            clim = HPolDesert
        
    return clim

Clim_func['Holdridge'] = (Holdridge_Data, Holdridge_Param, Holdridge_Alg)


##Thornthwaite-Feddema

def Thornthwaite_Data(dat):
    tas, adjust = Get_nc_adjust(dat, 'tas', 'grnz')
    tas += opt('temp_adjust') #convert from K to C

    pr = Get_nc(dat, 'pr')

    pet = Get_pet(dat,{'tas': tas})
    
    all_data = dict(
        pr=pr * opt('precip_adjust'),  #convert from m/s to mm/month
        tas = tas,  #include for sea zones for chart
        pet = pet,
        adjust=adjust
        )
    
    if opt('interp_scale') and opt('topo_map'):
        all_data['adjust'] = adjust
    

        
    return all_data

def Thornthwaite_Param(data):
    pet = data['pet']   #potential evapotranspiration in mm/month
    pr = data['pr']     #precipitation in mm/month



    if opt('land_subtype') == 'full':
        pet = np.maximum(np.mean(pet, 0), 0.001)    #prevent div0 errors
        pr_m = np.maximum(np.mean(pr, 0), 0.001)
        
        indm = np.where(pr_m > pet, 1-pet/pr_m, pr_m/pet-1)   #moisture indicator
        indt = pet*12/300   #thermal indicator

        all_param = dict(
            indm=(indm + 1)*3,  #convert to 0-6 scale
            indt=indt
            )
    else:
        pet = np.maximum(pet, 0.001)    #prevent div0 errors
        pr = np.maximum(pr, 0.001)
        indm = np.where(pr > pet, 1-pet/pr, pr/pet-1)

        indm_r = np.ptp(indm,0)  #range of variability indicator

        pet_r = np.ptp(pet,0)
        pr_r = np.ptp(pr,0)
        type_r = pet_r / pr_r   #type of variability indicato
        

        all_param = dict(
            indm_r = indm_r,   #convert to 0-4 scale
            type_r = type_r
            )

    return all_param

def Thornthwaite_Alg(par):
    clim=0
    if opt('land_subtype') == 'full':
        indm = par['indm']
        indt = par['indt']

        if indt > 5:        #simple square-grid indexing by indicators
            if indm > 5:
                clim = THTorSat
            elif indm > 4:
                clim = THTorWet
            elif indm > 3:
                clim = THTorMoist
            elif indm > 2:
                clim = THTorDry
            elif indm > 1:
                clim = THTorSemi
            else:
                clim = THTorArid
        elif indt > 4:
            if indm > 5:
                clim = THHotSat
            elif indm > 4:
                clim = THHotWet
            elif indm > 3:
                clim = THHotMoist
            elif indm > 2:
                clim = THHotDry
            elif indm > 1:
                clim = THHotSemi
            else:
                clim = THHotArid
        elif indt > 3:
            if indm > 5:
                clim = THWarmSat
            elif indm > 4:
                clim = THWarmWet
            elif indm > 3:
                clim = THWarmMoist
            elif indm > 2:
                clim = THWarmDry
            elif indm > 1:
                clim = THWarmSemi
            else:
                clim = THWarmArid
        elif indt > 2:
            if indm > 5:
                clim = THCoolSat
            elif indm > 4:
                clim = THCoolWet
            elif indm > 3:
                clim = THCoolMoist
            elif indm > 2:
                clim = THCoolDry
            elif indm > 1:
                clim = THCoolSemi
            else:
                clim = THCoolArid
        elif indt > 1:
            if indm > 5:
                clim = THColdSat
            elif indm > 4:
                clim = THColdWet
            elif indm > 3:
                clim = THColdMoist
            elif indm > 2:
                clim = THColdDry
            elif indm > 1:
                clim = THColdSemi
            else:
                clim = THColdArid
        else:
            if indm > 5:
                clim = THFrigSat
            elif indm > 4:
                clim = THFrigWet
            elif indm > 3:
                clim = THFrigMoist
            elif indm > 2:
                clim = THFrigDry
            elif indm > 1:
                clim = THFrigSemi
            else:
                clim = THFrigArid
                
    else:
        indm_r = par['indm_r']
        type_r = par['type_r']
        if type_r > 2:
            if indm_r > 1.5:
                clim = THTempExt
            elif indm_r > 1:
                clim = THTempHigh
            elif indm_r > 0.5:
                clim = THTempMed
            else:
                clim = THTempLow
        elif type_r > 0.5:
            if indm_r > 1.5:
                clim = THComboExt
            elif indm_r > 1:
                clim = THComboHigh
            elif indm_r > 0.5:
                clim = THComboMed
            else:
                clim = THComboLow
        else:
            if indm_r > 1.5:
                clim = THPrecipExt
            elif indm_r > 1:
                clim = THPrecipHigh
            elif indm_r > 0.5:
                clim = THPrecipMed
            else:
                clim = THPrecipLow
                
    return clim

Clim_func['Thornthwaite'] = (Thornthwaite_Data, Thornthwaite_Param, Thornthwaite_Alg)

##Whittaker Biomes

def Whittaker_Data(dat):
    tas, adjust = Get_nc_adjust(dat, 'tas', 'grnz')

    pr = Get_nc(dat, 'pr')

    all_data = dict(
        tas=tas + opt('temp_adjust'),    #convert from K to C
        pr=pr*opt('precip_adjust'),     #convert from m/s to mm/month
        adjust=adjust
        )

    return all_data

def Whittaker_Param(data):
    tas = data['tas']   #2-meter air temperature in C
    pr = data['pr']     #precipitaiton in mm/month

    Avg_Temp = np.mean(tas,0)
    Total_Precip = np.mean(pr,0)*12    #sum and convert to mm/year

    all_param = dict(
        Avg_Temp=Avg_Temp,
        Total_Precip=Total_Precip
        )

    return all_param

def Whittaker_Alg(par):
    clim=0
    ta = par['Avg_Temp']
    pr = par['Total_Precip'] / 10   #convert to cm/year

    if ta > -5 and (pr > 300 or pr > 0.016*(ta**3) - 0.7938*(ta**2) + 14.736*ta + 129.67):
        hum = 4
    elif pr > 150 or pr > -0.0025*(ta**3) + 0.0495*(ta**2) + 4.5373*ta + 39.952:
        hum = 3
    elif pr > 100 or pr > 0.0027*(ta**3) - 0.054*(ta**2) + 1.6589*ta + 23.452:
        hum = 2
    else:
        hum = 1

    if ta > 20 or ta > (6.71e-7)*(pr**3) - 0.000399*(pr**2) + 0.0761*pr + 14.068:
        if hum == 4:
            clim = WTropRain
        elif hum > 1:
            clim = WTropSeas
        else:
            clim = WSubtropDes
    elif ta > 5 or ta > 1.51e-6*(pr**3) - 0.000624*(pr**2) + 0.0942*pr - 2.248:
        if hum == 4:
            clim = WTempRain
        elif hum == 3:
            clim = WTempSeas
        elif hum == 2:
            clim = WWood
        else:
            clim = WTempDes
    elif ta > -4 or ta > -0.000114*(pr**2) + 0.0403*pr - 7.24:
        if hum > 2:
            clim = WBor
        elif hum == 2:
            clim = WWood
        else:
            clim = WTempDes
    else:
        if hum > 2 or ta < -10:
            clim = WTun
        elif hum == 2:
            clim = WWood
        else:
            clim = WTempDes

    return clim

Clim_func['Whittaker'] = (Whittaker_Data, Whittaker_Param, Whittaker_Alg)

##Two-Parameter Koppen-Geiger

def TwoParamKG_Alg(par):
    clim=0
    ta = par['Avg_Temp']
    pr = par['Total_Precip'] / 10   #convert to cm/year

    #determine threshold for arid zone
    
    Arid_threshold = 2*ta + 14
    if ta < 5:
        Arid_threshold = min(Arid_threshold + 14, 24)

    #Groups
        
    if pr < Arid_threshold:
        clim = B
    else:
        if pr < 28:
            if ta < min(-10, -0.2*pr - 8.33):
                clim = E
        elif ta < min(1, 0.144*pr - 18.1):
            clim = E
    if clim == 0:
        if ta < max(6, -0.0245*pr + 11.1):
            clim = D
        elif ta < max(20, -0.0169*pr + 24):
            clim = C
        else:
            clim = A

    if opt('land_subtype') == 'groups': #finish here for groups only
        return clim
                
    #Full set
        
    else:
        if clim == B:
            if pr < Arid_threshold/2:    #desert/steppe test
                clim = BW
            else:
                clim = BS
            Bxh = ta > 18     #hot/cold test   
            if clim == BW:
                if Bxh:
                    clim = HotDesert
                else:
                    clim = ColdDesert
            else:
                if Bxh:
                    clim = HotSteppe
                else:
                    clim = ColdSteppe
        elif clim == E:
            clim = Tundra
            if pr < 72:
                if ta < max(-20, 0.14*pr - 20.8):
                    clim = IceCap
            elif ta < min(-5, 0.0333*pr - 13.1):
                clim = IceCap
        elif clim == D:
            if ta < max(1, 0.0409*pr - 2.82):
                clim = Subarctic
            else:
                clim = Continental
        elif clim == C:
            if ta < 14:
                clim = Oceanic
            else:
                clim = Subtropical
        else:
            if pr > 200:
                clim = TropRainforest
            else:
                clim = TropSavanna
                
    return clim

Clim_func['TwoParamKG'] = (Whittaker_Data, Whittaker_Param, TwoParamKG_Alg)
    

##Woodward Vegetation Types

def Woodward_Data(dat):

    tas, adjust = Get_nc_adjust(dat, 'tas', 'grnz')

    if opt('temp_adjust_ts'):
        ts = Get_nc(dat, 'ts', adjust=adjust)
        tdif = tas-ts
    else:
        tdif = 0
    mint = Get_nc(dat, 'mint', adjust=adjust+tdif, bin_ext = -1)


    pr = Get_nc(dat, 'pr')
    
    all_data = dict(
        tas=tas+opt('temp_adjust'),
        mint=mint+opt('temp_adjust'),
        pr=pr*opt('precip_adjust'),    #convert from m/s to mm/month
        adjust=adjust
        )
    
    pet = Get_pet(dat,all_data)

    all_data['pet'] = pet

    return all_data

def Woodward_Param(data):
    tas = data['tas']   #2-meter air temperature in C
    pr = data['pr']     #precipitaiton in mm/month
    pet = data['pet']   #PET in mm/month
    mint = data['mint']

    Avg_Temp = np.mean(tas,0)   #include for charts and WCR
    Abs_Min = np.amax(mint,0)
    
    Total_Precip = np.mean(pr,0)*12    #sum and convert to mm/year
    Arid_f = np.mean(pr,0) / np.maximum(np.mean(pet,0),0.001) #ratio of total precip to total PET
    Deg_Month = Calc_GDD(tas, base=0)/30    #growing degree months
    Max_Temp = np.amax(tas,0)   #include for IPCC and WCR

    all_param = dict(
        Avg_Temp=Avg_Temp,
        Abs_Min=Abs_Min,
        Total_Precip=Total_Precip,
        Arid_f=Arid_f,
        Deg_Month=Calc_GDD_total(Deg_Month),
        Max_Temp=Max_Temp
        )

    return all_param

def Woodward_Alg(par):
    clim=0
    Abs_Min = par['Abs_Min']
    Arid_f = par['Arid_f']
    Total_Precip = par['Total_Precip']
    Deg_Month = par['Deg_Month']

    th_con = -32
    th_decid = 0
    th_ever = 10
    th_chill = 15

    if opt('temp_tunings') == 'tclim':
        th_con = -25
        th_decid = 2
        th_ever = 9
        th_chill = 15

    elif opt('temp_tunings') == 'tavg':
        th_con = -19
        th_decid = 5
        th_ever = 15.5
        th_chill = 22

    if Deg_Month < 50:
        clim = WoTundra
    elif Abs_Min < th_con:
        if Total_Precip < 400:
            clim = WoDry
        else:
            clim = WoConifer
    elif Total_Precip < 600:
        clim = WoDry
    elif Abs_Min < th_decid:
        clim = WoDecidWin
    elif Abs_Min < th_ever:
        clim = WoEverFrost
    else:
        if Arid_f > 1:
            if Abs_Min > th_chill:
                clim = WoEverTrop
            else:
                clim = WoEverChill
        else:
            clim = WoDecidDry
    return clim

Clim_func['Woodward'] = (Woodward_Data, Woodward_Param, Woodward_Alg)

def IPCC_Alg(par):
    Avg_Temp = par['Avg_Temp']
    Max_Temp = par['Max_Temp']
    Arid_f = par['Arid_f']

    if Avg_Temp > 18:
        Total_Precip = par['Total_Precip']
        if Total_Precip > 2000:
            clim = ITropWet
        elif Total_Precip > 1000:
            clim = ITropMoist
        else:
            clim = ITropDry
    elif Avg_Temp > 10:
        if Arid_f > 1:
            clim = IWarmMoist
        else:
            clim = IWarmDry
    elif Avg_Temp > 0:
        if Arid_f > 1:
            clim = ICoolMoist
        else:
            clim = ICoolDry
    elif Max_Temp > 10:
        if Arid_f > 1:
            clim = IBorMoist
        else:
            clim = IBorDry
    else:
        if Arid_f > 1:
            clim = IPolMoist
        else:
            clim = IPolDry
    return clim

Clim_func['IPCC'] = (Woodward_Data, Woodward_Param, IPCC_Alg)

def WCR_Alg(par):
    Avg_Temp = par['Avg_Temp']
    Max_Temp = par['Max_Temp']
    Arid_f = par['Arid_f']

    if Avg_Temp > 24:
        if Arid_f > 0.65:
            clim = WCRTropMoist
        elif Arid_f > 0.05:
            clim = WCRTropDry
        else:
            clim = WCRTropDes
    elif Avg_Temp > 18:
        if Arid_f > 0.65:
            clim = WCRSubTropMoist
        elif Arid_f > 0.05:
            clim = WCRSubTropDry
        else:
            clim = WCRSubTropDes
    elif Avg_Temp > 10:
        if Arid_f > 0.65:
            clim = WCRWarmTempMoist
        elif Arid_f > 0.05:
            clim = WCRWarmTempDry
        else:
            clim = WCRWarmTempDes
    elif Avg_Temp > 0:
        if Arid_f > 0.65:
            clim = WCRCoolTempMoist
        elif Arid_f > 0.05:
            clim = WCRCoolTempDry
        else:
            clim = WCRCoolTempDes
    elif Max_Temp > 10:
        if Arid_f > 0.65:
            clim = WCRBorMoist
        elif Arid_f > 0.05:
            clim = WCRBorDry
        else:
            clim = WCRBorDes
    else:
        if Arid_f > 0.65:
            clim = WCRPolMoist
        elif Arid_f > 0.05:
            clim = WCRPolDry
        else:
            clim = WCRPolDes

    return clim

Clim_func['WCR'] = (Woodward_Data, Woodward_Param, WCR_Alg)

## Common data function for Prentice, Pasta, and unproxied KG

def Biome_Data(dat):

    tas, adjust = Get_nc_adjust(dat, 'tas', 'grnz', t_unadjust=True)    #get unadjusted tas at first

    rss = Get_nc(dat, 'rss')
    ps = Get_nc(dat, 'ps')

    if opt('temp_tunings') == 'tavg':
        verb('    Skipping maxt and mint extracting, just using average temp')
        maxt = - opt('temp_adjust')
        mint = - opt('temp_adjust')
        tdif = 0
    else:
        if opt('temp_adjust_ts'):
            verb('    Finding adjustment from surface temp to 2-meter air temp')
            ts = Get_nc(dat, 'ts')
            tdif = tas - ts #adjustment from surface to 2-meter to apply to maxt and mint
        else:
            tdif = 0
        maxt = Get_nc(dat, 'maxt', adjust=tdif, bin_ext = 1)
        mint = Get_nc(dat, 'mint', adjust=tdif, bin_ext = -1)
    pr = Get_nc(dat, 'pr', no_interp = True)

    if opt('estimate_evap') == 'all' or (opt('interp_scale') and opt('estimate_evap') == 'sea'):
        verb('    Gathering data to estimate evapotranspiration')
        pet_small = Get_pet(dat, {}, no_interp = True)  #use original resolution
        evap = Estimate_evap(pet_small, pr*opt('precip_adjust'))  #estimate evaporation from pet and pr
        if opt('estimate_evap') == 'sea':
            verb('    Combining evap data for land with estimated evap for sea')
            evap_land = Get_nc(dat, 'evap', no_interp = True)
            evap_land *= -opt('precip_adjust')   #convert from m/s to mm/month (and flip sign)
            evap = np.where(get_mask(dat, 'lsm', convert=True), evap_land, evap)    #apply evap estimation only to sea areas
        if opt('interp_scale'):
            evap = Interp(evap)
    else:
        evap = Get_nc(dat, 'evap')
        evap *= -opt('precip_adjust')   #convert from m/s to mm/month (and flip sign)
    
    pet= Get_pet(dat,
                    {
                        'tas':tas + opt('temp_adjust'),
                        'rss':rss,
                        'ps':ps,
                        'maxt':maxt + opt('temp_adjust'),
                        'mint':mint + opt('temp_adjust')
                    })
    
    pet = np.maximum(pet,evap)  #ensure pet is always >= evaporation

    tas += adjust   #apply temperature adjustments at end, so they don't affect PET
    maxt += adjust  #because we can't apply any adjustments to evap, so adjusting only PET would have weird effects on Ar
    mint += adjust

    if opt('interp_scale'):
        pr = Interp(pr)
    
    ssru = Get_nc(dat, 'ssru')
    rin = rss - ssru
    
    all_data = dict(
        tas=tas + opt('temp_adjust'),    #convert from K to C
        mint=mint + opt('temp_adjust'),  #convert from K to C
        maxt=maxt + opt('temp_adjust'),  #convert from K to C
        evap=evap,
        pet=pet,
        pr=pr * opt('precip_adjust'),
        ps=ps,
        rin=rin,
        adjust=adjust,
        tdif=tdif,
        )
    
    if opt('land_type') != 'Prentice':
        if opt('pas_ice_def') in ('ice', 'ice_noadj'):
            snd = Get_nc(dat, 'snd', no_interp = True, bin_ext = -1)
            ice = np.minimum(snd,0.4)   #helps get smoother interpolation
            if opt('interp_scale'):
                ice = Interp(ice)
                if opt('pas_ice_def') == 'ice':
                    all_data['tice'] = maxt - tdif + opt('temp_adjust') #unadjusted max surface temperature for ice adjustments with interpolation
            all_data['ice'] = ice
        elif opt('pas_ice_def') == 'maxt':
            all_data['tice'] = maxt - tdif + opt('temp_adjust')


    return all_data

def Biome_Param(data):
    tas = data['tas']   #2-meter air temperature in C
    evap = data['evap']     #evaporation in mm/month
    pr = data['pr']         #precipitation in mm/month
    pet = data['pet']   #potential evapotranspiration in mm/year

    if opt('temp_tunings') == 'tavg':
        maxt = tas
        mint = tas
    else:
        maxt = data['maxt']     #absolute maximum temperature in C
        mint = data['mint']     #absolute minimum temperature in C

    verb('    Finding temperature extremes')
    Min_Abs = np.amin(mint, 0)
    Max_Abs = np.amax(maxt, 0)
    Max_Avg = np.amax(tas, 0)
    Min_Avg = np.amin(tas, 0)

    Min_Abs = np.minimum(Min_Abs, Min_Avg)  #ensure absolute extremes are always beyond average extremes
    Max_Abs = np.maximum(Max_Abs, Max_Avg)

    ##GDD config

    GInt_th = opt('pas_gint_thresh')  #threshold of growth interruption to interrupt GDD accumulation; 1250 by default

    GDDbase = 5     #base temp for GDD (C)
    GDDplats = 25   #'plateau' start; point of maximum GDD
    GDDplate = 40   #'plateau' end; GDD starts declining again
    GDDtop = 50     #maximum growth temp

    GDDzbase = 0    #values for GDDz
    GDDzplats = 20
    GDDzplate = 40
    GDDztop = 60

    GDDlbase = 10 / opt('gdd_par_ratio')   #values for light-limited GDD maximum (W/m^2)
    GDDlplats = 110 / opt('gdd_par_ratio')  #corresponds to 20 and 220 for sunlike PAR ratio of 0.5

    GDDlzbase = 0   #values for light-limited GDDz maximum
    GDDlzplats = 100 / opt('gdd_par_ratio') #corresponds to 200 for sunlike PAR ratio of 0.5
        #note that these are separately defined for sea zones as a backup if GDDlz not produced here (but to the same default values)

    verb('    Calculating GDD parameters')
    GDD = Calc_GDD(tas, base=GDDbase, plat_start=GDDplats, plat_end=GDDplate, comp=GDDtop)
    GDDz = Calc_GDD(tas, base=GDDzbase, plat_start=GDDzplats, plat_end=GDDzplate, comp=GDDztop)

    GDDlz = None
    if opt('land_subtype') == 'earthlike' or opt('land_subtype') == 'earthlike_no_pluv':
        pass
    elif opt('gdd_limit_light'):
        rin = data['rin']
        GDDl = Calc_GDD(rin, base=GDDlbase, plat_start=GDDlplats) / 10
        GDDlz = Calc_GDD(rin, base=GDDlzbase, plat_start=GDDlzplats) / 10
        GDD = np.minimum(GDD, GDDl)
        GDDz = np.minimum(GDDz, GDDlz)
        GDDlz = Calc_GDD_total(GDDlz)   #save for use with sea zones

    GDDp = (GDD)/np.maximum(np.sum(GDD,0),0.001)    #portion of total GDD in each month

    verb('    Calculating moisture parameters')

    Ar = np.mean(evap,0) / np.maximum(np.mean(pet,0), 0.001)  #avoid div0 errors
    GAr = np.sum(GDDp*evap,0) / np.maximum(np.sum(GDDp*pet,0), 0.001)

    Evr = np.mean(evap,0) / np.maximum(np.mean(pr,0), 0.001)
    GrS = np.sum(GDDp*pr,0) / np.maximum(np.mean(evap, 0), 0.001)


    verb('    Totaling GDD')
    if opt('land_type') == 'Prentice':
        GDD = Calc_GDD_total(GDD)
        GInt = None
    else:
        GInt = opt('month_length') * 450-GDDz
        GInt = np.maximum(GInt, 0)      #clamp to positive
        GDD, GInt = Calc_GDD_total(GDD, GInt, th_gi=GInt_th)

    GDDz = Calc_GDD_total(GDDz)

    all_param = dict(
        Min_Abs=Min_Abs,
        Max_Abs=Max_Abs,
        Max_Avg=Max_Avg,
        Min_Avg=Min_Avg,
        GDD=GDD,
        GDDz=GDDz,
        GDDlz=GDDlz,
        Ar=Ar,
        Evr=Evr,
        GrS=GrS,
        GAr=GAr,
        GInt=GInt
        )
    
    if opt('land_type') != 'Prentice':
        if opt('pas_ice_def') in ('ice', 'ice_noadj'):
            verb('    Finding ice cover')
            ice = data['ice']   #ice and snow thickness in m
            Min_Ice = np.amin(ice, 0)
            if opt('interp_scale') and opt('topo_map') and opt('pas_ice_def') == 'ice':     #hacks for interpolation with temperature adjustments
                tice = data['tice']
                Max_tice = np.amax(tice,0)
                Min_Ice = np.where(Max_tice > 0, 0, Min_Ice)
                Min_Ice = np.where(Min_Ice > 0.1, np.where(Max_tice > 0, 0, Min_Ice),
                                np.where(Max_tice < 0, np.where(np.amax(ice, 0) > 0.1, 0.2, Min_Ice), Min_Ice))
            all_param['Min_Ice_Land'] = Min_Ice
        elif opt('pas_ice_def') == 'maxt':
            all_param['Max_tice'] = np.amax(data['tice'], 0)
        if opt('land_type') == 'Pasta' and opt('pas_boil_pres'):
            verb('    Determining boiling point')
            ps = data['ps']  #surface pressure in hpa
            boilp = np.where(maxt > 108.3, 10**(10.26509-1810.94/(244.485+maxt)), 10**(10.19621-1730.63/(233.426+maxt)))     #antoine equation
            boil = np.where(maxt > 0, np.where(ps*100 < boilp, 1, 0), 0)
            boil = np.amax(boil,0)
            all_param['boil'] = boil

    
    return all_param

## Prentice et al. 1992 biomes

def Prentice_Alg(par):
    clim=0
    Min_Abs = par['Min_Abs']
    Max_Avg = par['Max_Avg']
    GDD = par['GDD']
    GDDz = par['GDDz']
    Ar = par['Ar']+0.05 #adjustment
    
    TropEver=False
    TropRain=False
    WarmTempEver=False
    TempSummer=False
    CoolTempCon=False
    BorEverCon=False
    BorSummer=False
    Sclero=False
    WarmGrass=False
    CoolGrass=False
    ColdGrass=False
    HotShrub=False
    ColdShrub=False
    NoneType=True

    th_trop = 10        #min for tropical, max for temperate summergreen
    th_warm = 0         #min for warm temperate and sclerophill, max for cool temperate and boreal summergreen
    th_temp = -25       #min for temperate summergreen
    th_cool = -32       #min for cool temperate
    th_bor_hi = -10     #max for boreal evergreen
    th_bor_lo = -45     #min for boreal evergreen

    if opt('temp_tunings') == 'tclim':
        th_trop = 9
        th_warm = 2
        th_temp = -20
        th_cool = -25
        th_bor_hi = -5
        th_bor_lo = -39

    elif opt('temp_tunings') == 'tavg':
        th_trop = 15.5
        th_warm = 5
        th_temp = -15
        th_cool = -19
        th_bor_hi = -2
        th_bor_lo = -35

    if Min_Abs > th_trop:
        if Ar > 0.8:
            TropEver = True
        if Ar > 0.45 and Ar < 0.95:
            TropRain = True
    if Min_Abs > th_warm:
        if Ar > 0.65:
            WarmTempEver = True
        if Ar > 0.28:
            Sclero = True
    if Min_Abs > th_temp:
        if Min_Abs < th_trop and Ar > 0.65 and GDD > 1200:
            TempSummer = True
    if Min_Abs > th_cool:
        if Min_Abs < th_warm and Ar > 0.65 and GDD > 900:
            CoolTempCon = True
    if Min_Abs > th_bor_lo:
        if Min_Abs < th_bor_hi and Ar > 0.75 and GDD > 350:
            BorEverCon = True
    if Min_Abs < th_warm and Ar > 0.65 and GDD > 350:
        BorSummer = True
    if Max_Avg > 22:
        HotShrub = True
        if Ar > 0.18:
            WarmGrass = True
    if Ar > 0.33:
        if GDD > 500:
            CoolGrass = True
        if GDDz > 100:
            ColdGrass = True
    if GDDz > 100:
        ColdShrub = True

    if TropEver or TropRain:
        WarmTempEver=False
        TempSummer=False
        CoolTempCon=False
        BorEverCon=False
        BorSummer=False
        Sclero=False
        WarmGrass=False
        CoolGrass=False
        ColdGrass=False
        HotShrub=False
        ColdShrub=False
        NoneType=False
    elif WarmTempEver:
        TempSummer=False
        CoolTempCon=False
        BorEverCon=False
        BorSummer=False
        Sclero=False
        WarmGrass=False
        CoolGrass=False
        ColdGrass=False
        HotShrub=False
        ColdShrub=False
        NoneType=False
    elif TempSummer or CoolTempCon or BorEverCon or BorSummer:
        Sclero=False
        WarmGrass=False
        CoolGrass=False
        ColdGrass=False
        HotShrub=False
        ColdShrub=False
        NoneType=False
    elif Sclero:
        WarmGrass=False
        CoolGrass=False
        ColdGrass=False
        HotShrub=False
        ColdShrub=False
        NoneType=False
    elif WarmGrass:
        CoolGrass=False
        ColdGrass=False
        HotShrub=False
        ColdShrub=False
        NoneType=False
    elif CoolGrass or ColdGrass:
        HotShrub=False
        ColdShrub=False
        NoneType=False
    elif HotShrub:
        ColdShrub=False
        NoneType=False
    elif ColdShrub:
        NoneType=False

    if TropEver:
        if TropRain:
            clim = PBTropSeason
        else:
            clim = PBTropRain
    elif TropRain:
        clim = PBTropDry
    elif WarmTempEver:
        clim = PBWarmMixed
    elif TempSummer:
        if CoolTempCon and BorSummer:
            if BorEverCon:
                clim = PBCoolMixed
            else:
                clim = PBTempDecid
    elif CoolTempCon:
        if BorSummer:
            if BorEverCon:
                clim = PBCoolCon
            else:
                clim = PBColdMixed
    elif BorEverCon:
        if BorSummer:
            clim = PBTaiga
    elif BorSummer:
        clim = PBColdDecid
    elif Sclero:
        clim = PBXero
    elif WarmGrass:
        clim = PBWarmGrass
    elif CoolGrass:
        if ColdGrass:
            clim = PBCoolGrass
    elif ColdGrass:
        clim = PBTundra
    elif HotShrub:
        clim = PBHotDesert
    elif ColdShrub:
        clim = PBSemiDesert
    elif NoneType:
        clim = PBIce

    return clim

Clim_func['Prentice'] = (Biome_Data, Biome_Param, Prentice_Alg)


## Pasta Bioclimate system

def Pasta_Alg(par):
    clim=0
    Min_Abs = par['Min_Abs']
    Max_Abs = par['Max_Abs']
    GDD = par['GDD']
    GDDz = par['GDDz']
    Ar = par['Ar']
    GAr = par['GAr']
    Evr = par['Evr']
    GrS = par['GrS']
    GInt = par['GInt']


    #Climate Thresholds:
    # GDD:
    th_XG = 50     #barren; measured in GDDz
    th_XF = 350     #marginal
    th_CE = 1300    #boreal

    # GInt
    th_XT = opt('pas_gint_thresh')    #peritropical; 1250 by default

    #winter minimum
    th_cool = 10    #cool winter
    th_cold = -10   #cold winter
    th_frigid = -40 #frigid winter

    #summer maximum
    th_hot = 50     #hot summer
    th_torrid = 70  #torrid summer
    th_boil = 100   #boiling summer (if not calculated from pressure)

    #Aridity
    th_XXr = 0.9    #rainforest
    th_XXf = 0.75   #forest
    th_XXs = 0.5    #moist; measured in GAr
    th_XA = 0.2     #semiarid
    th_Ad = 0.06    #semidesert

    #Growth supply
    th_XM = 0.8     #mediterranean
 
    #evaporation ratio
    th_XXp = 0.45   #pluvial
    th_TXrp = 0.4   #hyperpluvial rainforest

    #ice thickness
    th_CI = 0.1     #ice
    th_CI_t = 0     #temp threshold for alternative ice cover definitions

    if opt('temp_tunings') == 'tclim':    #alternate tuning for terraclim Earth data
        th_cool = 10
        th_cold = -4
        th_frigid = -35
        th_XM = 1.15
    
    elif opt('temp_tunings') == 'tavg': #alternate tuning for average monthly temperatures
        th_cool = 17
        th_cold = 0
        th_frigid = -30
        th_hot = 40
        th_torrid = 60
        th_boil = 90
    
    if opt('pas_med_thresh') > 0:   #override set med thresholds
        th_XM = opt('pas_med_thresh')

    if opt('land_subtype') in ('no_pluv', 'earthlike_no_pluv'):  #disable pluvial zones
        th_XXp = -1
        th_TXrp = -1
    
    if opt('land_subtype') in ('earthlike', 'earthlike_no_pluv'):    #disable non-earthlike zones
        boiling = False
        th_hot = 1000
        th_torrid = 1000
        if Min_Abs > th_cool:
            GInt = 0
    elif opt('pas_boil_pres'):
        boiling = par['boil']   #boiling summer test
    else:
        boiling = Max_Abs > th_boil
    
    if opt('pas_ice_def') in ('ice', 'ice_noadj'):
        ice = par['Min_Ice_Land'] > th_CI
    elif opt('pas_ice_def') == 'maxt':
        ice = par['Max_tice'] < th_CI_t
    else:
        ice = par['Max_Avg'] < th_CI_t



    XG = GDDz < th_XG     #barren test
    cold = Min_Abs < th_cold     #cold winter test
    torrid = boiling or Max_Abs > th_torrid      #torrid summer test

    if ice:     #ice test
        clim = CI
    elif Ar < th_XA and not XG:    #arid test
        Ahx = Ar < th_Ad     #hyperarid test
        if cold:
            if torrid:
                if Ahx:
                    clim = Ahe
                else:
                    clim = Ade
            else:
                if Ahx:
                    clim = Ahc
                else:
                    clim = Adc
        elif torrid:
            if Ahx:
                clim = Ahh
            else:
                clim = Adh
        else:
            if Ahx:
                clim = Aha
            else:
                clim = Ada
    else:
        XT = GInt < th_XT    #peritropical test
        XF = not XT and GDD < th_XF  #marginal test
        mild = Min_Abs > th_cool     #mild winter test
        warm = not boiling and Max_Abs < th_hot     #warm summer test
        XA = GAr < th_XXs  #semiarid test
        Xs = Ar < th_XXf  #savanna test
        Xxp = Evr < th_XXp  #pluvial test
        XM = GrS < th_XM    #med test
        if mild and warm:    #tropical test
            if XG:
                clim = TG
            elif XF:
                clim = TF
            elif XT:
                if XA:
                    if Xxp:
                        clim = TUAp
                    else:
                        clim = TUA
                elif Xs:
                    if Xxp:
                        clim = TUsp
                    else:
                        clim = TUs
                elif Ar < th_XXr:  #rainforest/seasonal test
                    if Xxp:
                        clim = TUfp
                    else:
                        clim = TUf
                else:
                    if Evr < th_TXrp: #hyperpluvial rainforest test
                        clim = TUrp
                    else:
                        clim = TUr
            else:
                if XA:
                    if Xxp:
                        clim = TQAp
                    else:
                        clim = TQA
                elif Xs:
                    if Xxp:
                        clim = TQsp
                    else:
                        clim = TQs
                else:
                    if Xxp:
                        clim = TQfp
                    else:
                        clim = TQf
        elif warm:  #cold test
            if XG:
                clim = CG
            elif XF:
                if cold:
                    clim = CFb
                else:
                    clim = CFa
            elif XA:
                if XM:
                    if cold:
                        clim = CAMb
                    else:
                        clim = CAMa
                elif cold:
                    if Xxp:
                        clim = CAbp
                    else:
                        clim = CAb
                else:
                    if Xxp:
                        clim = CAap
                    else:
                        clim = CAa
            elif XM:
                if cold:
                    clim = CMb
                else:
                    clim = CMa
            elif XT and not cold:
                if Xs:
                    if Xxp:
                        clim = CTsp
                    else:
                        clim = CTs
                else:
                    if Xxp:
                        clim = CTfp
                    else:
                        clim = CTf
            elif GDD > th_CE:    #temperate test
                if cold:
                    if Xxp:
                        clim = CDbp
                    else:
                        clim = CDb
                else:
                    if Xxp:
                        clim = CDap
                    else:
                        clim = CDa
            else:
                if Min_Abs < th_frigid:    #frigid winter test
                    if Xxp:
                        clim = CEcp
                    else:
                        clim = CEc
                elif cold:
                    if Xxp:
                        clim = CEbp
                    else:
                        clim = CEb
                else:
                    if Xxp:
                        clim = CEap
                    else:
                        clim = CEa
        elif mild:  #hot test
            if XG:
                clim = HG
            elif XF:
                if boiling:
                    clim = HFc
                elif torrid:
                    clim = HFb
                else:
                    clim = HFa
            elif XA:
                if XM:
                    if boiling:
                        clim = HAMc
                    elif torrid:
                        clim = HAMb
                    else:
                        clim = HAMa
                elif boiling:
                    if Xxp:
                        clim = HAcp
                    else:
                        clim = HAc
                elif torrid:
                    if Xxp:
                        clim = HAbp
                    else:
                        clim = HAb
                else:
                    if Xxp:
                        clim = HAap
                    else:
                        clim = HAa
            elif XM:
                if boiling:
                    clim = HMc
                elif torrid:
                    clim = HMb
                else:
                    clim = HMa
            elif XT and not torrid:
                if Xs:
                    if Xxp:
                        clim = HTsp
                    else:
                        clim = HTs
                else:
                    if Xxp:
                        clim = HTfp
                    else:
                        clim = HTf
            else:
                if boiling:
                    if Xxp:
                        clim = HDcp
                    else:
                        clim = HDc
                elif torrid:
                    if Xxp:
                        clim = HDbp
                    else:
                        clim = HDb
                else:
                    if Xxp:
                        clim = HDap
                    else:
                        clim = HDa
        else:  #extraseasonal test
            hyper = torrid or cold     #hyperseasonal test
            if XG:
                clim = EG
            elif XF:
                if hyper:
                    clim = EFb
                else:
                    clim = EFa
            elif XA:
                if XM:
                    if hyper:
                        clim = EAMb
                    else:
                        clim = EAMa
                elif hyper:
                    if Xxp:
                        clim = EAbp
                    else:
                        clim = EAb
                else:
                    if Xxp:
                        clim = EAap
                    else:
                        clim = EAa
            elif XM:
                if hyper:
                    clim = EMb
                else:
                    clim = EMa
            elif XT and not hyper:
                if Xs:
                    if Xxp:
                        clim = ETsp
                    else:
                        clim = ETs
                else:
                    if Xxp:
                        clim = ETfp
                    else:
                        clim = ETf
            else:
                if hyper:
                    if Xxp:
                        clim = EDbp
                    else:
                        clim = EDb
                else:
                    if Xxp:
                        clim = EDap
                    else:
                        clim = EDa
    return clim
                    

Clim_func['Pasta'] = (Biome_Data, Biome_Param, Pasta_Alg)

##Koppen-Geiger Unproxied

def Unproxied_Alg(par):
    clim=0
    Min_Abs = par['Min_Abs']
    GDD = par['GDD']
    Ar = par['Ar']
    GAr = par['GAr']
    Evr = par['Evr']
    GrS = par['GrS']

    #Climate Thresholds:
    # GDD:
    th_E = 350
    th_Xxa = 2400
    th_Xxb = 1300

    #winter minimum
    th_cool = 10    #cool winter
    th_cold = -10   #cold winter
    th_frigid = -40 #frigid winter

    #Aridity
    th_Tf = 0.92
    th_Tm_Ar = 0.85
    th_Tm_Evr = 0.45    #evaporation ratio
    th_B = 0.32
    th_BW = 0.14

    #Growth Supply
    th_Xs = 0.8

    #ice thickness
    th_EF = 0.1     #ice
    th_EF_t = 0     #temp threshold for alternative ice cover definitions

    if opt('temp_tunings') == 'tclim':    #alternate tuning for terraclim Earth data
        th_cool = 10
        th_cold = -4
        th_frigid = -35
        th_Xs = 1.15
    
    elif opt('temp_tunings') == 'tavg': #alternate tuning for average monthly temperatures
        th_cool = 17
        th_cold = 0
        th_frigid = -30
    
    if opt('pas_med_thresh') > 0:   #override set med thresholds
        th_Xs = opt('pas_med_thresh')
    
    if opt('pas_ice_def') in ('ice', 'ice_noadj'):
        ice = par['Min_Ice_Land'] > th_EF
    elif opt('pas_ice_def') == 'maxt':
        ice = par['Max_tice'] < th_EF_t
    else:
        ice = par['Max_Avg'] < th_EF_t

    #Groups
        
    if ice:
        clim = E
    elif Ar < th_B:
        clim = B
    elif GDD < th_E:
        clim = E
    elif Min_Abs > th_cool:
        clim = A
    elif Min_Abs > th_cold:
        clim = C
    else:
        clim = D

    if opt('land_subtype') == 'groups': #finish here for groups only
        return clim
                
    #Full Koppen set
        
    else:
        
        if clim == B:
            if Ar < th_BW:    #desert/steppe test
                clim = BW
            else:
                clim = BS
            if opt('land_subtype') != 'two_letter':
                Bxh = Min_Abs > th_cold   
                if clim == BW:
                    if Bxh:
                        clim = BWh
                    else:
                        clim = BWk
                else:
                    if Bxh:
                        clim = BSh
                    else:
                        clim = BSk
        elif clim == E:
            if ice:
                clim = EF
            else:
                clim = ET
        elif clim == A:
            if Ar > th_Tf:
                clim = Af
            elif Ar > th_Tm_Ar or Evr < th_Tm_Evr:
                clim = Am
            else:
                if GAr < Ar:
                    clim = As
                else:
                    clim = Aw
        else:
            Xs = GrS < th_Xs
            Xw = GAr > Ar * 1.02
            if Xs and Xw:
                if opt('kg_wet_summer_priority'):   #determine med or wet-summer where they overlap
                    Xs = False
                else:
                    Xw = False
            if GDD < th_Xxb:
                if Min_Abs > th_frigid:
                    tlet = 3
                else:
                    tlet = 4
            else:
                if GDD > th_Xxa:
                    tlet = 1
                else:
                    tlet = 2
            if opt('land_subtype') != 'reduced':    #only bother with this if not doing reduced sets
                if opt('land_subtype') == 'two_letter':
                    if clim == C:
                        if Xs:
                            clim = Cs
                        elif Xw:
                            clim = Cw
                        else:
                            clim = Cf
                    else:
                        if Xs:
                            clim = Ds
                        elif Xw:
                            clim = Dw
                        else:
                            clim = Df
                elif clim == C:       #choose within C and D by overlap of subtypes
                    if Xs:
                        if tlet == 1:
                            clim = Csa
                        elif tlet == 2:
                            clim = Csb
                        else:
                            clim = Csc
                    elif Xw:
                        if tlet == 1:
                            clim = Cwa
                        elif tlet == 2:
                            clim = Cwb
                        else:
                            clim = Cwc
                    else:
                        if tlet == 1:
                            clim = Cfa
                        elif tlet == 2:
                            clim = Cfb
                        else:
                            clim = Cfc
                else:
                    if Xs:
                        if tlet == 1:
                            clim = Dsa
                        elif tlet == 2:
                            clim = Dsb
                        elif tlet == 3:
                            clim = Dsc
                        else:
                            clim = Dsd
                    elif Xw:
                        if tlet == 1:
                            clim = Dwa
                        elif tlet == 2:
                            clim = Dwb
                        elif tlet == 3:
                            clim = Dwc
                        else:
                            clim = Dwd
                    else:
                        if tlet == 1:
                            clim = Dfa
                        elif tlet == 2:
                            clim = Dfb
                        elif tlet == 3:
                            clim = Dfc
                        else:
                            clim = Dfd
                            
        if opt('land_subtype') == 'reduced':
            if clim == Af:      #replace climates with reduced equivalents
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
            elif clim == C:      #C and D use more simple reduced breakdown
                if Xs:
                    clim = Med
                elif tlet == 1:
                    clim = Subtropical
                else:
                    clim = Oceanic
            elif clim == D:
                if tlet < 3:
                    clim = Continental
                else:
                    clim = Subarctic
            elif clim == ET:
                clim = Tundra
            elif clim == EF:
                clim = IceCap
                
    return clim                    

Clim_func['KG_unproxied'] = (Biome_Data, Biome_Param, Unproxied_Alg)


## Sea zones

#Functions for sea zones work slightly differently, taking outputs from land functions as input;
# new data is only retrieved if not already present in the land data, and any new output is added to the same dictionary as the land outputs

def Sea_Data(dat, data):
    if opt('sea_type') == 'sea_none' or opt('sea_subtype') == 'flat':    
        return data
    
    elif opt('sea_subtype') == 'full' or opt('sea_ice_use_temp'):
        if opt('sea_use_ts'):
            tkey = 'ts'
        else:
            tkey = 'tas'
        if tkey not in data:
            if opt('interp_scale') and opt('topo_map'):
                if 'adjust' not in data:
                    tkey_dat, adjust = Get_nc_adjust(dat, tkey, 'grnz')
                    data[tkey] = tkey_dat = adjust + opt('temp_adjust')    #convert to C and add ajustment
                else:
                    tkey_dat, was_in = Get_nc_if(dat, tkey, data, tkey)
                    if not was_in:
                        tkey_dat += data['adjust'] + opt('temp_adjust')    #convert to C and add ajustment
                        data[tkey] = tkey_dat
            else:
                tkey_dat, was_in = Get_nc_if(dat, tkey, data, tkey)
                if not was_in:
                    tkey_dat += opt('temp_adjust')    #convert to C and add ajustment
                    data[tkey] = tkey_dat
            
    if not opt('sea_ice_use_temp'):
        get_mask(dat, 'lsm', convert=True)
        data['sic'], was_in = Get_nc_if(dat, 'sic', data, 'sic', dummy_ice = opt('dummy_ice'))
        
    return data

def Sea_Param(data, par):
    if opt('sea_type') == 'sea_none' or opt('sea_subtype') == 'flat':
        return par
    
    elif opt('sea_subtype') == 'full' or opt('sea_ice_use_temp'):
        if opt('sea_use_ts'):
            ts = data['ts']     #surface temperature in C
        else:
            ts = data['tas']    #2-meter air temperature in C
        par['Min_Seatemp'] = np.amin(ts, axis=0)
    if opt('sea_ice_use_temp'):
        par['Max_Seatemp'] = np.amax(ts, axis=0)
        
    else:
        sic = data['sic']   #proportional sea ice cover
        if opt('seasonless'):
            par['Avg_Ice'] = np.mean(sic, axis=0)
        else:
            par['Max_Ice'] = np.amax(sic, axis=0)
            par['Min_Ice'] = np.amin(sic, axis=0)

    return par

def Sea_Alg(par):
    clim=0
    if opt('sea_type') == 'sea_none' or opt('sea_subtype') == 'flat':
        return SeaFlat

    else:
        clim = SeaTemp  #SeaTemp is default
        
        if opt('sea_subtype') == 'full' or opt('sea_ice_use_temp'):
            Min_Seatemp = par['Min_Seatemp']
        
            if opt('sea_subtype') == 'full' and Min_Seatemp > 18:
                clim = SeaTrop
                
            elif opt('sea_ice_use_temp'):
                Max_Seatemp = par['Max_Seatemp']
                if Min_Seatemp < -2:
                    if Max_Seatemp < -2:
                        clim = SeaPermIce
                    else:
                        clim = SeaSeasonalIce
                    
        if not opt('sea_ice_use_temp'):
            if opt('seasonless'):
                Avg_Ice = par['Avg_Ice']
                if Avg_Ice > 0.5:
                    clim = SeaPermIce
            else:
                Max_Ice = par['Max_Ice']
                Min_Ice = par['Min_Ice']
                if Max_Ice > 0.2:
                    if Min_Ice > 0.8:
                        clim = SeaPermIce
                    else:
                        clim = SeaSeasonalIce

    return clim  

Clim_func['sea_standard'] = (Sea_Data, Sea_Param, Sea_Alg)
Clim_func['sea_none'] = Clim_func['sea_standard']

#Alternative efficient algorithm that works on whole array at once:

def Sea_Alg_efficient(par):
    for k, v in par.items():
        clim = np.full_like(v, SeaFlat) #copy first available parameter for shape
        break
    if opt('sea_type') == 'sea_none' or opt('sea_subtype') == 'flat':
        return clim

    else:
        clim = np.full_like(clim, SeaTemp)  #SeaTemp is default
        
        if opt('sea_subtype') == 'full' or opt('sea_ice_use_temp'):
            Min_Seatemp = par['Min_Seatemp']
            if opt('sea_subtype') == 'full':
                clim = np.where(Min_Seatemp > 18, SeaTrop, clim)
                
            elif opt('sea_ice_use_temp'):
                Max_Seatemp = par['Max_Seatemp']
                clim = np.where(Min_Seatemp < -2,
                                np.where(Max_Seatemp < -2, SeaPermIce,
                                         SeaSeasonalIce),
                                clim)
                    
        if not opt('sea_ice_use_temp'):
            Max_Ice = par['Max_Ice']
            Min_Ice = par['Min_Ice']
            clim = np.where(Max_Ice > 0.2,
                            np.where(Min_Ice > 0.8, SeaPermIce,
                                     SeaSeasonalIce),
                            clim)

    return clim

Clim_func['sea_standard_efficient'] = (Sea_Data, Sea_Param, Sea_Alg_efficient)


## Pasta Ocean Bioclimate system

def Pasta_Sea_Data(dat, data):
    if opt('sea_subtype') != 'no_trop' or opt('sea_ice_use_temp'):
        if opt('sea_use_ts'):
            tkey = 'ts'
        else:
            tkey = 'tas'
        if opt('interp_scale') and opt('topo_map'):
            if 'adjust' not in data:
                tkey_dat, adjust = Get_nc_adjust(dat, tkey, 'grnz')
                data[tkey] = tkey_dat = adjust + opt('temp_adjust')    #convert to C and add ajustment
            else:
                tkey_dat, was_in = Get_nc_if(dat, tkey, data, tkey)
                if not was_in:
                    tkey_dat += data['adjust'] + opt('temp_adjust')    #convert to C and add ajustment
                data[tkey] = tkey_dat
        else:
            tkey_dat, was_in = Get_nc_if(dat, tkey, data, tkey)
            if not was_in:
                tkey_dat += opt('temp_adjust')    #convert to C and add ajustment
            data[tkey] = tkey_dat

        if opt('sea_subtype') == 'full' and opt('gdd_limit_light'):
            rss, was_in = Get_nc_if(dat, 'rss', data, 'rin')
            if not was_in:
                ssru = Get_nc(dat, 'ssru')
                data['rin'] = rss - ssru

    if not opt('sea_ice_use_temp'):
        if opt('pas_ice_def') in ('ice', 'ice_noadj'):
            sic = Get_nc(dat, 'sic', dummy_ice = opt('dummy_ice'))
            if opt('interp_scale'):
                if 'ice' in data:
                    ice = data['ice']
                else:
                    snd = Get_nc(dat, 'snd', no_interp = True)
                    ice = np.minimum(snd, 0.4)
                    ice = Interp(ice)
                get_mask(dat, 'lsm', True)
                mask = common('mask_big')
                sic = np.where(mask, np.maximum(sic, ice*10), sic)  #on land areas treat 1 cm of ice as equivalent to 10% sea ice cover
            data['sic'] = sic
        else:
            data['sic'], was_in = Get_nc_if(dat, 'sic', data, 'sic', dummy_ice = opt('dummy_ice'))
            
    return data

def Pasta_Sea_Param(data, par):

    ##GDD config

    GDDlzbase = 0   #values for light-limited GDDz maximum; used only if GDDLz not already calculated
    GDDlzplats = 100 / opt('gdd_par_ratio')

    if opt('sea_subtype') != 'no_trop' or opt('sea_ice_use_temp'):
        if opt('sea_use_ts'):
            ts = data['ts']     #surface temperature in C
        else:
            ts = data['tas']    #2-meter air temperature in C
        par['Min_Seatemp'] = np.amin(ts, axis=0)
        par['Max_Seatemp'] = np.amax(ts, axis=0)

        if opt('sea_subtype') == 'full' and opt('gdd_limit_light'):
            if 'GDDlz' not in data:
                rin = data['rin']
                GDDlz = Calc_GDD(rin, base=GDDlzbase, plat_start=GDDlzplats) / 10
                par['GDDlz'] = Calc_GDD_total(GDDlz)
        
    if not opt('sea_ice_use_temp'):
        sic = data['sic']   #ice thickness (or proportional cover)
        par['Max_Ice'] = np.amax(sic, axis=0)
        par['Min_Ice'] = np.amin(sic, axis=0)

    return par

def Pasta_Sea_Alg(par):
    clim=0
    if opt('sea_subtype') == 'no_trop' and not opt('sea_ice_use_temp'):
        Min_Seatemp = 10
        Max_Seatemp = 10
    else:
        Min_Seatemp = par['Min_Seatemp']
        Max_Seatemp = par['Max_Seatemp']
    if not opt('sea_ice_use_temp'):
        Max_Ice = par['Max_Ice']
        Min_Ice = par['Min_Ice']
    else:
        Max_Ice = np.where(Min_Seatemp < -2, 1, 0)
        Min_Ice = np.where(Max_Seatemp < -2, 1, 0)
    if opt('sea_subtype') == 'full' and opt('gdd_limit_light'):
        GDDlz = par['GDDlz']
    else:
        GDDlz = 10000

    #Thresholds

    th_f = 0.2  #threshold for seasonal ice cover
    th_fi = 0.8 #threshold for permanent ice cover
    th_g = 50   #threshold for dark oceans (GDDlz)
    th_t = 18   #threshold for tropical ocean (min C)
    th_h = 40   #threshold for hot ocean (max C)
    th_r = 60   #threshold for torrid ocean (max C)

    if opt('sea_subtype') != 'full':
        th_h = 1000
        th_r = 1000
        if opt('sea_subtype') == 'no_trop':
            th_t = 1000
    
    if opt('seasonless'):
        th_f = 0.5
        th_fi = 0.5
    
    if Max_Ice > th_f:
        if Min_Ice > th_fi:
            clim = Ofi
        elif GDDlz < th_g:
            clim = Ofg
        else:
            clim = Ofd
    elif GDDlz < th_g:
        clim = Og
    elif Max_Seatemp > th_r:
        clim = Or
    elif Max_Seatemp > th_h:
        if Min_Seatemp > th_t:
            clim = Oh
        else:
            clim = Oe
    elif Min_Seatemp > th_t:
        clim = Ot
    else:
        clim = Oc

    return clim  

Clim_func['sea_Pasta'] = (Pasta_Sea_Data, Pasta_Sea_Param, Pasta_Sea_Alg)

## Special functions to attempt to ensure that data for uses outside algorithms is always available

def Extra_Data(dat, data):
    if opt('sea_type') != 'none':
        try:
            mask = get_mask(dat, 'lsm', convert=True)
            if 'mask' not in data:
                data['mask'] = mask
        except:
            pass
    if opt('make_chart'):
        try:
            if 'tas' not in data:
                if opt('interp_scale') and opt('topo_map'):
                    if 'adjust' not in data:
                        lat, lon = coords_from_file(dat,'lat','lon')
                        coords = (lat,lon)
                        tkey_dat, adjust = Get_nc_adjust(dat, 'tas', 'grnz', coords)
                        data['tas'] = tkey_dat = adjust + opt('temp_adjust')    #convert to C and add ajustment
                    else:
                        tkey_dat, was_in = Get_nc_if(dat, 'tas', data, 'tas')
                        if not was_in:
                            tkey_dat += data['adjust'] + opt('temp_adjust')    #convert to C and add ajustment
                            data['tas'] = tkey_dat
                else:
                    tkey_dat, was_in = Get_nc_if(dat, 'tas', data, 'tas')
                    if not was_in:
                        tkey_dat += opt('temp_adjust')    #convert to C
                        data['tas'] = tkey_dat
            pr, was_in = Get_nc_if(dat, 'pr', data, 'pr')
            if not was_in:
                if np.amax(pr) < 1: #try to avoid converting if not m/s as presumed
                    pr*=opt('precip_adjust')  #convert from m/s to mm/month
                data['pr'] = pr
        except:
            pass
    return data

def Extra_Param(data, par):
    if opt('make_chart'):
        try:
            if 'Avg_Temp' not in par:
                par['Avg_Temp'] = np.mean(data['tas'], 0)
            if 'Total_Precip' not in par:
                par['Total_Precip'] = np.mean(data['pr'],0)*12  #convert to mm/year
        except:
            pass
    return par


## Runs main routines if run directly

if __name__ == "__main__":

    files, in_opts = Get_input()

    Make_map(files, in_opts)

    while True:
        print('')
        if Prompt_bool('Run again with same files? (y/n): '):
            fi, in_opts = Get_input(False)
            Make_map(files, in_opts)
        else:
            break
    

