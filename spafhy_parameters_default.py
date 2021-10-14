# -*- coding: utf-8 -*-
"""
PARAMETERS OF SPAFHY

Created on Mon Jun 25 18:34:12 2018

@author: slauniai
"""
#pathlib.Path converts '/' to '\' on windows,
#i.e. one can always use '/' in file path
# import pathlib
#SpaFHy  v. 250618 parameter files for different catchments

def parameters(catchment, scenario):
    pgen = {'catchment_id': None,
            'prefix': 'pur',
            'gis_folder': r'/content/PURUVESI_input_data/'+catchment+r'/',
            'gis_scenarios': r'/content/PURUVESI_input_data/'+catchment+r'/' +scenario + r'/',
            'forcing_file': r'/content/PURUVESI_input_data/puruvesi_saa_19810101-2020062.csv',
            'output_folder': r'/content/output/'+catchment+r'/'+ scenario +r'/',
            'runoff_file': '', #
            'ncf_file': r'/content/output/'+catchment+r'/'+ scenario +r'/'+'ch.nc',
            'start_date': '2018-01-01',
            'end_date': '2028-12-31',
            'spinup_end': '2013-12-31',
            'dt': 86400.0,
            'spatial_cpy': True,
            'spatial_soil': True     
           }


    # canopygrid
    pcpy = {'loc': {'lat': 61.4, 'lon': 23.7},
            'flow' : { # flow field
                     'zmeas': 2.0,
                     'zground': 0.5,
                     'zo_ground': 0.01
                     },
            'interc': { # interception
                        'wmax': 1.5, # storage capacity for rain (mm/LAI)
                        'wmaxsnow': 4.5, # storage capacity for snow (mm/LAI),
                        },
            'snow': {
                    # degree-day snow model
                    'kmelt': 2.8934e-05, # melt coefficient in open (mm/s)
                    'kfreeze': 5.79e-6, # freezing coefficient (mm/s)
                    'r': 0.05 # maximum fraction of liquid in snow (-)
                    },
                    
            'physpara': {
                        # canopy conductance
                        'amax': 10.0, # maximum photosynthetic rate (umolm-2(leaf)s-1)
                        'g1_conif': 2.1, # stomatal parameter, conifers
                        'g1_decid': 3.5, # stomatal parameter, deciduous
                        'q50': 50.0, # light response parameter (Wm-2)
                        'kp': 0.6, # light attenuation parameter (-)
                        'rw': 0.20, # critical value for REW (-),
                        'rwmin': 0.02, # minimum relative conductance (-)
                        # soil evaporation
                        'gsoil': 1e-2 # soil surface conductance if soil is fully wet (m/s)
                        },
            'phenopara': {
                        #seasonal cycle of physiology: smax [degC], tau[d], xo[degC],fmin[-](residual photocapasity)
                        'smax': 18.5, # degC
                        'tau': 13.0, # days
                        'xo': -4.0, # degC
                        'fmin': 0.05, # minimum photosynthetic capacity in winter (-)
                        # deciduos phenology
                        'lai_decid_min': 0.1, # minimum relative LAI (-)
                        'ddo': 45.0, # degree-days for bud-burst (5degC threshold)
                        'ddur': 23.0, # duration of leaf development (days)
                        'sdl': 9.0, # daylength for senescence start (h)
                        'sdur': 30.0, # duration of leaf senescence (days),
                         },
            'state': {
                       'lai_conif': 3.5, # conifer 1-sided LAI (m2 m-2)
                       'lai_decid_max': 0.5, # maximum annual deciduous 1-sided LAI (m2 m-2): 
                       'hc': 16.0, # canopy height (m)
                       'cf': 0.6, # canopy closure fraction (-)
                       #initial state of canopy storage [mm] and snow water equivalent [mm]
                       'w': 0.0, # canopy storage mm
                       'swe': 0.0, # snow water equivalent mm
                       }
            }

        
    # BUCKET
    pbu = {'depth': 0.4,  # root zone depth (m)
           # following soil properties are used if spatial_soil = False
           'poros': 0.43, # porosity (-)
           'fc': 0.33, # field capacity (-)
           'wp': 0.13,	 # wilting point (-)
           'ksat': 2.0e-6, 
           'beta': 4.7,
           #organic (moss) layer
           'org_depth': 0.04, # depth of organic top layer (m)
           'org_poros': 0.9, # porosity (-)
           'org_fc': 0.3, # field capacity (-)
           'org_rw': 0.24, # critical vol. moisture content (-) for decreasing phase in Ef
           'maxpond': 0.0, # max ponding allowed (m)
           #initial states: rootzone and toplayer soil saturation ratio [-] and pond storage [m]
           'rootzone_sat': 0.6, # root zone saturation ratio (-)
           'org_sat': 1.0, # organic top layer saturation ratio (-)
           'pond_sto': 0.0 # pond storage
           }
    
    # TOPMODEL
    ptop = {'dt': 86400.0, # timestep (s)
            'm': 0.01, # scaling depth (m)
            'ko': 0.001, # transmissivity parameter (ms-1)
            'twi_cutoff': 99.5,  # cutoff of cumulative twi distribution (%)
            'so': 0.05 # initial saturation deficit (m)
           }
    
    return pgen, pcpy, pbu, ptop

def soil_properties():
    """
    Defines 4 soil types: Fine, Medium and Coarse textured + organic Peat
    and Humus.
    Currently SpaFHy needs following parameters: soil_id, poros, dc, wp, wr,
    n, alpha, Ksat, beta
    """
    psoil = {
             'FineTextured': 
                 {'airentry': 34.2,
                  'alpha': 0.018,
                  'beta': 7.9,
                  'fc': 0.34,
                  'ksat': 1e-06,
                  'n': 1.16,
                  'poros': 0.5,
                  'soil_id': 3.0,
                  'wp': 0.25,
                  'wr': 0.07,
                 },

             'MediumTextured': 
                 {'airentry': 20.8,
                  'alpha': 0.024,
                  'beta': 4.7,
                  'fc': 0.33,
                  'ksat': 1e-05,
                  'n': 1.2,
                  'poros': 0.43,
                  'soil_id': 2.0,
                  'wp': 0.13,
                  'wr': 0.05,
                 },

            'CoarseTextured':
                 {'airentry': 14.7,
                  'alpha': 0.039,
                  'beta': 3.1,
                  'fc': 0.21,
                  'ksat': 0.0001,
                  'n': 1.4,
                  'poros': 0.41,
                  'soil_id': 1.0,
                  'wp': 0.1,
                  'wr': 0.05,
                 },

             'Peat':
                 {'airentry': 29.2,
                  'alpha': 0.123,
                  'beta': 6.0,
                  'fc': 0.414,
                  'ksat': 5e-05,
                  'n': 1.28,
                  'poros': 0.9,
                  'soil_id': 4.0,
                  'wp': 0.11,
                  'wr': 0.0,
                 },
              'Humus':
                 {'airentry': 29.2,
                  'alpha': 0.123,
                  'beta': 6.0,
                  'fc': 0.35,
                  'ksat': 8e-06,
                  'n': 1.28,
                  'poros': 0.85,
                  'soil_id': 5.0,
                  'wp': 0.15,
                  'wr': 0.01,
                 },
            }
    return psoil

def nutpara(pgen, catchment, scenario):
    catchn = pgen['catchment_id']
    #print (catchn, type(catchn))
    npara={
        'spafhyfile': pgen['output_folder'] +r'/ch.nc',                           #Spathy output in netCDF-format
        'mottisims' : r'/content/nutspafhy_operandum/stands.p',                                                  #Motti simulation parameters for growth and yield and biomass
        'nutspafhyfile': r'/nutspafhy.nc',   #NutSpahy nutrient balance data in netCDF-file
        'ofolder' : pgen['output_folder'] ,                           #outputfolder
        'npr_out' : pgen['output_folder'] + 'nutrient_export_load_'  +  catchment + '_'+scenario+'.csv',
    }    
    return npara

def immpara(catchment):
    #These are a priori calculated for SVE catchments 
    #input catchment number as string
    #ID	(iN_peat, iN_miner), (iP_peat,iP_miner)

    calibrated_imm ={
            '2':{'N':(0.868,0.903),'P':(0.909,0.931)},
            '10':{'N':(0.878,0.913),'P':(0.932,0.954)},
            '13':{'N':(0.834,0.878),'P':(0.848,0.882)},
            '14':{'N':(0.906,0.924),'P':(0.930,0.943)},
            '21':{'N':(0.834,0.877),'P':(0.872,0.899)},
            '22':{'N':(0.861,0.903),'P':(0.926,0.892)},
            '24':{'N':(0.807,0.855),'P':(0.757,0.807)},
            '25':{'N':(0.867,0.908),'P':(0.927,0.952)},
            '27':{'N':(0.908,0.937),'P':(0.963,0.949)},
            '31':{'N':(0.871,0.906),'P':(0.954,0.934)},
            '32':{'N':(0.889,0.919),'P':(0.750,0.800)},
            '33':{'N':(0.842,0.881),'P':(0.887,0.915)}
            }

    return calibrated_imm[catchment]['N'], calibrated_imm[catchment]['P']

def immpara_rmse():
    rmse={
            'N':(0.019, 0.019), 'P':(0.067, 0.051)}


    return rmse['N'], rmse['P']
