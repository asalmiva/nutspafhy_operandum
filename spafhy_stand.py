# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:57:55 2016
@author: slauniai


soil.biology -subpackage function

Contains functions and parameters related to biological activity of soil

"""

import numpy as np
from scipy.interpolate import interp1d   
from itertools import product
from pyproj import Proj, transform

def soilRespiration(GisData, soildata, Ts,  Wliq, imm_n, imm_p, decopara=None, soiltype='Yolo', limitpara=None, dt=1):
    """
    computes soil respiration rate (CO2-flux) based on Pumpanen et al. (2003) Soil.Sci.Soc.Am
    Restricts respiration as in Skopp et al. (1990),Soil.Sci.Soc.Am
    IN:
        minmask - minerla soil mask (array where mineral soil is 1)        
        Ts - soil temperature (degC)
        Wliq - soil vol. moisture content (m3 m-3)
        poros - soil porosity  (m3 m-3)
        decopara - [R10, Q10]
        limitpara - [a,b,d,g] of Skopp -model
    OUT:
        rsoil - soil respiration rate
        fm - relative modifier (Skopp et al.)
        co2mi - array of co2 efflux from mineral soil, kg/ha/day
    """
    
    nrows = int(((GisData['info'][0]).split()[1]))   #shape c,r
    ncols = int(((GisData['info'][1]).split()[1]))
    #ixm = np.where(np.ravel(GisData['smc'])==1)
    ixm = np.equal(GisData['smc'],1)

    ix_1 = np.where(np.logical_and(np.equal(GisData['sfc'],1),np.equal(GisData['smc'],1)))   # Lehto fertility
    ix_2 = np.where(np.logical_and(np.equal(GisData['sfc'],2),np.equal(GisData['smc'],1)))   # Lehtomainen fertility
    ix_3 = np.where(np.logical_and(np.equal(GisData['sfc'],3),np.equal(GisData['smc'],1)))   # Tuore fertility
    ix_4 = np.where(np.logical_and(np.equal(GisData['sfc'],4),np.equal(GisData['smc'],1)))   # Kuivahko fertility
    ix_5 = np.where(np.logical_and(np.greater_equal(GisData['sfc'],5),np.equal(GisData['smc'],1)))   # Kuiva fertility

    ix={'sfc_1': ix_1,'sfc_2': ix_2, 'sfc_3': ix_3, 'sfc_4': ix_4, 'sfc_5': ix_5}    
    N =  {'sfc_1': 2.4,'sfc_2': 2.2, 'sfc_3': 1.8, 'sfc_4': 1.6, 'sfc_5': 1.4}       # Tamminen 1991 folia forestalia 777, page 18 Table 16: N cont in OM % dm
    P =  {'sfc_1': 0.17,'sfc_2': 0.15, 'sfc_3': 0.13, 'sfc_4': 0.11, 'sfc_5': 0.1}            

    Ts = 16. if Ts>16. else Ts
    Ts=-45. if  Ts<-5. else Ts  
    co2mi = np.empty((ncols,nrows)); co2mi[:]=np.nan
    fm = np.empty((ncols,nrows)); fm[:]=np.nan
    Nrel = np.empty((ncols,nrows)); Nrel[:]=np.nan
    Prel = Nrel.copy()
    poros= soildata['poros'] #GisData['params']['poros']

    sp={'Yolo':[3.83, 4.43, 1.25, 0.854], 
        'Valentine': [1.65,6.15,0.385,1.03]} #Skopp param [a,b,d,g], Table 1 
    if decopara is None:
        Q10=2.3; R10=4.0; # (-), umol m-2s-1
    else:
        decopara=np.array(decopara, ndmin=1)
        Q10=decopara[:,1]; R10=decopara[:,0]
    
    if limitpara is None:
        p=sp[soiltype]
    else:
        p=limitpara
        
    #unrestricted rate    
    rs0=R10*np.power(Q10, (Ts-10.0)/10.0) #umol m-2 s-1
    rs0 = rs0 *1e-6 * 44.0 * 1e-3 *  1e4  * 86400  * 0.4   #, -> mol -> g (44 g/mol) -> kg -> ha ->day, to heterotrophic respiration 
    if Ts<-10. : rs0 = 0.0    
    #fm=np.minimum(p[0]*Wliq**p[2], p[1]*afp**p[3]) #]0...1]
    #fm=np.minimum(fm, np.ones(np.shape(fm)))
    fm=np.minimum(p[0]*Wliq[ixm]**p[2], p[1]*(poros[ixm]-Wliq[ixm])**p[3]) #]0...1]
    fm=np.minimum(fm, np.ones(np.shape(fm)))
    co2mi[ixm]=rs0 * fm * dt
    
    #N release
    C_in_OM = 0.55                                                              # C content in OM kg kg-1
    CO2_to_C = 12./44.
    Nmicrob = imm_n                                                               # microbial immobilisation    
    Pmicrob = imm_p                                                             # microbial immobilisation
    for k in ix.keys():
        Nrel[ix[k]] = co2mi[ix[k]] * CO2_to_C / C_in_OM * N[k] / 100. * (1.-Nmicrob)    #Nrelease kg ha-1 day-1 
        Prel[ix[k]] = co2mi[ix[k]] * CO2_to_C / C_in_OM * P[k] / 100. * (1.-Pmicrob)    #Nrelease kg ha-1 day-1 
    
    return Nrel, Prel, ixm


#def understory_biomass(site, age, x=[]):
def understory_uptake(lat, lon, ts, gisdata, expected_yield, simtime):
    """
    Created on Wed Jun 18 12:07:47 2014

    @author: slauniai

    understory_biomass(site, age, x=[]):    
    Computes understory biomasses using models of Muukkonen & Makipaa, 2006 Bor. Env. Res.\n
    INPUT:
        lat - latitude in YKJ or EUREF equivalent 
        lon - longitude 
        ts - annual temperature sum in degree days 
        gisdata - includes np arrays of catchment stand and site properties
        expected_yield of stand during the simulation period m3 ha-1
        simtime - simulation time in years
        x - array of independent variables (optional, if not provided age-based model is used):
            x[0]=lat (degN, in decimal degrees)
            x[1]=lon (degE in decimal degrees) 
            x[2]=elev (m)
            x[3]=temperature sum (degC)
            x[4]=site nutrient level (-) 
            x[5]=stem vol. (m3 ha-1)
            x[6]=stem nr (ha-1)
            x[7]=basal area (m2 ha-1)
            x[8]=site drainage status,integer
    OUTPUT:
        y - dry biomasses (kg ha-1) of different groups\n
    SOURCE:
        Muukkonen & Makipaa, 2006. Bor.Env.Res. 11, 355-369.\n
    AUTHOR:
        Samuli Launiainen 18.06.2014, Modified for array operations by Ari Laurén 13.4.2020 \n
    NOTE:
         Multi-regression models not yet tested!
         In model equations independent variables named differently to M&M (2006): here x[0] = z1, x[1]=z2, ... x[7]=z8 and x[8]=z10\n
         \n
         Site nutrient level x[4] at upland sites:
             1: herb-rich forest 
             2: herb-rich heat f. 
             3: mesic heath f. 
             4: sub-xeric heath f.
             5: xeric heath f. 
             6: barren heath f.
             7: rock,cliff or sand f. 
         Site nutrient level x[4] at mires:\n
             1: herb-rich hw-spruce swamps, pine mires, fens, 
             2: V.myrtillus / tall sedge spruce swamps, tall sedge pine fens, tall sedge fens,
             3: Carex clobularis / V.vitis-idaea swamps, Carex globularis pine swamps, low sedge (oligotrophic) fens,
             4: Low sedge, dwarf-shrub & cottongrass pine bogs, ombo-oligotrophic bogs,
             5: S.fuscum pine bogs, ombotrophic and S.fuscum low sedge bogs.
         Drainage status x[8] at mires (Paavilainen & Paivanen, 1995):
             1: undrained
             2: Recently draines, slight effect on understory veg., no effect on stand
             3: Transforming drained mires, clear effect on understory veg and stand
             4: Transforming drained mires, veget. resembles upland forest site type, tree-stand forest-like.
  
    """
        
    
    #------------- classify and map pixels-------------------------------------------------------- 
    ix_pine_upland = np.where( np.logical_and( np.logical_and(np.greater_equal( gisdata['p_vol'],gisdata['s_vol']), 
                                                              np.greater_equal( gisdata['p_vol'],gisdata['b_vol']))   ,
                                                              np.equal(gisdata['smc'], 1)))
    ix_spruce_upland = np.where( np.logical_and( np.logical_and(np.greater_equal( gisdata['s_vol'],gisdata['p_vol']), 
                                                              np.greater_equal( gisdata['s_vol'],gisdata['b_vol']))   ,
                                                              np.equal(gisdata['smc'], 1)))
    ix_broadleaved_upland = np.where( np.logical_and( np.logical_and(np.greater_equal( gisdata['b_vol'],gisdata['p_vol']), 
                                                              np.greater_equal( gisdata['b_vol'],gisdata['s_vol']))   ,
                                                              np.equal(gisdata['smc'], 1)))
    ix_spruce_mire = np.where(np.equal(gisdata['smc'], 2))
    ix_pine_bog = np.where(np.equal(gisdata['smc'], 3))
    ix_open_peat = np.where(np.equal(gisdata['smc'], 4))
    
    #---------------------------------------
    #latitude = 0.0897*lat/10000. + 0.3462                                       #approximate conversion to decimal degrees within Finland,  N
    #longitude = 0.1986*lon/10000. + 17.117                                      #approximate conversion to decimal degrees within Finland in degrees E
    inProj = Proj(init='epsg:3067')
    outProj = Proj(init='epsg:4326')
    longitude,latitude = transform(inProj,outProj,lon,lat)
    Nstems = 900.   # x6 number of stems -ha, default 900
    drain_s =4      # x8 drainage status, default value 4

    #---------------------------------------

    def gv_biomass_and_nutrients(gisdata, ix_pine_upland, ix_spruce_upland, ix_broadleaved_upland, ix_spruce_mire, ix_pine_bog,
                ix_open_peat, latitude, longitude, dem, ts, sfc, vol, Nstems, ba, drain_s, age):   
        #--------------- nutrient contents in vegetation-----------------------
        """
        Computed:
           - total biomass and bottom layer; field layer is gained as a difference of tot and bottom layer (more cohrent results)
           - N and P storage in the each pixel
           - annual use of N and P due to litterfall
        Muukkonen Mäkipää 2005 upland sites: field layer contains dwarf shrubs and (herbs + grasses), see Fig 1
            share     dwarf shrubs     herbs 
            - Pine       91%            9%
            - Spruce     71%            29%
            - broad l    38%            62%
        Peatland sites (assumption):
            share      dwarf shrubs    herbs
            - Pine bogs    90%          10%
            - Spruce mires 50%          50%
        Palviainen et al. 2005 Ecol Res (2005) 20: 652–660, Table 2
        Nutrient concentrations for
                                N              P
            - Dwarf shrubs      1.2%         1.0mg/g
            - herbs & grasses   1.8%         2.0mg/g
            - upland mosses     1.25%        1.4 mg/g
        Nutrient concentrations for sphagna (FIND):
                                N              P     for N :(Bragazza et al Global Change Biology (2005) 11, 106–114, doi: 10.1111/j.1365-2486.2004.00886.x)
            - sphagnum          0.6%           1.4 mg/g     (Palviainen et al 2005)   
        Annual litterfall proportions from above-ground biomass (Mälkönen 1974, Tamm 1953):
            - Dwarf shrubs          0.2
            - herbs & grasses        1
            - mosses                0.3
            Tamm, C.O. 1953. Growth, yield and nutrition in carpets of a forest moss (Hylocomium splendens). Meddelanden från Statens Skogsforsknings Institute 43 (1): 1-140.
        We assume retranslocation of N and P away from senescing tissues before litterfall:
                                N           P
            - Dwarf shrubs     0.5         0.5
            - Herbs & grasses  0.5         0.5
            - mossess          0.0         0.0
        
        Turnover of total biomass including the belowground biomass is assumed to be 1.2 x above-ground biomass turnover
        
        """

        fl_share = {'description': 'share of dwarf shrubs (ds) and herbas & grasses (h) from field layer biomass, kg kg-1',
                    'pine_upland':{'ds': 0.87, 'h': 0.13}, 'spruce_upland':{'ds': 0.71, 'h': 0.29}, 
                    'broadleaved_upland':{'ds': 0.38, 'h': 0.62}, 'spruce_mire':{'ds': 0.90, 'h': 0.10}, 
                    'pine_bog':{'ds': 0.50, 'h': 0.50}}
        nut_con ={'description': 'nutrient concentration of dwarf shrubs (ds), herbs & grasses (h), upland mosses (um), and sphagna (s), unit mg/g',
                  'ds':{'N':12.0, 'P':1.0}, 'h':{'N':18.0, 'P':2.0}, 'um':{'N':12.5, 'P':1.4}, 's':{'N':6.0, 'P':1.4}}
        lit_share = {'description': 'share of living biomass that is lost as litter annually for dwarf shrubs (ds), herbs & grasses (h), upland mosses (um), and sphagna (s), unit: kg kg-1',
                   'ds': 0.2, 'h': 0.5, 'um': 0.3, 's': 0.3}
        retrans ={'description': 'share of nutrients retranslocated before litterfallfor dwarf shrubs (ds), herbs & grasses (h), upland mosses (um), and sphagna (s), unit: kg kg-1',
                  'ds': {'N':0.5, 'P':0.5},'h': {'N':0.5, 'P':0.5}, 
                  'um': {'N':0.0, 'P':0.0},'s': {'N':0.0, 'P':0.7}}
        fl_to_total_turnover = 1.2   # converts the turnover of above-ground bionmass to total including root turnover
        fl_above_to_total = 1.7   # converts aboveground biomass to total biomass 
        
        #--------- create output arrays -----------------------------
        gv_tot = np.zeros(np.shape(gisdata['cmask']))                           # Ground vegetation mass kg ha-1
        gv_field = np.zeros(np.shape(gisdata['cmask']))                         # Field layer vegetation mass
        gv_bot = np.zeros(np.shape(gisdata['cmask']))                           # Bottom layer vegetation mass
        ds_litterfall = np.zeros(np.shape(gisdata['cmask']))                    # dwarf shrub litterfall kg ha-1 yr-1
        h_litterfall = np.zeros(np.shape(gisdata['cmask']))                     # herbs and grasses litterfall kg ha-1 yr-1
        um_litterfall = np.zeros(np.shape(gisdata['cmask']))                    # upland mosses litterfall kg ha-1 yr-1
        s_litterfall = np.zeros(np.shape(gisdata['cmask']))                     # sphagnum mosses litterfall kg ha-1 yr-1
        nup_litter = np.zeros(np.shape(gisdata['cmask']))                       # N uptake due to litterfall kg ha-1 yr-1
        pup_litter = np.zeros(np.shape(gisdata['cmask']))                       # P uptake due to litterfall kg ha-1 yr-1
        n_gv = np.zeros(np.shape(gisdata['cmask']))                             # N in ground vegetation kg ha-1
        p_gv = np.zeros(np.shape(gisdata['cmask']))                             # P in ground vegetation kg ha-1
        zeromask = np.zeros(np.shape(gisdata['cmask']))
        """------ Ground vegetation models from Muukkonen & Mäkipää 2006 BER vol 11, Tables 6,7,8"""    
        #***************** Pine upland ***************************************
        ix = ix_pine_upland
        # dependent variable is sqrt(biomass -0.5), final output in kg ha-1
        gv_tot[ix] = np.square(22.523 + 0.084*sfc[ix]*latitude + 0.01*longitude*age[ix] -0.031*sfc[ix]*age[ix] \
                            -7e-4*np.square(age[ix]) -3e-4*np.square(vol[ix]) +6e-4*vol[ix]* age[ix]) -0.5 +231.56 #Total
        gv_field[ix] = np.square(13.865 +0.013*latitude*longitude -2.969*sfc[ix] +3e-5*ts*age[ix]) - 0.5 + 96.72 #Field layer total. >last term correction factor
        gv_bot[ix] = np.square(8.623 +0.09*sfc[ix]*latitude +0.004*longitude*age[ix] +3e-5*dem[ix]*ts  \
                            -3e-4*np.square(vol[ix]) -5e-4*np.square(age[ix]) + 8e-4*vol[ix]*age[ix] )-0.5 + 355.13 #Bottom layer total
        # removing inconsistent values
        gv_field[ix] = np.minimum(gv_tot[ix], gv_field[ix])
        gv_bot[ix] = np.minimum(gv_tot[ix], gv_bot[ix])        
        gv_field[ix] = np.maximum(gv_field[ix], gv_tot[ix] - gv_bot[ix])
        
        #------------------------------------------------------------------
        #annual litterfall rates and nutrient uptake due to litterfall
        ds_litterfall[ix] = fl_share['pine_upland']['ds']*gv_field[ix]*lit_share['ds']*fl_to_total_turnover
        h_litterfall[ix] = fl_share['pine_upland']['h']*gv_field[ix]*lit_share['h']*fl_to_total_turnover
        um_litterfall[ix] = gv_bot[ix]*lit_share['um']
        n_gv[ix] = gv_field[ix] * fl_share['pine_upland']['ds']*nut_con['ds']['N']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['pine_upland']['h']*nut_con['h']['N']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['um']['N']*1e-3
        p_gv[ix] = gv_field[ix] * fl_share['pine_upland']['ds']*nut_con['ds']['P']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['pine_upland']['h']*nut_con['h']['P']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['um']['P']*1e-3                        
        nup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['N']*1e-3 * (1.0 -retrans['ds']['N']) \
                        +h_litterfall[ix] * nut_con['h']['N']*1e-3 * (1.0 -retrans['h']['N']) \
                        +um_litterfall[ix] * nut_con['um']['N']*1e-3 * (1.0 -retrans['um']['N'])
        pup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['P']*1e-3 * (1.0 -retrans['ds']['P']) \
                        +h_litterfall[ix] * nut_con['h']['P']*1e-3 * (1.0 -retrans['h']['P']) \
                        +um_litterfall[ix] * nut_con['um']['P']*1e-3 * (1.0 -retrans['um']['P'])

        #***************** Spruce upland ***************************************                
        ix = ix_spruce_upland
        gv_tot[ix] = np.square(22.522 +0.026*sfc[ix]*age[ix] +0.11*sfc[ix]*latitude -0.003*sfc[ix]*ts) -0.5 + 206.67  #Total
        gv_field[ix] = np.square(-42.593 + 0.981*latitude -0.008*np.square(ba[ix]) +0.002*ba[ix]*age[ix])-0.5 + 67.15  #Field layer total
        # removing inconsistent values
        gv_field[ix] = np.minimum(gv_tot[ix], gv_field[ix])        
        #gv_bot_ei orig laskennasssa, mutta mosses vastaa 75% gv biomassasta joten hyvä olla mukana
        gv_bot[ix] = np.square(9.672+0.029*sfc[ix]*age[ix]+0.078*sfc[ix]*latitude+0.186*ba[ix]) -0.5 +264.82 #mosses (not bottom layer total available )
        # removing inconsistent values, tarvitaanko näitä enää kun tarkistettu
        gv_bot[ix] = np.minimum(gv_tot[ix], gv_bot[ix]) # taa ei orig laskennassa

        
        
        #annual litterfall rates
        ds_litterfall[ix] = fl_share['spruce_upland']['ds']*gv_field[ix]*lit_share['ds']*fl_to_total_turnover
        h_litterfall[ix] = fl_share['spruce_upland']['h']*gv_field[ix]*lit_share['h']*fl_to_total_turnover
        um_litterfall[ix] = gv_bot[ix]*lit_share['um']
        n_gv[ix] = gv_field[ix] * fl_share['spruce_upland']['ds']*nut_con['ds']['N']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['spruce_upland']['h']*nut_con['h']['N']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['um']['N']*1e-3
        p_gv[ix] = gv_field[ix] * fl_share['spruce_upland']['ds']*nut_con['ds']['P']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['spruce_upland']['h']*nut_con['h']['P']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['um']['P']*1e-3
        nup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['N']*1e-3 * (1.0 -retrans['ds']['N']) \
                        +h_litterfall[ix] * nut_con['h']['N']*1e-3 * (1.0 -retrans['h']['N']) \
                        +um_litterfall[ix] * nut_con['um']['N']*1e-3 * (1.0 -retrans['um']['N'])
        pup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['P']*1e-3 * (1.0 -retrans['ds']['P']) \
                        +h_litterfall[ix] * nut_con['h']['P']*1e-3 * (1.0 -retrans['h']['P']) \
                        +um_litterfall[ix] * nut_con['um']['P']*1e-3 * (1.0 -retrans['um']['P'])
        
        #***************** Broadleaved upland ***************************************        
        ix = ix_broadleaved_upland
        gv_tot[ix] = np.square(19.8 +0.691*sfc[ix]*latitude -38.578*sfc[ix])-0.5 + 156.51      #Total
        gv_field[ix] = np.square(-95.393 +0.094*latitude*longitude -1e-6*Nstems*ts -0.106*longitude**2 +5e-4*latitude*ts)-0.5 + 55.40    #Field layer total
        #alla ei orig laskennssa gv_bot
        gv_bot[ix] = np.square(20.931+0.096*sfc[ix]*latitude-0.0006*longitude*ts)-0.5 + 236.6 #mosses (not bottom layer total available )
        # removing inconsistent values#samoin kuinn spruce upland?
        gv_field[ix] = np.minimum(gv_tot[ix], gv_field[ix])        
        gv_bot[ix] = np.minimum(gv_tot[ix], gv_bot[ix]) # taa ei orig laskennassa

        
        #annual litterfall rates
        ds_litterfall[ix] = fl_share['broadleaved_upland']['ds']*gv_field[ix]*lit_share['ds']*fl_to_total_turnover
        h_litterfall[ix] = fl_share['broadleaved_upland']['h']*gv_field[ix]*lit_share['h']*fl_to_total_turnover
        um_litterfall[ix] = gv_bot[ix]*lit_share['um']
        n_gv[ix] = gv_field[ix] * fl_share['broadleaved_upland']['ds']*nut_con['ds']['N']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['broadleaved_upland']['h']*nut_con['h']['N']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['um']['N']*1e-3
        p_gv[ix] = gv_field[ix] * fl_share['broadleaved_upland']['ds']*nut_con['ds']['P']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['broadleaved_upland']['h']*nut_con['h']['P']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['um']['P']*1e-3
        nup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['N']*1e-3 * (1.0 -retrans['ds']['N']) \
                        +h_litterfall[ix] * nut_con['h']['N']*1e-3 * (1.0 -retrans['h']['N']) \
                        +um_litterfall[ix] * nut_con['um']['N']*1e-3 * (1.0 -retrans['um']['N'])
        pup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['P']*1e-3 * (1.0 -retrans['ds']['P']) \
                        +h_litterfall[ix] * nut_con['h']['P']*1e-3 * (1.0 -retrans['h']['P']) \
                        +um_litterfall[ix] * nut_con['um']['P']*1e-3 * (1.0 -retrans['um']['P'])

     
        #***************** Spruce mire ***************************************
        ix = ix_spruce_mire        
        gv_tot[ix] = np.square(35.52 +0.001*longitude*dem[ix] -1.1*drain_s**2 -2e-5*vol[ix]*Nstems \
                                +4e-5*Nstems*age[ix] +0.139*longitude*drain_s) -0.5 + 116.54 #Total
        gv_bot[ix] =  np.square(-3.182 + 0.022*latitude*longitude +2e-4*dem[ix]*age[ix] \
                                -0.077*sfc[ix]*longitude -0.003*longitude*vol[ix] + 2e-4*np.square(vol[ix]))-0.5 + 98.10  #Bottom layer total
        gv_field[ix] =  np.square(23.24 -1.163*drain_s**2 +1.515*sfc[ix]*drain_s -2e-5*vol[ix]*Nstems\
                                +8e-5*ts*age[ix] +1e-5*Nstems*dem[ix])-0.5 +  162.58   #Field layer total
        # removing inconsistent values
        gv_field[ix] = np.minimum(gv_tot[ix], gv_field[ix])
        gv_bot[ix] = np.minimum(gv_tot[ix], gv_bot[ix])        
        gv_field[ix] = np.maximum(gv_field[ix], gv_tot[ix] - gv_bot[ix])

        #annual litterfall rates
        ds_litterfall[ix] = fl_share['spruce_mire']['ds']*(gv_tot[ix]-gv_bot[ix])*lit_share['ds']*fl_to_total_turnover
        h_litterfall[ix] = fl_share['spruce_mire']['h']*(gv_tot[ix]-gv_bot[ix])*lit_share['h']*fl_to_total_turnover
        s_litterfall[ix] = gv_bot[ix]*lit_share['s']
        n_gv[ix] = gv_field[ix] * fl_share['spruce_mire']['ds']*nut_con['ds']['N']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['spruce_mire']['h']*nut_con['h']['N']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['s']['N']*1e-3
        p_gv[ix] = gv_field[ix] * fl_share['spruce_mire']['ds']*nut_con['ds']['P']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['spruce_mire']['h']*nut_con['h']['P']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['s']['P']*1e-3
        nup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['N']*1e-3 * (1.0 -retrans['ds']['N']) \
                        +h_litterfall[ix] * nut_con['h']['N']*1e-3 * (1.0 -retrans['h']['N']) \
                        +s_litterfall[ix] * nut_con['s']['N']*1e-3 * (1.0 -retrans['s']['N'])
        pup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['P']*1e-3 * (1.0 -retrans['ds']['P']) \
                        +h_litterfall[ix] * nut_con['h']['P']*1e-3 * (1.0 -retrans['h']['P']) \
                        +s_litterfall[ix] * nut_con['s']['P']*1e-3 * (1.0 -retrans['s']['P'])

        
       #***************** Pine bogs ***************************************
        #ix = ix_pine_bog            
        #gv_tot[ix] =  np.square(50.098 +0.005*longitude*dem[ix] -1e-5*vol[ix]*Nstems +0.026*sfc[ix]*age[ix] \
        #            -1e-3*dem[ix]*ts -0.014*vol[ix]*drain_s) - 0.5 + 167.40                #Total           
        #gv_bot[ix] =  np.square(31.809 +0.008*longitude*dem[ix] -3e-4*Nstems*ba[ix] \
        #                        +6e-5*Nstems*age[ix] -0.188*dem[ix]) -0.5 + 222.22                #Bottom layer total
        #gv_field[ix] =  np.square(48.12 -1e-5*ts**2 +0.013*sfc[ix]*age[ix] -0.04*vol[ix]*age[ix] \
        #                        +0.026*sfc[ix]*vol[ix]) - 0.5 +133.26                                        #Field layer total
              #***************** Pine bogs ***************************************
        ix = ix_pine_bog            
        #Nstems=1555. #mean muukkonen makipaa artikkelista		
        #gv_tot[ix] =  np.square(50.098 +0.005*longitude*dem[ix] -1e-5*vol[ix]*Nstems +0.026*sfc[ix]*age[ix] \
        #            -1e-3*dem[ix]*ts -0.014*vol[ix]*drain_s) - 0.5 + 167.40                #Total           
        gv_tot[ix] =  np.square(50.098 +0.005*longitude*dem[ix] -1e-5*vol[ix]*Nstems +0.026*sfc[ix]*age[ix] \
                    -1e-4*dem[ix]*ts -0.014*vol[ix]*drain_s) - 0.5 + 167.40                #Total           
        gv_bot[ix] =  np.square(31.809 +0.008*longitude*dem[ix] -3e-4*Nstems*ba[ix] \
                                +6e-5*Nstems*age[ix] -0.188*dem[ix]) -0.5 + 222.22                #Bottom layer total
        #gv_field[ix] =  np.square(48.12 -1e-5*ts**2 +0.013*sfc[ix]*age[ix] -0.04*vol[ix]*age[ix] \
        #                        +0.026*sfc[ix]*vol[ix]) - 0.5 +133.26                                        #Field layer total
        gv_field[ix] =  np.square(48.12 -1e-5*ts**2 +0.013*sfc[ix]*age[ix] -0.04*vol[ix]*drain_s \
                                +0.026*sfc[ix]*vol[ix]) - 0.5 +133.26                                        #Field layer total
 



        # removing inconsistent values
        gv_field[ix] = np.minimum(gv_tot[ix], gv_field[ix])
        gv_bot[ix] = np.minimum(gv_tot[ix], gv_bot[ix])        
        gv_field[ix] = np.maximum(gv_field[ix], gv_tot[ix] - gv_bot[ix])
              
        #annual litterfall rates
        ds_litterfall[ix] = fl_share['pine_bog']['ds']*(gv_tot[ix]-gv_bot[ix])*lit_share['ds']*fl_to_total_turnover
        h_litterfall[ix] = fl_share['pine_bog']['h']*(gv_tot[ix]-gv_bot[ix])*lit_share['h']*fl_to_total_turnover
        s_litterfall[ix] = gv_bot[ix]*lit_share['s']
        n_gv[ix] = gv_field[ix] * fl_share['pine_bog']['ds']*nut_con['ds']['N']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['pine_bog']['h']*nut_con['h']['N']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['s']['N']*1e-3
        p_gv[ix] = gv_field[ix] * fl_share['pine_bog']['ds']*nut_con['ds']['P']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['pine_bog']['h']*nut_con['h']['P']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['s']['P']*1e-3
        nup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['N']*1e-3 * (1.0 -retrans['ds']['N']) \
                        +h_litterfall[ix] * nut_con['h']['N']*1e-3 * (1.0 -retrans['h']['N']) \
                        +s_litterfall[ix] * nut_con['s']['N']*1e-3 * (1.0 -retrans['s']['N'])
        pup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['P']*1e-3 * (1.0 -retrans['ds']['P']) \
                        +h_litterfall[ix] * nut_con['h']['P']*1e-3 * (1.0 -retrans['h']['P']) \
                        +s_litterfall[ix] * nut_con['s']['P']*1e-3 * (1.0 -retrans['s']['P'])

        #**************** Open peatlands**********************************
        # Not in Mäkipää & Muukkonen, apply Pine bogs
        ix = ix_open_peat            
        age[ix] = 10.
        vol[ix] = 5.
        ba[ix] = 1.
        Nstems=100.

        gv_bot[ix] =  np.square(31.809 +0.008*longitude*dem[ix] -3e-4*Nstems*ba[ix] \
                                +6e-5*Nstems*age[ix] -0.188*dem[ix]) -0.5 + 222.22                #Bottom layer total
        gv_field[ix] =  np.square(48.12 -1e-5*ts**2 +0.013*sfc[ix]*age[ix] -0.04*vol[ix]*age[ix] \
                                +0.026*sfc[ix]*vol[ix]) - 0.5 +133.26                                        #Field layer total
        gv_tot[ix] =  np.square(50.098 +0.005*longitude*dem[ix] -1e-5*vol[ix]*Nstems +0.026*sfc[ix]*age[ix] \
                    -1e-3*dem[ix]*ts -0.014*vol[ix]*drain_s) - 0.5 + 167.40                #Total           
        
        # removing inconsistent values
        gv_field[ix] = np.minimum(gv_tot[ix], gv_field[ix])
        gv_bot[ix] = np.minimum(gv_tot[ix], gv_bot[ix])        
        gv_field[ix] = np.maximum(gv_field[ix], gv_tot[ix] - gv_bot[ix])
        
        #annual litterfall rates
        ds_litterfall[ix] = fl_share['pine_bog']['ds']*(gv_tot[ix]-gv_bot[ix])*lit_share['ds']*fl_to_total_turnover
        h_litterfall[ix] = fl_share['pine_bog']['h']*(gv_tot[ix]-gv_bot[ix])*lit_share['h']*fl_to_total_turnover
        s_litterfall[ix] = gv_bot[ix]*lit_share['s']
        n_gv[ix] = gv_field[ix] * fl_share['pine_bog']['ds']*nut_con['ds']['N']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['pine_bog']['h']*nut_con['h']['N']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['s']['N']*1e-3
        p_gv[ix] = gv_field[ix] * fl_share['pine_bog']['ds']*nut_con['ds']['P']*1e-3*fl_above_to_total \
                        +gv_field[ix] * fl_share['pine_bog']['h']*nut_con['h']['P']*1e-3*fl_above_to_total \
                        +gv_bot[ix] *nut_con['s']['P']*1e-3
        nup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['N']*1e-3 * (1.0 -retrans['ds']['N']) \
                        +h_litterfall[ix] * nut_con['h']['N']*1e-3 * (1.0 -retrans['h']['N']) \
                        +s_litterfall[ix] * nut_con['s']['N']*1e-3 * (1.0 -retrans['s']['N'])
        pup_litter[ix] = ds_litterfall[ix] * nut_con['ds']['P']*1e-3 * (1.0 -retrans['ds']['P']) \
                        +h_litterfall[ix] * nut_con['h']['P']*1e-3 * (1.0 -retrans['h']['P']) \
                        +s_litterfall[ix] * nut_con['s']['P']*1e-3 * (1.0 -retrans['s']['P'])
        
        
        
        #------------Change clear-cut areas: reduce to 1/3 of modelled ---------------------------------------------------
        to_cc = 0.33
        #ix_cc = np.where(np.logical_and(gisdata['age']<5.0, gisdata['smc']!=4))  #small stands excluding open peatlands
        ix_cc = np.where(gisdata['age']<5.0)
        n_gv[ix_cc] = n_gv[ix_cc] * to_cc 
        p_gv[ix_cc] = p_gv[ix_cc] * to_cc
        nup_litter[ix_cc] = nup_litter[ix_cc] * to_cc
        pup_litter[ix_cc] = pup_litter[ix_cc] * to_cc 
        gv_bot[ix_cc] = gv_bot[ix_cc] * to_cc

        if np.nanmin(n_gv) < 0.0:
            print ('NEGATIVE UPTAKE')
            print (np.nanmin(n_gv))
            import sys; sys.exit()
        return n_gv, p_gv, nup_litter, pup_litter, gv_bot

    dem = gisdata['dem'].copy()  # x2 surface elevation m asl
    sfc = np.where(gisdata['sfc']>5, 5, gisdata['sfc'])
    vol = gisdata['vol'].copy()  # x5 stand volume m3 ha-1
    ba = gisdata['ba'].copy()   # x7 basal area m2 ha-1
    age = gisdata['age'].copy()  # x9 stand age in yrs

    # initial N and P mass, kg ha-1    
    n_gv, p_gv, nup_litter, pup_litter, gv_tot = gv_biomass_and_nutrients(gisdata, ix_pine_upland, ix_spruce_upland, ix_broadleaved_upland, ix_spruce_mire, ix_pine_bog,
                ix_open_peat, latitude, longitude, dem, ts, sfc, vol, Nstems, ba, drain_s, age)
    
    # ground vegetation mass at the end of simulation, kg ha-1    
    vol = vol + expected_yield
    age = age + simtime
    n_gv_end, p_gv_end, nup_litter_end, pup_litter_end, gv_tot = gv_biomass_and_nutrients(gisdata, ix_pine_upland, ix_spruce_upland, ix_broadleaved_upland, ix_spruce_mire, ix_pine_bog,
                ix_open_peat, latitude, longitude, dem, ts, sfc, vol, Nstems, ba, drain_s, age)
    
    # nutrient uptake due to net change in gv biomass, only positive values accepted, negative do not associate to nutrient uptake
    nup_net = np.where(n_gv_end - n_gv > 0.0, n_gv_end - n_gv, 0.0)
    pup_net = np.where(p_gv_end - p_gv > 0.0, p_gv_end - p_gv, 0.0)
    
    nup_litter_mean = np.mean([nup_litter, nup_litter_end], axis = 0)
    pup_litter_mean = np.mean([pup_litter, pup_litter_end], axis = 0)
    
    nup = nup_net + nup_litter_mean*simtime         # total N uptake kg ha-1 simulation time (in yrs) -1
    pup = pup_net + pup_litter_mean*simtime         # total P uptake kg ha-1 simulation time (in yrs) -1
    
    popt = False
    if popt:
        import matplotlib.pylab as plt
        fig = plt.figure(num='nutvege')
        plt.subplot(331); plt.imshow(n_gv); plt.colorbar()
        plt.subplot(332); plt.imshow(nup_net); plt.colorbar()
        plt.subplot(333); plt.imshow(nup_litter_mean); plt.colorbar()
        plt.subplot(334); plt.imshow(p_gv); plt.colorbar()
        plt.subplot(335); plt.imshow(pup_net); plt.colorbar()
        plt.subplot(336); plt.imshow(pup_litter); plt.colorbar()
        plt.subplot(337); plt.imshow(gv_tot); plt.colorbar()
        plt.subplot(338); plt.imshow(nup); plt.colorbar()
        plt.subplot(339); plt.imshow(pup); plt.colorbar()
    
    
    return nup, pup


def TAmr(T):
    #Ojanen et al. 2010 Forest Ecology and Management 260:411-421            
    #compute mean growing season air temperature
    import datetime
    import pandas as pd
    dfOut = pd.DataFrame(T)    
    dfOut['doy']= dfOut.index.dayofyear
    start= dfOut.index[0]
    end = dfOut.index[-1]    
    summer = dfOut.groupby(dfOut['doy']).mean()
    start =  datetime.datetime(2010,5,1).timetuple().tm_yday
    end =  datetime.datetime(2010,10,31).timetuple().tm_yday
    TAmr= summer[start:end].mean()
    return TAmr[0]    
    
def peatRespiration(GisData, T, imm_n, imm_p, TASmr=10.4, gwl=None, dt=1):
    """
    Returns potential soil oxygen consumption through soil respiration (CO2) flux through the soil surface according to Ojanen et al 2010 ForEco, eq1
    Boundary conditions: max soil temperature 16 deg C. No minimum value. \n
    Input:
        peatmask, this model for peat only \n
        vol volume of the growung stock m3ha-1 \n       
        Tsoil soil temperature in 5 cm depth, deg C \n
        TAmr mean air temperature May-October, deg C \n        
        gwl ground water level, m, negative down  \n
        bd bulk density of topsoil \n
    Output:
        Heterotrophic respiration, CO2 efflux caused by decomposition: Potential CO2 efflux from soil as kg m-2 s-1 \n
    """
    import numpy as np
    #Ojanen et al. 2010 Forest Ecology and Management 260:411-421            
    
    nrows = int(((GisData['info'][0]).split()[1]))   #shape c,r
    ncols = int(((GisData['info'][1]).split()[1]))
    peatm = GisData['peatm']


    co2 = np.empty((ncols,nrows)); co2[:]=np.nan
    Rref = np.empty((ncols,nrows)); Rref[:]=np.nan
    Nrel = np.empty((ncols,nrows)); Nrel[:]=np.nan
    Prel = Nrel.copy()
    B = Nrel.copy()
    
    #ixp = np.where(peatm==True) 
    #ix_2 = np.where(np.logical_and(np.less_equal(GisData['sfc'],2),np.equal(peatm,True)))   # Rhtkg fertility
    #ix_3 = np.where(np.logical_and(np.equal(GisData['sfc'],3),np.equal(peatm,True)))   # Mtkg fertility
    #ix_4 = np.where(np.logical_and(np.equal(GisData['sfc'],4),np.equal(peatm,True)))   # Ptkg fertility
    #ix_5 = np.where(np.logical_and(np.greater_equal(GisData['sfc'],5),np.equal(peatm,True)))   # Vtkg fertility
    #ix={'sfc_2': ix_2, 'sfc_3': ix_3, 'sfc_4': ix_4, 'sfc_5': ix_5}    
    
    
    ixp = np.greater_equal(GisData['smc'],2)
    #print("smc 6 20 cr:", GisData['smc'][20,6])
    #print(ixp)
    ix_2 = np.where(np.logical_and(np.less_equal(GisData['sfc'],2),np.greater_equal(GisData['smc'],2)))   # Rhtkg fertility
    ix_3 = np.where(np.logical_and(np.equal(GisData['sfc'],3),np.greater_equal(GisData['smc'],2)))   # Mtkg fertility
    ix_4 = np.where(np.logical_and(np.equal(GisData['sfc'],4),np.greater_equal(GisData['smc'],2)))   # Ptkg fertility
    ix_5 = np.where(np.logical_and(np.greater_equal(GisData['sfc'],5),np.greater_equal(GisData['smc'],2)))   # Vtkg fertility
    ix={'sfc_2': ix_2, 'sfc_3': ix_3, 'sfc_4': ix_4, 'sfc_5': ix_5}    
    
    
    
    V = GisData['vol']
    #V=100.
    if T >16. : T=16.
    #if T<-5.: T=-45.    
    T5ref=10.0; T50=-46.02 #; wtm =80.0                                           # wtm mean water table in cm May-Oct, maximum 80 in the Ojanen data, max used becase later cut with the current wt
    pt = 99.0
    bd = {'sfc_2': 0.14, 'sfc_3': 0.11, 'sfc_4': 0.10, 'sfc_5': 0.08}           # Mese study: bulk densities in different fertility classes                                                                 # peat layer thickness, cm            
    N =  {'sfc_2': 1.9, 'sfc_3': 1.6, 'sfc_4': 1.4, 'sfc_5': 1.2}               # Mese study: N cont in OM % dm
    P =  {'sfc_2': 0.1, 'sfc_3': 0.08, 'sfc_4': 0.06, 'sfc_5': 0.05}            # Mese study: P cont in OM % dm

    for k in ix.keys():
        #Rref[ix[k]] = 0.0695 + V[ix[k]]*3.7e-4 + bd[k]*1000.0 * 5.4e-4 +wtm * 1.2e-3       # parameters: Table 3 RHet
        #Rref[ix[k]] = 0.0695 + V*3.7e-4 + bd[k]*1000.0 * 5.4e-4 +wtm * 1.2e-3       # parameters: Table 3 RHet
        Rref[ix[k]] = 0.0695 + V[ix[k]]*3.7e-4 + bd[k]*1000.0 * 5.4e-4 +gwl[ix[k]]*(-100.0) * 1.2e-3       # parameters: Table 3 RHet
    
    #********** Different bd for each fertility class
    for k in ix.keys():
        B[ix[k]]= 156.032 + 16.5*TASmr - 0.196*pt + 0.354*bd[k]*1000.             # Ojanen 2010 et al. Table 4    
    
   #Computes momentary Heterotrophic CO2 flux as a function of soil temperature and peat bulk density    
    for k in ix.keys():
        co2[ix[k]] = Rref[ix[k]]*np.exp(B[ix[k]]*(1.0/(T5ref-T50)-1.0/(T-T50)))                       #g m-2 h-1           
        co2[ix[k]] = co2[ix[k]] *10000. /1000. * 24. * dt                                   # Conversion to kg ha day-1            

    #N release
    C_in_OM = 0.55                                                              # C content in OM kg kg-1
    CO2_to_C = 12./44.
    Nmicrob = imm_n                                                               # microbial immobilisation    
    Pmicrob = imm_p                                                               # microbial immobilisation    
    
    #Restrict CO2 evolution by gwl -> no co2 efflux below the water table    
    #May 21 2020 removed
    #if gwl is not None: co2[ixp]= z_distrib_decomposition(gwl[ixp])*co2[ixp]   

    for k in ix.keys():
        Nrel[ix[k]] = co2[ix[k]] * CO2_to_C / C_in_OM * N[k] / 100. * (1.-Nmicrob)   # Nrelease kg ha-1 day-1
        Prel[ix[k]] = co2[ix[k]] * CO2_to_C / C_in_OM * P[k] / 100. * (1.-Pmicrob)   # Prelease kg ha-1 day-1
    
    return Nrel, Prel, ixp

def z_distrib_decomposition(gwl, Qdistrib_beta = 0.96):
    """
    Distribution Gale & Grigal Canadian Journal of Forest Research, 1987, 17(8): 829-834, 10.1139/x87-131  \n
    Input:
        gwl m, negative or positive (absolute value)
   Output:
        share of actual decomposition from the potential
    """
    qd = Qdistrib_beta**0.0 - Qdistrib_beta**(np.abs(gwl*100.))
    qd=np.where(qd > 0.1, qd, 0.1)
    return qd 
    
def uptake(GisData, motti, simtime):

    nrows = int(((GisData['info'][0]).split()[1]))   #shape c,r
    ncols = int(((GisData['info'][1]).split()[1]))

    #********** Species ************
    vol = GisData['vol']                                                        # total volume in pixel
    #spe = np.zeros(np.shape(GisData['vol']))                                    # initialize species array                                                 
    spe = np.ones(np.shape(GisData['vol']))
    pine = GisData['p_vol']; spruce = GisData['s_vol']; birch = GisData['b_vol']
    
    pidx=np.greater_equal(pine, spruce+birch); sidx=np.less(pine, spruce+birch) # nutrient accumulation in birch is close to birch
    spe[pidx]=1; spe[sidx]=2

    #******* Soil main class **********
    smc = GisData['smc'].copy()
    midx = np.where(GisData['smc']<2); smc[midx] = 1    #mineral
    pidx = np.where(GisData['smc']>1); smc[pidx] = 2    #peat
    oidx = np.where(GisData['smc']==4)  #indices for open peatland

    #******* Soil fertility class, age, stand height
    sfc = GisData['sfc']
    sfc = sfc / GisData['sfc'] *3 # So far all fertility class 3
    age = GisData['age'].copy() #; age[age<1.]= 1.0
    height = GisData['hc'].copy() #; height[height<0.5]=0.5
    
    #**** the soil fertility & main classes present *********
    fC = np.unique(sfc[~np.isnan(sfc)])
    mC = np.unique(smc[~np.isnan(smc)])
    sp = np.unique(spe[~np.isnan(spe)])
    
    #******** creating parameter arrays *****************
    si = np.empty((ncols,nrows)); si[:]=np.nan         #site index
    iage = np.empty((ncols,nrows)); iage[:]=np.nan      #index age
    b1 = np.empty((ncols,nrows)); b1[:]=np.nan          #parameters
    b2 = np.empty((ncols,nrows)); b2[:]=np.nan          

    ysi = np.empty((ncols,nrows)); ysi[:]=np.nan      #parameters for yield
    yiage = np.empty((ncols,nrows)); yiage[:]=np.nan
    yb1 = np.empty((ncols,nrows)); yb1[:]=np.nan
    yb2 = np.empty((ncols,nrows)); yb2[:]=np.nan
    #******* locating to map *************************
    lis = list(product(fC,mC,sp))
    for li in lis:
        s,minpe,spec = li        
        speci = 'Pine' if spec==1 else 'Spruce'
        idx = np.where(np.logical_and(np.logical_and(sfc==s, smc==minpe), spe==spec))
     
        si[idx]=motti[s][minpe][speci]['hpara']['para'][0]
        iage[idx]=motti[s][minpe][speci]['hpara']['para'][1]
        b1[idx]=motti[s][minpe][speci]['hpara']['para'][2]
        b2[idx]=motti[s][minpe][speci]['hpara']['para'][3]

        ysi[idx]=motti[s][minpe][speci]['ypara']['para'][0]
        yiage[idx]=motti[s][minpe][speci]['ypara']['para'][1]
        yb1[idx]=motti[s][minpe][speci]['ypara']['para'][2]
        yb2[idx]=motti[s][minpe][speci]['ypara']['para'][3]

    #********** computing yield for the simulation time ********************        
    f = lambda age, si, iage, b1, b2: si*((1.0-np.exp(-1.0*b1*age))/(1.0-np.exp(-1.0*b1*iage)))**b2 
    hcalc = f(age, si, iage, b1,b2)
    rel = height/hcalc    #relative heigth growth performance    
    rel[rel>2.]=2.
    rel[rel<0.5]=0.5
    del si, iage, b1, b2
    
    #*******  expected yield ************************
    dt = simtime
    y = (f(age+dt, ysi, yiage, yb1, yb2) - f(age, ysi, yiage, yb1, yb2))*rel
    del ysi, yiage, yb1, yb2
    
    
    #N assimilation functions by Raija Laiho 1997 (peat soils)
    nUp = np.empty((ncols,nrows)); nUp[:]=np.nan; pUp = nUp.copy()              # initiate N and P uptake arrays
    idxP = np.where(spe==1)                                                     # indices for pine dominated pixels
    idxS = np.where(spe==2)                                                     # indices for spruce dominated pixels
    #idxB = np.where(spe==3)                                                    # indices for birch dominated pixels   
    
   #********* Nutrient net uptake computation ***********************                                                                   
    MarjoNut = lambda vol, lna, b, k: np.exp(lna + b*np.log(vol) + k)
    par = {                                                                 # Palviainen & Finer, 2012 Eur J For Res 131:945-964, eq 2, Table 7
    'N':{'pine':[1.856,0.631,0.050], 'spruce': [2.864,0.557,0.051], 'birch':[1.590,0.788,0.044]},
    'P':{'pine':[-2.387,0.754,0.158], 'spruce':[-2.112,0.773,0.070], 'birch':[-3.051,1.114,0.102]}
        }
    #******** Nitrogen **************    
    lna,b,k=par['N']['pine']                                                # pine
    nUp[idxP]= MarjoNut(vol[idxP]+y[idxP], lna, b, k) - MarjoNut(vol[idxP], lna, b, k)
    lna,b,k=par['N']['spruce']                                              # spruce
    nUp[idxS]= MarjoNut(vol[idxS]+y[idxS], lna, b, k) - MarjoNut(vol[idxS], lna, b, k)
    #******** Phosphorus **************    
    lna,b,k=par['P']['pine']                                                # pine
    pUp[idxP]= MarjoNut(vol[idxP]+y[idxP], lna, b, k) - MarjoNut(vol[idxP], lna, b, k)
    lna,b,k=par['P']['spruce']                                              # spruce
    pUp[idxS]= MarjoNut(vol[idxS]+y[idxS], lna, b, k) - MarjoNut(vol[idxS], lna, b, k)
 
    #****** from leafarea back to leaf mass
    sla={'pine':5.54, 'spruce': 5.65, 'decid': 18.46}                          # m2/kg, Kellomäki et al. 2001 Atm. Env.    
    llp = {'pine': 3.0, 'spruce':5.0}                                          # leaf longevity, pine, yrs
    lnutcont = {'N':{'pine': 1.0e-2, 'spruce': 1.0e-2},                        # nitrogen content, kg kg-1
                'P': {'pine': 1.0e-3, 'spruce': 1.0e-3}}   
    retrans = {'N':{'pine': 0.5, 'spruce': 0.5},                               # retranslocation
               'P':{'pine': 0.5, 'spruce': 0.5}}
    
    #nUp_gv = 12.0*dt; pUp_gv = 2.0*dt                                          # Palviainen väikkäri n and p uptake by ground vegetation kg/yr
    #nUp_gv = 4.0*dt; pUp_gv = 0.6*dt                                          
    
    # *****N uptake due to changing of leaves ****************
    #ATTN! Multiply with 1e4, not 1e3
    #nleafup = (GisData['LAI_pine'] / sla['pine']*1e3 * lnutcont['N']['pine'] * (1.-retrans['N']['pine']) / llp['pine']  \
    #    + GisData['LAI_spruce'] / sla['spruce']*1e3 * lnutcont['N']['spruce'] * (1.-retrans['N']['spruce']) / llp['spruce']) *dt
    nleafup = (GisData['LAI_pine'] / sla['pine']*1e4 * lnutcont['N']['pine'] * (1.-retrans['N']['pine']) / llp['pine']  \
        + GisData['LAI_spruce'] / sla['spruce']*1e4 * lnutcont['N']['spruce'] * (1.-retrans['N']['spruce']) / llp['spruce']) *dt
    to_gross = 1.0
    Nup_tot = nUp*to_gross + nleafup 

    # *****P uptake due     to changing of leaves *************************    
    #pleafup = (GisData['LAI_pine'] / sla['pine']*1e3 * lnutcont['P']['pine'] * (1.-retrans['P']['pine']) / llp['pine']  \
    #    + GisData['LAI_spruce'] / sla['spruce']*1e3 * lnutcont['P']['spruce'] * (1.-retrans['P']['spruce']) / llp['spruce']) *dt
    pleafup = (GisData['LAI_pine'] / sla['pine']*1e4 * lnutcont['P']['pine'] * (1.-retrans['P']['pine']) / llp['pine']  \
        + GisData['LAI_spruce'] / sla['spruce']*1e4 * lnutcont['P']['spruce'] * (1.-retrans['P']['spruce']) / llp['spruce']) *dt

    to_gross = 1.0    
    Pup_tot = pUp*to_gross + pleafup 
    
    # open peatlands: insert small uptake for stand to avoid empty array 
    Nup_tot[oidx] = 0.5 * simtime # 0.5 kg/yr/ha
    Pup_tot[oidx] = 0.05 * simtime # 0.05 kg/yr/ha
    y[oidx] = 0.1   #m3/ha/yr
    
    return Nup_tot , Pup_tot, y

def nutBalance(sto, up_s, up_gv, rel):
    sto = sto + rel - up_s - up_gv
    sto[sto<0.] = 0.0
    return sto
    

def nConc(cmask, D, Wliq, nsto):
    c=np.empty(np.shape(cmask))
    c[:,:] = np.NaN
    ix = np.isfinite(cmask)
    c[ix] = (nsto[ix]*1e6) / (Wliq[ix] * D[ix] * 1e4 *1e3) #mg/l
    return c
    
def ddgrowth(ddsm, meandd, yrs):
    #distributes nutrient uptake to days, scales with temp sum, 
    #returns fraction on uptake for each day
    dd=np.array(ddsm / (yrs*meandd))
    uprate = np.gradient(dd)
    return uprate
    
def my_Raijan_ravinnef(vol, p1, p2, p3):
    #returns kg/ha        
    return (p1*vol + p2*np.log(vol)**p3)*10.0

def vdataQuery(alue, alku, loppu, kysely, fname=None):
    """
    Runs Vesidata html standard queries    
    
    IN:
        alue - alueid (int)
        alku - '2015-05-25', (str)
        loppu -'2015-06-01', (str)
        kysely -'raw', 'wlevel', 'saa', (str)
        fname - filename for saving ascii-file
    OUT:
        dat - pd DataFrame; index is time and keys are from 1st line of query
    Samuli L. 25.4.2016; queries by Jukka Pöntinen & Anne Lehto
    
    käyttöesim1: https://taimi.in.metla.fi/cgi/bin/12.vesidata_haku.pl?id=3&alku=2016-01-25&loppu=2016-02-10&kysely=wlevel 
    käyttöesim2: https://taimi.in.metla.fi/cgi/bin/12.vesidata_haku.pl?id=Porkkavaara&alku=2016-01-25&loppu=2016-02-10&kysely=raw
    käyttöesim3: https://taimi.in.metla.fi/cgi/bin/12.vesidata_haku.pl?id=Porkkavaara&alku=2016-01-25&loppu=2016-02-10&kysely=saa
    
    vaaditaan parametrit:
    id = Alueen nimi tai numero, esim Kivipuro tai 33, joka on Kivipuron aluenumero, 
         annual-ryhmän kyselyyn voi antaa id=all, jolloin haetaan kaikki alueet
    
    alku = päivä,josta lähtien haetaan 2016-01-25
    
    loppu = päivä,johon saakka haetaan 2016-02-10
    
    kysely: 
    'wlevel' = haetaan vedenkorkeusmuuttujat tietyssä järjestyksessä
    'raw' = haetaan näiden lisäksi kaikki 'raw'-ryhmän muuttujat
    'saa' = haetaan päivittäinen sää, eli sademäärä ja keskilämpötila 
    'annual' = haetaan vuoden lasketut tulokset, päivämäärän alkupäivä on vuoden 1. päivä, loppupäivä esim. vuoden toinen päivä
    'craw'= haetaan kaikki tämän ryhmän muuttujat 
    'dload'= haetaan kaikki tämän ryhmän muuttujat 
    'roff'= haetaan kaikki tämän ryhmän muuttujat 
    'wquality'= haetaan kaikki tämän ryhmän muuttujat 

    """

    import urllib2, os, shutil
    import pandas as pd
    #addr='https://taimi.in.metla.fi/cgi/bin/12.vesidata_haku.pl?id=all&alku=2014-01-01&loppu=2014-10-25&kysely=annual' KAIKKI ANNUAL-MUUTTUJAT
    
    #addr='https://taimi.in.metla.fi/cgi/bin/vesidata_haku.pl?id=Liuhapuro&alku=2015-05-25&loppu=2015-06-10&kysely=raw'
    
    addr='https://taimi.in.metla.fi/cgi/bin/vesidata_haku.pl?id=%s&alku=%s&loppu=%s&kysely=%s' %(str(alue), alku, loppu, kysely)
    ou='tmp.txt'
    
    f=urllib2.urlopen(addr) #open url, read to list and close
    r=f.read().split("\n")
    f.close()
    
    g=open(ou, 'w') #open tmp file, write, close
    g.writelines("%s\n" % item for item in r)
    g.close()
    
    #read  'tmp.txt' back to dataframe
    if kysely is 'annual': #annual query has different format
        dat=pd.read_csv(ou)
        f=dat['v_alue_metodi']
        yr=[]; alue=[]; mtd=[]        
        for k in range(0, len(f)):
            yr.append(float(f[k].split('a')[0]))
            mtd.append(int(f[k].split('d')[1]))
            x=f[k].split('m')[0]
            alue.append(int(x.split('e')[1]))
        dat=dat.drop('v_alue_metodi',1)
        dat.insert(0,'alue_id', alue); dat.insert(1, 'vuosi',yr); dat.insert(2,'mtd', mtd)
        
    else: #...than the other queries
        dat=pd.read_csv(ou,index_col=0)
        dat.index=dat.index.to_datetime() #convert to datetime
    
    
    if kysely is 'wlevel': #manipulate column names
        cols=list(dat.columns.values)
        h=[]  
        for item in cols:
            h.append(item.split("=")[1])
        dat.columns=h
    
    if fname is not None: #copy to fname, remove ou
        shutil.copy(ou, fname)
    os.remove(ou)

    return dat    
