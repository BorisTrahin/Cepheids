'''
Python 3 code
Author: Boris Trahin
Last update: December 27th, 2020
Idea : Transform raw data from Dictionnaire_data to usable data in SPIPS and add some useful informations for SPIPS.
Functions : - readDb_obs : read data from .dpy (created with data_original.py) file for a specified star.
            - writeDb_obs : write data in .dpy file for a specified star.
            - GetData_SPIPS: Function to run. Use dictionnaire_data function to create a directionnary (dic) with data and information.
'''

import pickle
import numpy as np
import math
from numpy import ma
import warnings
from data_original import Dictionnaire_data
import os
import copy

warnings.filterwarnings('ignore')

savedir = "./"
_dBfile_original = './_data_original.dpy'
_dBfile_spips = './_data_spips.dpy'

##############################################################################################################################
def readDb_data(star):
    global _dBfile_original
    if not os.path.exists(_dBfile_original):
        return None
    f = open(_dBfile_original)
    db = pickle.load(f)
    f.close()
    if star in list(db.keys()):
        return db[star]
    else:
        print(list(db.keys()))
        return None

##############################################################################################################################
def writeDb_spips(star, data):
    global _dBfile_spips
    if not os.path.exists(_dBfile_spips):
        f1 = open(_dBfile_spips,'wb')
        f1.close()
        db={star:data}
    else:
        f1 = open(_dBfile_spips,'rb')
        db = pickle.load(f1)
        f1.close()
        db[star] = data
    f2 = open(_dBfile_spips,'wb')
    pickle.dump(db, f2)
    f2.close()
    return

##############################################################################################################################
def GetData_SPIPS(star):
    '''
    Parameters: - star : {str} (Star name : with space (example : FM Aql))
    Returns: - Dictionnary 'dic_i' with two part : one is a list of tuple with data named 'Data' (SPIPS format) and
    the other is 'Info' with different first guess for SPIPS (distance [kpc], color excess, period, Teff, etc.).
    '''
    global _dBfile_original
    global _dBfile_spips
    # Use Dictionnaire_data function to get data in the web database
    dic  = readDb_data(star)
    teff = dic['Data']['Teff']
    diam = dic['Data']['Diameters']
    phot = dic['Data']['Photometry']
    vit = dic['Data']['Radial Velocity']
    info = dic['Info']
    obs = []
    if math.isnan(info['MJD0']):
        info['MJD0'] = 35000.
    if math.isnan(info['Period']):
        info['Period'] = 10.
    
    ###############################################
    # Look for RV data with the reference name and extract data, MJD, uncertainties (format can be different fo all reference, depending of the source)
    ref_vit = list(vit.keys())
    liste_vit = []
    for x in ref_vit:
        if type(vit[x]) == float:
            pass
        else:
            err = None
            if (x == 'Other Data'):
                for i in range(len(vit[x])):
                    obs.append(vit[x][i])
            #References from Vizier
            elif (x == 'Groenewegen+ 2013'):
                MJD = vit[x]['JD'].astype(float)-2400000.5
                vitesse = vit[x]['RV'].astype(float)
                err = np.zeros(len(vitesse)).astype(float) + 0.5
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Barnes+ 2005'):
                MJD = vit[x]['JD'].astype(float)-2400000.5
                vitesse = vit[x]['RV'].astype(float)
                err = np.zeros(len(vitesse)).astype(float) + 0.4
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Storm+ 2011'):
                MJD = vit[x]['JD'].astype(float)-2400000.5
                vitesse = vit[x]['RV'].astype(float)
                err = vit[x]['err_RV'].astype(float)
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Gorynya+ 1992-1998'):
                MJD = vit[x]['JD'].astype(float)-2400000.5
                vitesse = vit[x]['RV'].astype(float)
                err = vit[x]['err_RV'].astype(float)
                legend = ['Reference','MJD','RV','err_RV']
                tab = ['Gorynya+ 1992-1998', MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Petterson+ 2005'):
                legend = list(vit[x].keys())
                ntest=0
                MJD = vit[x]['JD'].astype(float)-2400000.5
                for elt in legend:
                    if ('RV' in elt):
                        vitesse = vit[x]['RV']
                        err = np.zeros(len(vitesse)) + 0.5
                        ntest=1
                    if ntest==1:
                        legend2 = ['Reference','MJD','RV','err_RV']
                        tab = [x, MJD, vitesse,err]
                        dic = dict(list(zip(legend2,tab)))
                        liste_vit.append(dic)
                    else:
                        pass
            elif (x == 'Nardetto+ 2009'):
                phase = vit[x]['Phase'].astype(float)
                vitesse = vit[x]['RV'].astype(float)
                err = vit[x]['err_RV'].astype(float)
                legend = ['Reference','Phase','RV','err_RV']
                tab = [x, phase, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Borgniet+ 2019'):
                MJD = vit[x]['MJD']
                vitesse = vit[x]['RV']
                err = vit[x]['e_RV']
                Inst = vit[x]['Inst']
                legend = ['Reference','MJD','RV','err_RV', 'Inst']
                tab = [x, MJD, vitesse,err, Inst]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Anderson+ 2014') or (x == 'Anderson+ 2016'):
                MJD = vit[x]['JD'].astype(float)-0.5
                vitesse = vit[x]['RV']
                err = vit[x]['err_RV']
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            #References from McMaster
            elif (x == 'Barnes+ 1987') or (x == 'Barnes+ 1988'):
                MJD = vit[x]['HJD-2400000'].astype(float)-0.5
                vitesse = vit[x]['Radial velocity [km/s]'].astype(float)
                err = np.zeros(len(vitesse)) + 3.8
                legend = ['Reference','MJD','RV','err_RV']
                tab = ['Barnes+ 1987-1988', MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Wilson+ 1989'):
                MJD = vit[x]['HJD-2400000'].astype(float)-0.5
                vitesse = vit[x]['Radial velocity [km/s]']
                err = np.zeros(len(vitesse)) + 2.
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Bersier+ 1994'):
                MJD = vit[x]['HJD - 2400000'].astype(float)-0.5
                vitesse = vit[x]['Heliocentric radial velocity [km/s]'].astype(float)
                err = vit[x]['Radial Velocity uncertainty [km/s]'].astype(float)
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Kiss+ 1998'):
                try:
                    MJD = vit[x]['HJD-2400000'].astype(float)-0.5
                    vitesse = vit[x]['Radial velocity v1 [km/sec]'].astype(float)
                    err = np.zeros(len(vitesse)) + 0.2
                    legend = ['Reference','MJD','RV','err_RV']
                    tab = [x, MJD, vitesse,err]
                    dic = dict(list(zip(legend,tab)))
                    liste_vit.append(dic)
                except: pass
            elif (x == 'Lloyd+ 1980'):
                MJD = vit[x]['HJD-2400000'].astype(float)-0.5
                vitesse = vit[x]['Radial Velocity [km/s]'].astype(float)
                err = np.zeros(len(vitesse)) + 2.8
                legend = ['Reference','MJD','RV','err_RV']
                tab = ['Lloyd Evans+ 1980', MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Metzger+ 1992') or (x == 'Metzger+ 1993'):
                MJD = vit[x]['HJD - 2400000'].astype(float)-0.5
                vitesse = vit[x]['Heliocentric radial velocity [km/s]'].astype(float)
                err = vit[x]['Radial Velocity uncertainty [km/s]'].astype(float)
                legend = ['Reference','MJD','RV','err_RV']
                tab = ['Metzger+ 1992-1993', MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Pont+ 1994') or (x == 'Pont+ 1996'):
                MJD = vit[x]['HJD - 2400000'].astype(float)-0.5
                vitesse = vit[x]['Heliocentric radial velocity [km/s]'].astype(float)
                err = vit[x]['Radial Velocity uncertainty [km/s]'].astype(float)
                legend = ['Reference','MJD','RV','err_RV']
                tab = ['Pont+ 1994-1996', MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Coulson+ 1985'):
                MJD = vit[x]['HJD-2400000'].astype(float)-0.5
                if ('Radial velocity [km/s]' in vit[x]):
                    vitesse = vit[x]['Radial velocity [km/s]'].astype(float)
                if ('Radial velocity [km/sec]' in vit[x]):
                    vitesse = vit[x]['Radial velocity [km/sec]'].astype(float)
                err = np.zeros(len(vitesse)) + 1.5
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Gieren+ 1981') or (x == 'Gieren+ 1985') or (x == 'Gieren+ 1989'):
                MJD = vit[x]['HJD-2400000'].astype(float)-0.5
                if ('Radial velocity [km/s]' in vit[x]):
                    vitesse = vit[x]['Radial velocity [km/s]']
                if ('Radial Velocity [km/s]' in vit[x]):
                    vitesse = vit[x]['Radial Velocity [km/s]']
                vitesse2=[]
                a=0
                for i in range(len(vitesse)):
                    if str(vitesse[i])=='-':
                        a=i
                if a!=0:
                    for i in range(len(vitesse)):
                        if i<a:
                            vitesse2.append(float(vitesse[i]))
                        if i==a:
                            vitesse2.append(float(vit[x]['Period: 3.671253 days'][a]))
                            vitesse2[i]=-vitesse2[i]
                        if i>a:
                            vitesse2.append(float(vitesse[i]))
                else:
                    vitesse2=vitesse.astype(float)
                err = np.zeros(len(vitesse2)) + 2.
                legend = ['Reference','MJD','RV','err_RV']
                tab = ['Gieren+ 1981-1994', MJD, vitesse2,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Gieren+ 1994'):
                MJD = vit[x]['HJD - 2400000'].astype(float)-0.5
                vitesse = vit[x]['Heliocentric radial velocity [km/s]'].astype(float)
                err = vit[x]['Radial Velocity uncertainty [km/s]'].astype(float)
                legend = ['Reference','MJD','RV','err_RV']
                tab = ['Gieren+ 1981-1994', MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Metzger+ 1992') or (x == 'Metzger+ 1993') or (x == 'Metzger+ 1991'):
                MJD = vit[x]['HJD - 2400000'].astype(float)-0.5
                vitesse = vit[x]['Heliocentric radial velocity [km/s]'].astype(float)
                err = vit[x]['Radial Velocity uncertainty [km/s]'].astype(float)
                legend = ['Reference','MJD','RV','err_RV']
                tab = ['Metzger+ 1991-1993', MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Evans+ 1993'):
                MJD = vit[x]['HJD-2400000'].astype(float)-0.5
                vitesse = vit[x]['Radial velocity [km/sec] '].astype(float)
                err = np.zeros(len(vitesse)) + 2.
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            elif (x == 'Evans+ 1983'):
                MJD = vit[x]['HJD-2400000 [days]'].astype(float)-0.5
                vitesse = vit[x]['Radial velocity [km/s]'].astype(float)
                err = np.zeros(len(vitesse)) + 2.
                legend = ['Reference','MJD','RV','err_RV']
                tab = [x, MJD, vitesse,err]
                dic = dict(list(zip(legend,tab)))
                liste_vit.append(dic)
            else:
                legend = list(vit[x].keys())
                legend2=[]
                tab=[]
                print(x,vit[x])
                for elt in legend:
                    if ('JD' in elt) and ('2400000' in elt):
                        MJD = vit[x][elt].astype(float)-0.5
                        legend2.append('MJD')
                        tab.append(MJD)
                    if (('velocity' in elt) or ('Velocity' in elt)) and not (('uncertainty' in elt) or ('err' in elt)):
                        RV = vit[x][elt].astype(float)
                        legend2.append('RV')
                        tab.append(RV)   
                    if ('uncertainty' in elt) or ('err' in elt):
                        err = vit[x][elt]
                        legend2.append('err_RV')
                        tab.append(err)
                    if not (('uncertainty' in elt) or ('err' in elt)):
                        err = np.zeros(len(RV)) + 2.
                        legend2.append('err_RV')
                        tab.append(err)
                legend2.append('Reference')
                tab.append(x)
                dic = dict(list(zip(legend2,tab)))
                liste_vit.append(dic)
    for x in liste_vit:
        if ('Nardetto' in x['Reference']):
            for i in range(len(x['Phase'])):
                if ((math.isnan(x['err_RV'][i])) or (math.isnan(x['RV'][i]))):
                    pass
                else:
                    obs.append([x['Phase'][i], 'vrad; '+ x['Reference'], x['RV'][i] , x['err_RV'][i]])
        else:
            for i in range(len(x['MJD'])):
                if ((math.isnan(x['err_RV'][i])) or (math.isnan(x['RV'][i]))):
                    pass
                elif 'Inst' in list(x.keys()):
                    obs.append([x['MJD'][i], 'vrad; '+ x['Reference']+' (%s)'%x['Inst'][i], x['RV'][i] , x['err_RV'][i]])
                else:
                    obs.append([x['MJD'][i], 'vrad; '+ x['Reference'], x['RV'][i] , x['err_RV'][i]])

    ###############################################
    ###############################################
    ###############################################
    # Same for photometry
    ref_phot= list(phot.keys())
    liste_phot = []
    for x in ref_phot:
        if type(phot[x]) == float:
            pass
        else:
            if (x == 'Other Data'):
                for i in range(len(phot[x])):
                    obs.append(phot[x][i])
            #References from Vizier
            elif (x == 'Berdnikov+ 2008'):
                legend = list(phot[x].keys())
                Vmag = phot[x]['V'].astype(float)
                MJD = phot[x]['JD'].astype(float)-2400000.5
                legend2 = ['Reference','MJD','V_GCPD_Johnson']
                tab = [x,MJD,Vmag]
                for elt in legend:
                    if ('B_V' in elt):
                        B_V = phot[x][elt]
                        legend2.append('B_GCPD_Johnson-V_GCPD_Johnson')
                        tab.append(B_V)
                    if ('Rc_Ic' in elt):
                        R_I = phot[x][elt]
                        legend2.append('R_GCPD_Cousins-I_GCPD_Cousins')
                        tab.append(R_I)
                    if ('V_Ic' in elt):
                        V_I = phot[x][elt].astype(float)
                        legend2.append('V_GCPD_Johnson-I_GCPD_Cousins')
                        I = phot[x]['V']-phot[x][elt]
                        legend2.append('I_GCPD_Cousins')
                        tab.append(V_I)
                        tab.append(I)
                    if ('V_Rc' in elt):
                        V_R = phot[x][elt]
                        legend2.append('V_GCPD_Johnson-R_GCPD_Cousins')
                        tab.append(V_R)
                err = np.zeros(len(Vmag)) + 0.02
                legend2.append('err')
                tab.append(err)
                dic = dict(list(zip(legend2,tab)))
                liste_phot.append(dic)
            elif (x =='Berdnikov+ 2014'):
            	pass
            elif (x == 'Feast+ 2008'):
                MJD = phot[x]['JD'].astype(float)-2400000.5
                Jmag = phot[x]['J'] - 0.134*(phot[x]['J']-phot[x]['K']) - 0.001
                Hmag = phot[x]['H'] - 0.022*(phot[x]['J']-phot[x]['K']) + 0.004
                Kmag = phot[x]['K'] - 0.027*(phot[x]['J']-phot[x]['K']) - 0.003
                legend = ['Reference','MJD','J_CTIO_ANDICAM', 'H_CTIO_ANDICAM', 'K_CTIO_ANDICAM','err']
                err = np.zeros(len(MJD)) + 0.08
                tab = [x, MJD, Jmag, Hmag, Kmag, err]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'WISE'):
                MJD = np.array([phot[x]['MJD']])
                W1mag = np.array([phot[x]['W1mag']])
                W2mag = np.array([phot[x]['W2mag']])
                W3mag = np.array([phot[x]['W3mag']])
                W4mag = np.array([phot[x]['W4mag']])
                e_W1mag = np.array([phot[x]['e_W1mag']])
                e_W2mag = np.array([phot[x]['e_W2mag']])
                e_W3mag = np.array([phot[x]['e_W3mag']])
                e_W4mag = np.array([phot[x]['e_W4mag']])
                legend = ['Reference', 'MJD', 'W1_WISE', 'W2_WISE', 'W3_WISE', 'W4_WISE', 'e_W1_WISE', 'e_W2_WISE', 'e_W3_WISE', 'e_W4_WISE']
                tab = [x, MJD, W1mag, W2mag, W3mag, W4mag, e_W1mag, e_W2mag, e_W3mag, e_W4mag]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == '2MASS'):
                MJD = np.array([phot[x]['MJD']])
                Js = np.array([phot[x]['J']])
                Hs = np.array([phot[x]['H']])
                Ks = np.array([phot[x]['Ks']])
                e_Js = np.array([phot[x]['e_J']])
                e_Hs = np.array([phot[x]['e_H']])
                e_Ks = np.array([phot[x]['e_Ks']])
                legend = ['Reference', 'MJD', 'J_2MASS', 'H_2MASS', 'Ks_2MASS', 'e_J_2MASS', 'e_H_2MASS', 'e_Ks_2MASS']
                tab = [x, MJD, Js, Hs, Ks, e_Js, e_Hs, e_Ks]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Gaia eDR3'):
                MJD_G = phot[x]['MJD_G']
                MJD_Gbp = phot[x]['MJD_Gbp']
                MJD_Grp = phot[x]['MJD_Grp']
                G = phot[x]['G']
                Gbp = phot[x]['Gbp']
                Grp = phot[x]['Grp']
                err_G = phot[x]['err_G']
                err_Gbp = phot[x]['err_Gbp']
                err_Grp = phot[x]['err_Grp']
                legend = ['Reference', 'MJD_G', 'MJD_Gbp', 'MJD_Grp', 'G_GAIA_GAI3', 'Gbp_GAIA_GAIA3', 'Grp_GAIA_GAIA3', 'err_G', 'err_Gbp', 'err_Grp']
                tab = [x, MJD_G, MJD_Gbp, MJD_Grp, G, Gbp, Grp, err_G, err_Gbp, err_Grp]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Hipparcos'):
                hip_mag = phot[x]['HP']
                hip_err = phot[x]['err_HP']
                hip_MJD = phot[x]['MJD']
                legend = ['Reference', 'MJD', 'Hp_MvB_Hipparcos', 'err_HP']
                tab = [x, hip_MJD, hip_mag, hip_err]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Monson+ 2011'):
                MJD = phot[x]['JD'].astype(float)-2400000.5
                Jmag = phot[x]['J'].astype(float)
                Kmag = phot[x]['K'].astype(float)
                Hmag = phot[x]['H'].astype(float)
                errH = phot[x]['err_H'].astype(float)
                errJ = phot[x]['err_J'].astype(float)
                errK = phot[x]['err_K'].astype(float)
                legend = ['Reference','MJD','J','H','K', 'err_J', 'err_H', 'err_K']
                tab = [x,MJD,Jmag,Hmag,Kmag,errJ, errH, errK]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x =='Monson+ 2012'):
                mag_3_6 = phot[x]['3.6um'].astype(float)
                mag_4_5 = phot[x]['4.5um'].astype(float)
                err_3_6 = phot[x]['err_3.6'].astype(float)
                err_4_5 = phot[x]['err_4.5'].astype(float)
                MJD = phot[x]['MJD'].astype(float)
                legend = ['Reference', 'MJD', 'I1_Spitzer_IRAC', 'err_3.6', 'I2_Spitzer_IRAC', 'err_4.5']
                tab = [x, MJD, mag_3_6, err_3_6, mag_4_5, err_4_5]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Tycho'):
                tyc_mag_b = phot[x]['B']
                tyc_err_b = phot[x]['err_B']
                tyc_mag_v = phot[x]['V']
                tyc_err_v = phot[x]['err_V']
                tyc_MJD = phot[x]['MJD']
                legend = ['Reference', 'MJD', 'B_MvB_TYCHO', 'err_B', 'V_MvB_TYCHO', 'err_V']
                tab = [x, tyc_MJD, tyc_mag_b, tyc_err_b, tyc_mag_v, tyc_err_v]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x =='AAVSO'):
                legend = list(phot[x].keys())
                legend2=['Reference']
                tab=[x]
                for elt in legend:
                    if (elt == 'MJD_B'):
                        aavso_MJD_b_2 = phot[x]['MJD_B']
                        aavso_MJD_b=[]
                        for i in range(len(aavso_MJD_b_2)):
                            aavso_MJD_b.append(aavso_MJD_b_2[i])
                        legend2.append('MJD_B')
                        tab.append(aavso_MJD_b)
                    if (elt == 'MJD_V'):
                        aavso_MJD_v_2 = phot[x]['MJD_V']
                        aavso_MJD_v=[]
                        for i in range(len(aavso_MJD_v_2)):
                            aavso_MJD_v.append(aavso_MJD_v_2[i])
                        legend2.append('MJD_V')
                        tab.append(aavso_MJD_v)
                    if (elt == 'MJD_R'):
                        aavso_MJD_r_2 = phot[x]['MJD_R']
                        aavso_MJD_r=[]
                        for i in range(len(aavso_MJD_r_2)):
                            aavso_MJD_r.append(aavso_MJD_r_2[i])
                        legend2.append('MJD_R')
                        tab.append(aavso_MJD_r)
                    if (elt == 'MJD_I'):
                        aavso_MJD_i_2 = phot[x]['MJD_I']
                        aavso_MJD_i=[]
                        for i in range(len(aavso_MJD_i_2)):
                            aavso_MJD_i.append(aavso_MJD_i_2[i])
                        legend2.append('MJD_I')
                        tab.append(aavso_MJD_i)
                    if (elt == 'B'):
                        aavso_mag_b = phot[x]['B']
                        legend2.append('B')
                        tab.append(aavso_mag_b)
                    if (elt == 'V'):
                        aavso_mag_v = phot[x]['V']
                        legend2.append('V')
                        tab.append(aavso_mag_v)
                    if (elt == 'R'):
                        aavso_mag_r = phot[x]['R']
                        legend2.append('R')
                        tab.append(aavso_mag_r)
                    if (elt == 'I'):
                        aavso_mag_i = phot[x]['I']
                        legend2.append('I')
                        tab.append(aavso_mag_i)
                    if (elt == 'err_B'):
                        aavso_err_b = phot[x]['err_B']
                        legend2.append('err_B')
                        tab.append(aavso_err_b)
                    if (elt == 'err_V'):
                        aavso_err_v = phot[x]['err_V']
                        legend2.append('err_V')
                        tab.append(aavso_err_v)
                    if (elt == 'err_R'):
                        aavso_err_r = phot[x]['err_R']
                        legend2.append('err_R')
                        tab.append(aavso_err_r)
                    if (elt == 'err_I'):
                        aavso_err_i = phot[x]['err_I']
                        legend2.append('err_I')
                        tab.append(aavso_err_i)
                dic = dict(list(zip(legend2,tab)))
                liste_phot.append(dic)
            #References from McMaster
            elif (x == 'Kiss+ 1998'):
                try:
                    MJD = phot[x]['HJD-2400000'].astype(float)-0.5
                    Vmag = phot[x]['V [mag]'].astype(float)
                    B_V = phot[x]['B-V [mag]'].astype(float)
                    err = np.zeros(len(Vmag)) + 0.02
                    tab = [x, MJD, Vmag, B_V, err]
                    legend = ['Reference', 'MJD', 'V_GCPD_Johnson', 'B_GCPD_Johnson-V_GCPD_Johnson', 'err']
                    dic = dict(list(zip(legend,tab)))
                    liste_phot.append(dic)
                except:
                    MJD = phot[x]['Columns:'].astype(float)-0.5
                    Vmag = phot[x]['         2) phi'].astype(float)
                    B_V = phot[x]['         3) V [mag]'].astype(float)
                    err = np.zeros(len(Vmag)) + 0.02
                    tab = [x, MJD, Vmag, B_V, err]
                    legend = ['Reference', 'MJD', 'V_GCPD_Johnson', 'B_GCPD_Johnson-V_GCPD_Johnson', 'err']
                    dic = dict(list(zip(legend,tab)))
                    liste_phot.append(dic)
            elif (x == 'Coulson+ 1985'):
                legend = list(phot[x].keys())
                Vmag = phot[x]['V (Johnson) [mag]'].astype(float)
                err = np.zeros(len(Vmag)) + 0.02
                MJD = phot[x]['HJD-2400000'].astype(float)-0.5
                legend2 = ['Reference', 'MJD', 'V_GCPD_Johnson', 'err']
                tab = [x, MJD, Vmag, err]
                for elt in legend:
                    if ('B-V' in elt):
                        B_V = phot[x][elt]
                        ntest=0
                        for i in range(len(B_V)):
                            if '*' in str(B_V[i]):
                                ntest=1
                                pass
                        if ntest==0:
                            B_V = phot[x][elt].astype(float)
                            legend2.append('B_GCPD_Johnson-V_GCPD_Johnson')
                            tab.append(B_V)
                    if ('V-I' in elt):
                        V_I = phot[x][elt]
                        ntest=0
                        for i in range(len(V_I)):
                            if '*' in str(V_I[i]):
                                ntest=1
                                pass
                        if ntest==0:
                            V_I = phot[x][elt].astype(float)
                            I = phot[x]['V (Johnson) [mag]'].astype(float)-phot[x][elt].astype(float)
                            legend2.append('V_GCPD_Johnson-I_GCPD_Cousins')
                            legend2.append('I_GCPD_Cousins')
                            tab.append(V_I)
                            tab.append(I)
                    if ('V-R' in elt):
                        V_R = phot[x][elt]
                        ntest=0
                        for i in range(len(V_R)):
                            if '*' in str(V_R[i]):
                                ntest=1
                                pass
                        if ntest==0:
                            V_R = phot[x][elt].astype(float)
                            legend2.append('V_GCPD_Johnson-R_GCPD_Cousins')
                            tab.append(V_R)
                dic = dict(list(zip(legend2,tab)))
                liste_phot.append(dic)
            elif (x == 'Moffett+ 1984'):
                legend = list(phot[x].keys())
                Vmag = phot[x]['V (Johnson) [mag]'].astype(float)
                err = np.zeros(len(Vmag)) + 0.02
                for elt in legend:
                    if (('JD' in elt) or ('DJ' in elt)):
                        MJD = phot[x][elt].astype(float)-0.5
                    if ('B-V' in elt):
                        B_V = phot[x][elt].astype(float)
                legend2 = ['Reference', 'MJD', 'V_GCPD_Johnson', 'B_GCPD_Johnson-V_GCPD_Johnson', 'err']
                tab = [x, MJD, Vmag, B_V, err]
                dic = dict(list(zip(legend2,tab)))
                liste_phot.append(dic)
            elif (x == 'Barnes+ 1997'):
                MJD = phot[x]['HJD-2400000'].astype(float) - 0.5
                Jmag = phot[x]['(J-K) [mag]'].astype(float) + phot[x]['K [mag]'].astype(float)
                Hmag = phot[x]['(H-K) [mag]'].astype(float) + phot[x]['K [mag]'].astype(float)
                Kmag = phot[x]['K [mag]'].astype(float)
                err = np.zeros(len(Jmag)) + 0.02
                legend = ['Reference','MJD','J_CTIO_ANDICAM','H_CTIO_ANDICAM','K_CTIO_ANDICAM','err']
                tab = ['Barnes+ 1997', MJD, Jmag, Hmag, Kmag, err]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Barnes+ 1997_2'):
                legend = list(phot[x].keys())
                Vmag = phot[x]['V [mag]'].astype(float)
                err = np.zeros(len(Vmag)) + 0.02
                for elt in legend:
                    if ('JD' in elt):
                        MJD = phot[x][elt].astype(float)-0.5
                    if ('B-V' in elt):
                        B_V = phot[x][elt].astype(float)
                legend2 = ['Reference','MJD','V_GCPD_Johnson', 'B_GCPD_Johnson-V_GCPD_Johnson', 'err']
                tab = ['Barnes+ 1997', MJD, Vmag, B_V, err]
                dic = dict(list(zip(legend2,tab)))
                liste_phot.append(dic)
            elif (x == 'Pel+ 1976') or (x == 'Walraven+ 1964'):
                MJD = phot[x]['HJD-2400000'].astype(float)-0.5
                V_W = -2.5*phot[x]['Walraven V [mag]'].astype(float)
                V_B_W = phot[x]['Walraven V-B [mag]'].astype(float)
                err = np.zeros(len(V_W)) + 0.02
                legend = ['Reference','MJD','V_Walraven','V_Walraven-B_Walraven','err']
                tab = [x,MJD,V_W,V_B_W,err]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Schechter+ 1992'):
                MJD = phot[x]['HJD-2400000'].astype(float)-0.5
                Kmag = phot[x]['K (CIT) [mag]']
                Jmag = phot[x]['J-K (CIT) [mag]']+Kmag
                Hmag = phot[x]['H-K (CIT) [mag]']+Kmag
                errJ = np.zeros(len(Jmag)) + 0.011
                errH = np.zeros(len(Hmag)) + 0.008
                errK = np.zeros(len(Kmag)) + 0.009
                legend = ['Reference','MJD','J','H','K','err_J','err_H','err_K']
                tab = [x,MJD,Jmag,Hmag,Kmag,errJ,errH,errK]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Laney+ 1992'):
                MJD = phot[x]['HJD-2400000'].astype(float)-0.5
                JmagCarter = phot[x]['J (Carter) [mag]'].astype(float)
                KmagCarter = phot[x]['H (Carter) [mag]'].astype(float)
                KmagCarter = phot[x]['K (Carter) [mag]'].astype(float)
                JmagCIT = jmag - 0.134*(jmag-kmag) - 0.001
                HmagCIT = hmag - 0.022*(jmag-kmag) + 0.004
                KmagCIT = kmag - 0.027*(jmag-kmag) - 0.003
                errJ = np.zeros(len(JmagCIT)) + 0.02
                errH = np.zeros(len(HmagCIT)) + 0.02
                errK = np.zeros(len(KmagCIT)) + 0.02
                legend = ['Reference','MJD','J','H','K','err_J','err_H','err_K']
                tab = [x,MJD,JmagCIT,HmagCIT,KmagCIT,errJ,errH,errK]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Lloyd+ 1980'):
                MJD = phot[x]['HJD-2400000'].astype(float)+2400000 - 2400000.5
                Jmag = phot[x]['J (Glass) [mag]'].astype(float)
                Hmag = phot[x]['H (Glass) [mag]'].astype(float)
                Kmag = phot[x]['K (Glass) [mag]'].astype(float)
                JmagCarter = Jmag - 0.134*(Jmag-Kmag) - 0.001
                HmagCarter = Hmag - 0.022*(Jmag-Kmag) + 0.004
                KmagCarter = Kmag - 0.027*(Jmag-Kmag) - 0.003
                errJ = np.zeros(len(JmagCarter)) + phot[x]['Uncertainty in J [0.01 mag]'].astype(float)*0.01
                errH = np.zeros(len(JmagCarter)) + phot[x]['Uncertainty in H [0.01 mag]'].astype(float)*0.01
                errK = np.zeros(len(JmagCarter)) + phot[x]['Uncertainty in K [0.01 mag]'].astype(float)*0.01
                legend = ['Reference','MJD','J','H','K','err_J','err_H','err_K']
                tab = ['Lloyd Evans+ 1980',MJD,JmagCarter,HmagCarter,KmagCarter,errJ,errH,errK]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Gieren+ 1981') or (x == 'Gieren+ 1985'):
                legend = list(phot[x].keys())
                Vmag = phot[x]['V (Johnson) [mag]'].astype(float)
                err = np.zeros(len(Vmag)) + 0.02
                legend2 = ['Reference', 'V_GCPD', 'err']
                tab = ['Gieren+ 1981-1985', Vmag, err]
                for elt in legend:
                    if ('JD' in elt):
                        MJD = phot[x][elt].astype(float) + 2400000 - 2400000.5
                        legend2.append('MJD')
                        tab.append(MJD)
                    if ('B-V' in elt):
                        B_V = phot[x][elt].astype(float)
                        legend2.append('B_GCPD_Johnson-V_GCPD_Johnson')
                        tab.append(B_V)
                    if ('V-I' in elt):
                        V_I = phot[x][elt].astype(float)
                        I = phot[x]['V (Johnson) [mag]'].astype(float)-phot[x][elt].astype(float)
                        legend2.append('V_GCPD_Johnson-I_GCPD_Cousins')
                        legend2.append('I_GCPD_Cousins')
                        tab.append(V_I)
                        tab.append(I)
                    if ('V-R' in elt):
                        V_R = phot[x][elt].astype(float)
                        legend2.append('V_GCPD_Johnson-R_GCPD_Cousins')
                        tab.append(V_R)
                dic = dict(list(zip(legend2,tab)))
                liste_phot.append(dic)
            elif (x == 'Madore+ 1975'):
                legend = list(phot[x].keys())
                Vmag = phot[x]['V [mag]'].astype(float)
                err = np.zeros(len(Vmag)) + 0.02
                legend2 = ['Reference', 'V_GCPD_Johnson', 'err']
                tab = [x, Vmag, err]
                for elt in legend:
                    if ('JD' in elt):
                        MJD = phot[x][elt].astype(float)-0.5
                    if ('B-V' in elt):
                        B_V = phot[x][elt].astype(float)
                legend = ['Reference','MJD', 'V_GCPD_Johnson', 'B_GCPD_Johnson-V_GCPD_Johnson', 'err']
                tab = [x,MJD, Vmag, B_V, err]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Welch+ 1984'):
                MJD = phot[x]['Geocentric JD - 2400000'].astype(float)+2400000 - 2400000.5
                Jmag = phot[x]['J (CIT) [mag]'].astype(float)
                Hmag = phot[x]['H (CIT) [mag]'].astype(float)
                Kmag = phot[x]['K (CIT) [mag]'].astype(float)
                errJ = np.zeros(len(Jmag)) + 0.02
                errH = np.zeros(len(Hmag)) + 0.02
                errK = np.zeros(len(Kmag)) + 0.02
                legend = ['Reference','MJD','J','H','K','err_J','err_H','err_K']
                tab = [x,MJD,Jmag,Hmag,Kmag,errJ,errH,errK]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Szabados+ 1981') or (x == 'Szabados+ 1980') or (x == 'Szabados+ 1991') or (x == 'Szabados+ 1977'):
                MJD = phot[x]['HJD-2400000'].astype(float)+2400000 - 2400000.5
                B_V = phot[x]['B-V [mag]']
                Vmag = phot[x]['V [mag]']
                for i in range(len(Vmag)):
                    if (':' in str(Vmag[i])):
                        Vmag[i]=Vmag[i].replace(':','')
                for i in range(len(B_V)):
                    if (':' in str(B_V[i])):
                        B_V[i]=B_V[i].replace(':','')
                Vmag = Vmag.astype(float)
                B_V = B_V.astype(float)
                err = np.zeros(len(Vmag)) + 0.02
                legend = ['Reference','MJD', 'V_GCPD_Johnson', 'B_GCPD_Johnson-V_GCPD_Johnson', 'err']
                tab = ['Szabados+ 1977-1991', MJD, Vmag, B_V, err]
                dic = dict(list(zip(legend,tab)))
                liste_phot.append(dic)
            elif (x == 'Harris+ 1980'):
                try:
                    MJD = phot[x]['HJD-2400000'].astype(float)-0.5
                    Vmag = phot[x]['V [mag]'].astype(float)
                    err = np.zeros(len(Vmag)) + 0.02
                    legend = ['Reference','MJD','V_GCPD_Johnson','err']
                    tab = [x,MJD,Vmag,err]
                    dic = dict(list(zip(legend,tab)))
                    liste_phot.append(dic)
                except:
                    pass
            else:
                pass
    for x in liste_phot:
        for band in x:         
            if (band == 'I1_Spitzer_IRAC'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_3.6'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err_3.6'][i]])
            elif (band == 'I2_Spitzer_IRAC'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_4.5'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err_4.5'][i]])
            if (band == 'W1_WISE'):
                for i in range(len(x[band])):
                    if ((math.isnan(x['e_W1_WISE'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['e_W1_WISE'][i]])
            if (band == 'W2_WISE'):
                for i in range(len(x[band])):
                    if ((math.isnan(x['e_W2_WISE'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['e_W2_WISE'][i]])
            if (band == 'W3_WISE'):
                for i in range(len(x[band])):
                    if ((math.isnan(x['e_W3_WISE'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['e_W3_WISE'][i]])
            if (band == 'W4_WISE'):
                for i in range(len(x[band])):
                    if ((math.isnan(x['e_W4_WISE'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['e_W4_WISE'][i]])
            if (band == 'J_2MASS'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['e_J_2MASS'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['e_J_2MASS'][i]])
            if (band == 'H_2MASS'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['e_H_2MASS'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['e_H_2MASS'][i]])
            if (band == 'Ks_2MASS'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['e_Ks_2MASS'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['e_Ks_2MASS'][i]])
            elif (band == 'Hp_MvB_Hipparcos'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_HP'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err_HP'][i]])
            elif (band == 'B_MvB_TYCHO'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_B'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err_B'][i]])
            elif (band == 'V_MvB_TYCHO'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_V'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err_V'][i]])
            elif (band == 'G_GAIA_GAIA2'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_G'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD_G'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err_G'][i]])
            elif (band == 'Gbp_GAIA_GAIA2'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_Gbp'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD_Gbp'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err_Gbp'][i]])
            elif (band == 'Grp_GAIA_GAIA2'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_Grp'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD_Grp'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err_Grp'][i]])
            elif (band == 'B'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_B'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD_B'][i], 'mag; '+ x['Reference'], 'B_GCPD_Johnson', x[band][i] , x['err_B'][i]])
            elif (band == 'V'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_V'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD_V'][i], 'mag; '+ x['Reference'], 'V_GCPD_Johnson', x[band][i] , x['err_V'][i]])
            elif (band == 'R'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_R'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD_R'][i], 'mag; '+ x['Reference'], 'R_GCPD_Cousins', x[band][i] , x['err_R'][i]])
            elif (band == 'I'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_I'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD_I'][i], 'mag; '+ x['Reference'], 'I_GCPD_Cousins', x[band][i] , x['err_I'][i]])
            elif (band == 'B_GCPD_Johnson'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'V_GCPD_Johnson'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'R_GCPD_Cousins'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'I_GCPD_Cousins'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'V_GCPD_Johnson-R_GCPD_Cousins'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'color; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'V_GCPD_Johnson-I_GCPD_Cousins'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'color; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'B_GCPD_Johnson-V_GCPD_Johnson'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'color; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'R_GCPD_Cousins-I_GCPD_Cousins'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'color; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'J'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_J'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], 'J_CTIO_ANDICAM', x[band][i] , x['err_J'][i]])
            elif (band == 'H'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_H'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], 'H_CTIO_ANDICAM', x[band][i] , x['err_H'][i]])
            elif (band == 'K'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err_K'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], 'K_CTIO_ANDICAM', x[band][i] , x['err_K'][i]])
            elif (band == 'J_CTIO_ANDICAM'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'H_CTIO_ANDICAM'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'K_CTIO_ANDICAM'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'V_Walraven'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'mag; '+ x['Reference'], band, x[band][i] , x['err'][i]])
            elif (band == 'B_Walraven-V_Walraven'):
                c = [t for t in x[band] if not math.isnan(t)]
                for i in range(len(x[band])):
                    if ((math.isnan(x['err'][i])) or (math.isnan(x[band][i]))):
                        pass
                    else:
                        obs.append([x['MJD'][i], 'color; '+ x['Reference'], band, -x[band][i] , x['err'][i]])
            else:
                pass

    ###############################################
    ###############################################
    ###############################################
    # Same for Effective Temperature
    ref_teff = list(teff.keys())
    liste_teff = []
    for x in ref_teff:
        if type(teff[x]) == float:
            pass
        else:
            if (x == 'Other Data'):
                for i in range(len(teff[x])):
                    obs.append(teff[x][i])

    ###############################################
    ###############################################
    ###############################################
    # Same for Angular Diameters
    ref_diam = list(diam.keys())
    liste_diam = []
    for x in ref_diam:
        if type(diam[x]) == float:
            pass
        else:
            if (x == 'Other Data'):
                for i in range(len(diam[x])):
                    obs.append(diam[x][i])

    ###############################################
    if len(obs) == 0:
        obs = 'No Data'

    ###############################################
    # Useful information for SPIPS
    P = info['Period']
    PlxGaia = info['Plx eGDR3']
    MJD0 = info['MJD0']
    excess = info['E(B-V)']
    Teff = info['Teff']
    coords = info['Coordinates (RA,DEC) [J2000]']
    info_spips = {'Period' : P, 'MJD0' : MJD0, 'E(B-V)' : excess, 
    'Teff' : Teff, 'Plx eGDR3': PlxGaia, 'Coordinates': coords}

    ###############################################
    # Create dictionnary with observation and parameters 
    dic_i = {}
    dic_i['Data'] = obs
    dic_i['Info'] = info_spips
    # Save in a database
    creat_database = True
    if creat_database:
       writeDb_spips(star, dic_i) 
    return dic_i

##############################################################################################################################