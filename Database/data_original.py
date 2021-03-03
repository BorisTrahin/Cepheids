'''
Python 3 code
Author: Boris Trahin
Last update: December 27th, 2020
Idea : Create a dictionnary with data and information about the star collected with cepheids_data.py.
Functions : - readDb_obs : read data from .dpy file for a specified star.
            - writeDb_obs : write data in .dpy file for a specified star.
            - Dictionnaire_data : Create a dictionnary with data (from McMaster and Vizier) and information about the star (useful in SPIPS). 
'''

import numpy as np
from useful_functions import *
from cepheids_data import *
import warnings
import pickle
import os

warnings.filterwarnings('ignore')
savedir = "./"
_dBfile_original = './_data_original.dpy'

##############################################################################################################################
def readDb_obs(star):
    global _dBfile_original
    if not os.path.exists(_dBfile_original):
        a = False
    else:
    	f = open(_dBfile_original,'rb')
    	db = pickle.load(f)
    	f.close()
    	if star in list(db.keys()):
        	a = True
    	else:
        	a = False
    return a

##############################################################################################################################
def writeDb_obs(star, data):
    global _dBfile_original
    if not os.path.exists(_dBfile_original):
    	f1 = open(_dBfile_original,'wb')
    	f1.close()
    	db={star:data}
    else:
    	f1 = open(_dBfile_original,'rb')
    	db = pickle.load(f1)
    	f1.close()
    	db[star] = data
    f2 = open(_dBfile_original,'wb')
    pickle.dump(db, f2)
    f2.close()
    return

##############################################################################################################################
def Dictionnaire_data(star):
    global _dBfile_original
    tmp = readDb_obs(star)
    if tmp is False:
        print((' > Adding %s to the database (data_original.dpy)'%star))
        tab_final = GetData_McMaster(star)
        dic = {}
        dic_p = {}
        dic_v = {}
        dic_t = {}
        dic_d = {}
        data_p = tab_final['Photometry']
        data_v = tab_final['Radial Velocity']
        #McMaster data photometry
        crocus_p = True
        if (len(data_p) == 0):
            crocus_p = False
        if crocus_p:
            for x in data_p:
                elt = x
                dic2 = {}
                nom_ref = elt[0]
                legend = elt[2]
                data = elt[1]
                legend2 = legend[0:np.shape(data)[1]]
                for i in range(np.shape(data)[1]):
                    dic2[legend2[i]] = data[:,i]
                if not nom_ref in dic_p:
                    dic_p[nom_ref] = dic2
                else:
                    inew=2
                    for i in np.arange(2,10,1):
                        if ('%s_%s'%(nom_ref,i)) in list(dic_p.keys()):
                            inew=i+1
                    dic_p[('%s_%s'%(nom_ref,inew))] = dic2
        #McMaster data RV
        crocus_v = True
        if (len(data_v) == 0):
        	crocus_v = False  
        if crocus_v:                
            for x in data_v:
                elt = x
                dic2 = {}
                nom_ref = elt[0]
                legend = elt[2]
                data = elt[1]
                legend2 = legend[0:np.shape(data)[1]]
                for i in range(np.shape(data)[1]):
                    dic2[legend2[i]] = data[:,i]
                if not nom_ref in dic_v:
                    dic_v[nom_ref] = dic2
                else:
                    inew=2
                    for i in np.arange(2,10,1):
                        if ('%s_%s'%(nom_ref,i)) in list(dic_v.keys()):
                            inew=i+1
                    dic_v[('%s_%s'%(nom_ref,inew))] = dic2
        notavailable=[]
        #AAVSO data photometry
        try:
            AAVSO = GetData_AAVSO(star)
            if AAVSO is not None:
                dic_p.update(AAVSO)
            else:
                notavailable.append('AAVSO')
        except:
            pass
        #Vizier data photometry
        Berdnikov2008 = GetData_Berdnikov2008(star)
        if Berdnikov2008 is not None:
            dic_p.update(Berdnikov2008)
        else:
            notavailable.append('Berdnikov+ 2008')
        Berdnikov2014 = GetData_Berdnikov2014(star)
        if Berdnikov2014 is not None:
            dic_p.update(Berdnikov2014)
        else:
            notavailable.append('Berdnikov+ 2014')
        Feast2008 = GetData_Feast2008(star)
        if Feast2008 is not None:
            dic_p.update(Feast2008)
        else:
            notavailable.append('Feast+ 2008')
        Hipparcos = GetData_Hipparcos(star)
        if Hipparcos is not None:
            dic_p.update(Hipparcos)
        else:
            notavailable.append('Hipparcos')
        GaiaDR2 = GetData_GaiaDR2(star)
        if GaiaDR2 is not None:
            dic_p.update(GaiaDR2)
        else:
            notavailable.append('Gaia DR2')
        Monson2011 = GetData_Monson2011(star)
        if Monson2011 is not None:
            dic_p.update(Monson2011)
        else:
            notavailable.append('Monson+ 2011')
        Monson2012 = GetData_Monson2012(star)
        if Monson2012 is not None:
            dic_p.update(Monson2012)
        else:
            notavailable.append('Monson+ 2012')
        MASS = GetData_2MASS(star)
        if MASS is not None:
            dic_p.update(MASS)
        else:
            notavailable.append('2MASS')
        Tycho = GetData_Tycho(star)
        if Tycho is not None:
            dic_p.update(Tycho)
        WISE = GetData_WISE(star)
        if WISE is not None:
            dic_p.update(WISE)
        else:
            notavailable.append('WISE')
        #Vizier data RV
        Anderson2014 = GetData_Anderson2014(star)
        if Anderson2014 is not None:
            dic_v.update(Anderson2014)
        else:
            notavailable.append('Anderson+ 2014')
        Anderson2016 = GetData_Anderson2016(star)
        if Anderson2016 is not None:
            dic_v.update(Anderson2016)
        else:
            notavailable.append('Anderson+ 2016')
        Barnes2005 = GetData_Barnes2005(star)
        if Barnes2005 is not None:
            dic_v.update(Barnes2005)
        else:
            notavailable.append('Barnes+ 2005')
        Borgniet2019 = GetData_Borgniet2019(star)
        if Borgniet2019 is not None:
            dic_v.update(Borgniet2019)
        else:
            notavailable.append('Borgniet+ 2019')
        Gorynya1992 = GetData_Gorynya1992(star)
        if Gorynya1992 is not None:
            dic_v.update(Gorynya1992)
        else:
            notavailable.append('Gorynya+ 1992-1998')
        Groenewegen2013 = GetData_Groenewegen2013(star)
        if Groenewegen2013 is not None:
            dic_v.update(Groenewegen2013)
        else:
            notavailable.append('Groenewegen+ 2013')
        Nardetto2009 = GetData_Nardetto2009(star)
        if Nardetto2009 is not None:
            dic_v.update(Nardetto2009)
        else:
            notavailable.append('Nardetto+ 2009')
        Petterson2005 = GetData_Petterson2005(star)
        if Petterson2005 is not None:
            dic_v.update(Petterson2005)    
        else:
            notavailable.append('Petterson+ 2005')  
        Storm2011 = GetData_Storm2011(star)
        if Storm2011 is not None:
            dic_v.update(Storm2011)
        else:
            notavailable.append('Storm+ 2011')
        Storm2004 = GetData_Storm2004(star)
        if Storm2004 is not None:
            dic_v.update(Storm2004)
        else:
            notavailable.append('Storm+ 2004')
        # Other Data
        OtherData = GetData_OtherData(star)
        if OtherData is not None:
            if OtherData[0]=={}:
                dic_t['Teff'] = 'No data'
            else:
                dic_t.update(OtherData[0])
            if OtherData[1]=={}:
                dic_d['Diameters'] = 'No data'
            else:
                dic_d.update(OtherData[1])
            if OtherData[2]=={}:
                pass
            else:
                dic_v.update(OtherData[2])
            if OtherData[3]=={}:
                pass
            else:
                dic_p.update(OtherData[3])
        else:
            notavailable.append('Other Data')
        dic_data = {}
        dic['Photometry'] = dic_p
        dic['Radial Velocity'] = dic_v
        dic['Diameters'] = dic_d
        dic['Teff'] = dic_t
        info = info_cepheids(star)
        print('      --- DONE ---')
        dic_data['Info'] = info
        dic_data['Data'] = dic
        i, j = 0, 0
        for x in dic_p:
            if type(dic_p[x]) == float:
                i += 1
        for x in dic_v:
            if type(dic_v[x]) == float:
                j += 1
        if ((len(dic_p) == i) & (len(dic_v) == j)):
            print(('No data found for %s'%star))
        create_database = True        
        if create_database:
            writeDb_obs(star, dic_data)
    else:
        print((' > %s already in the database (data_original.dpy)'%star))
        f = open(_dBfile_original)
        db = pickle.load(f)
        f.close()
        for i in references:
            if i in list(db[star]['Data']['Photometry'].keys()) or i in list(db[star]['Data']['Radial Velocity'].keys()):
                pass
            else:
                pass
        dic_data = db[star]
    return dic_data

##############################################################################################################################