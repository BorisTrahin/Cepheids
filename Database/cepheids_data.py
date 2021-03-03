'''
Python 3 code
Author: Boris Trahin
Last update: December 27th, 2020
Idea : Get data some online databases as Vizier, Crocus and our personnal database
Functions : - coordinates: Get the coordinates of a star in SkyCoord form.
			- GetData_Anderson2014 : Get Anderson+ 2014 (J/A+A/566/L10, RV) data from Vizier for the specified star.
			- GetData_Anderson2016 : Get Anderson+ 2016 (J/ApJS/226/18, RV) data from Vizier for the specified star.
			- GetData_Barnes2005 : Get Barnes+ 2005 (J/ApJS/156/227, RV) data from Vizier for the specified star.
			- GetData_Berdnikov2008 : Get Berdnikov+ 2008 (II/285, V,B-V,V-R,R-I,V-I) data from Vizier for the specified star.
			- GetData_Berdnikov2014 : Get Berdnikov+ 2014 (J/PAZh/40/147, V,B,I) data from Vizier for the specified star.
			- GetData_Borgniet2019: Get Borgniet+ 2019 (J/A+A/631/A37, RV) data from Vizier for the specified star.
			- GetData_Feast2008 : Get Feast+ 2008 (J/MNRAS/386/2115, J,H,K,L,V-K) data from Vizier for the specified star.
			- GetData_GaiaDR2 : Get Gaia DR2 (ESA Archive, G, Gbp, Grp) data for the specified star.
			- GetData_GaiaeDR3 : Get Gaia eDR3 (I/350/gaiaedr3) parameters from Vizier for the specified star.
			- GetData_Gorynya1992 : Get Gorynya+ 1992-1998 (III/229, RV) data from Vizier for the specified star.
			- GetData_Groenewegen2013 : Get Groenewegen+ 2013 (J/A+A/550/A70, RV) data from Vizier for the specified star.
			- GetData_Hipparcos : Get Hipparcos (ESA Archive, HP) data for the specified star.
			- GetData_Monson2011 : Get Monson+ 2011 (J/ApJS/193/12/table4, J,H,K) data from Vizier for the specified star.
			- GetData_Monson2012 : Get Monson+ 2012 (J/ApJ/759/146, Spitzer I1 and I2) data from Vizier for the specified star.
			- GetData_Nardetto2009 : Get Nardetto+ 2009 (J/A+A/502/951, RV) data from Vizier for the specified star (HARPS data of 8 galactic cepheids with cross-correlation).
			- GetData_Petterson2005 : Get Petterson+ 2005 (J/MNRAS/362/1167, RV) data from Vizier for the specified star.
			- GetData_Storm2004 : Get Storm+ 2004 (J/A+A/415/531/, RV) data from Vizier for the specified star.
			- GetData_Storm2011 : Get Storm+ 2011 (J/A+A/534/A94, RV) data from Vizier for the specified star.
			- GetData_Tycho : Get Tycho (ESA Archive, B,V) data for the specified star.
			- GetData_WISE : Get WISE (II/311/wise, W1, W2, W3, W4) data from Vizier for the specified star.
			- GetData_2MASS : Get 'II/246/out (II/311/wise, J, H, Ks) data from Vizier for the specified star.
			- GetData_AAVSO : Get AAVSO WebObs data for the specified star.
			- Reference_McMaster : Create 2 arrays containing references of photometry and radial velocity available in McMaster for the specified star.
			- GetData_McMaster : Get McMaster data available in the references obtained with Reference_McMaster for the selected star.
			- GetData_OtherData : Other data from own observations or not in the electronic form.
			- GetData_Groenewegen2018 : Get Groenewegen+ 2018 (J/A+A/619/A8/) mean magnitudes in the 2MASS system from Vizier for the specified star.
			- info_cepheids : Get some informations of the selected star (period, distance, ...).
'''

from numpy import ma
from astropy.io import ascii
from astropy.table import Table
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
from OtherData import OtherData
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import SkyCoord
from astropy.utils.data import conf
import astropy.units as u
import pickle
from bs4 import BeautifulSoup
import urllib.request, urllib.parse, urllib.error
from urllib.request import urlopen
import requests
import numpy as np
import math
import os
from useful_functions import *
import sys
import importlib
importlib.reload(sys)
sys.setdefaultencoding('utf-8')

def coordinates(star):
	new_limit = 60
	with conf.set_temp('remote_timeout', new_limit):
		coord = SkyCoord.from_name(star)
	return coord

references = ['Anderson+ 2014', 'Anderson+ 2016', 'Barnes+ 2015', 'Berdnikov+ 2008', 'Berdnikov+ 2014', 'Borgniet+ 2019', 'Feast+ 2008', 'Gaia DR2', 'Gorynya+ 1992-1998', 'Groenewegen+ 2013', 'Hipparcos', 'Monson+ 2011', 'Monson+ 2012', 'Nardetto+ 2009', 'Petterson+ 2005', 'Storm+ 2011', 'Tycho', 'WISE', '2MASS', 'AAVSO']
##############################################################################################################################
#####################################################      VIZIER      #######################################################
##############################################################################################################################
def GetData_Anderson2014(star):
	v = Vizier(columns=["Name"])
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/A+A/566/L10')
	d={}
	if str(star2)=='Empty TableList':
		return
	else:
		print('      - downloading Anderson+ 2014 data from Vizier...')
		v = Vizier(catalog=["J/A+A/566/L10"])
		v.ROW_LIMIT = -1
		result = v.get_catalogs('J/A+A/566/L10')
		if star2[0][0][0] == 'RS Pup':
			JD = ma.getdata(result[4]['BJD'])
			RV = ma.getdata(result[4]['RV'])
			err_RV = ma.getdata(result[4]['e_RV'])
			legend = ['JD','RV','err_RV']
			row_table = [JD, RV, err_RV]
			d['Anderson+ 2014'] = dict(list(zip(legend,row_table)))      
		if star2[0][0][0] == 'l Car':
			JD = ma.getdata(result[3]['BJD'])
			RV = ma.getdata(result[3]['RV'])
			err_RV = ma.getdata(result[3]['e_RV'])
			legend = ['JD','RV','err_RV']
			row_table = [JD, RV,err_RV]
			d['Anderson+ 2014'] = dict(list(zip(legend,row_table)))       
		if star2[0][0][0] == 'QZ Nor':
			JD = ma.getdata(result[1]['BJD'])
			RV = ma.getdata(result[1]['RV'])
			err_RV = ma.getdata(result[1]['e_RV'])
			legend = ['JD','RV','err_RV']
			row_table = [JD, RV,err_RV]
			d['Anderson+ 2014'] = dict(list(zip(legend,row_table)))
		if star2[0][0][0] == 'V335 Pup':
			JD = ma.getdata(result[2]['BJD'])
			RV = ma.getdata(result[2]['RV'])
			err_RV = ma.getdata(result[2]['e_RV'])
			legend = ['JD','RV','err_RV']
			row_table = [JD, RV,err_RV]
			d['Anderson+ 2014'] = dict(list(zip(legend,row_table)))
		return d

##############################################################################################################################
def GetData_Anderson2016(star):
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/ApJS/226/18')
	d={}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/ApJS/226/18/Table2',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Anderson+ 2016 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/ApJS/226/18/Table2',Name=star2)[0]
			JD = ma.getdata(result['BJD'])
			RV = ma.getdata(result['RVel'])
			err_RV = ma.getdata(result['e_RVel'])
			legend = ['JD','RV','err_RV']
			row_table = [JD, RV, err_RV]
			d['Anderson+ 2016'] = dict(list(zip(legend,row_table)))      
		return d

##############################################################################################################################
def GetData_Barnes2005(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/ApJS/156/227')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/ApJS/156/227/table3',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Barnes+ 2005 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/ApJS/156/227/table3',Name=star2)[0]
			JD = ma.getdata(result['HJD'])
			E = ma.getdata(result['E'])
			RV = ma.getdata(result['RV'])
			Phase = ma.getdata(result['Phase'])
			legend = ['JD','E','Phase','RV']
			row_table = [JD, E, Phase, RV]
			d['Barnes+ 2005'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Berdnikov2008(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='II/285')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='II/285/photo',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Berdnikov+ 2008 data from Vizier...')
			result = Vizier.query_constraints(catalog='II/285/photo',Name=star2)[0]
			JD = ma.getdata(result['JD'])
			Vmag = ma.getdata(result['Vmag'])
			B_V = ma.getdata(result['B-V'])
			V_Rc = ma.getdata(result['V-Rc'])
			Rc_Ic = ma.getdata(result['Rc-Ic'])
			V_Ic = ma.getdata(result['V-Ic'])

			B_V2 = np.array([x for x in B_V if not math.isnan(x)])
			V_Rc2 = np.array([x for x in V_Rc if not math.isnan(x)])
			Rc_Ic2 = np.array([x for x in Rc_Ic if not math.isnan(x)]) 
			V_Ic2 = np.array([x for x in V_Ic if not math.isnan(x)]) 
			legend = ['JD','V']
			row_table = [JD,Vmag]
			if len(B_V2) != 0:
				legend.append('B_V')
				row_table.append(B_V)
			if len(V_Rc2) != 0:
				legend.append('V_Rc')
				row_table.append(V_Rc)
			if len(Rc_Ic2) != 0:
				legend.append('Rc_Ic')
				row_table.append(Rc_Ic)
			if len(V_Ic2) != 0:
				legend.append('V_Ic')
				row_table.append(V_Ic)
			d['Berdnikov+ 2008'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Berdnikov2014(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/PAZh/40/147')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/PAZh/40/147/table1',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Berdnikov+ 2014 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/PAZh/40/147/table1',Name=star2)[0]
			filt = ma.getdata(result['Filt'])
			JD = ma.getdata(result['HJD'])
			Vmag = ma.getdata(result['mag'][filt == 'V'])
			Bmag = ma.getdata(result['mag'][filt == 'B'])
			Icmag = ma.getdata(result['mag'][filt == 'Ic'])
			legend = ['JD','V','B','Ic']
			row_table = [JD, Vmag, Bmag, Icmag]
			d['Berdnikov+ 2014'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Borgniet2019(star):    
	Vizier.ROW_LIMIT = -1
	result = Vizier.query_region(coordinates(star), radius="2s",catalog='J/A+A/631/A37/')
	d = {}
	if str(result)=='Empty TableList':
		return
	else:
		star=result[0][0][0]
		MJD, RV, e_RV, Inst = [], [], [], []
		result_b = Vizier.query_constraints(catalog='J/A+A/631/A37/rv-bmc',Star=star)
		result_g = Vizier.query_constraints(catalog='J/A+A/631/A37/rv-gac',Star=star)
		result_r = Vizier.query_constraints(catalog='J/A+A/631/A37/rv-rmc',Star=star)
		if str(result_b)=='Empty TableList':
			pass
		else:
			for i in range(len(ma.getdata(result_b[0]['MJD']))):
				RV.append(ma.getdata(result_b[0]['rvc'])[i])
				MJD.append(ma.getdata(result_b[0]['MJD'])[i])
				e_RV.append(ma.getdata(result_b[0]['e_rvc'])[i])
				Inst.append(result_b[0]['Inst'][i])
		if str(result_g)=='Empty TableList':
			pass
		else:
			print('      - downloading Borgniet+ 2019 data from Vizier...')
			for i in range(len(ma.getdata(result_g[0]['MJD']))):
				MJD.append(ma.getdata(result_g[0]['MJD'])[i])
				RV.append(ma.getdata(result_g[0]['rvc'])[i])
				e_RV.append(ma.getdata(result_g[0]['e_rvc'])[i])
				Inst.append(result_g[0]['Inst'][i])
		if str(result_r)=='Empty TableList':
			pass
		else:
			for i in range(len(ma.getdata(result_r[0]['MJD']))):
				MJD.append(ma.getdata(result_r[0]['MJD'])[i])
				RV.append(ma.getdata(result_r[0]['rvc'])[i])
				e_RV.append(ma.getdata(result_r[0]['e_rvc'])[i])
				Inst.append(result_r[0]['Inst'][i])
		if len(MJD)==0:
			return
		else:
			legend = ['MJD','RV','e_RV', 'Inst']
			row_table = [MJD, RV, e_RV, Inst]
			d['Borgniet+ 2019'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Feast2008(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/MNRAS/386/2115')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/MNRAS/386/2115/tablea1',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Feast+ 2008 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/MNRAS/386/2115/tablea1',Name=star2)[0]
			JD = ma.getdata(result['HJD'])
			Kmag = ma.getdata(result['Kmag'])
			Jmag = ma.getdata(result['Jmag'])
			Hmag = ma.getdata(result['Hmag'])
			V_K = ma.getdata(result['V-K'])
			Lmag = ma.getdata(result['Lmag'])
			Phase = ma.getdata(result['Phase'])
			legend = ['JD','Phase', 'J', 'H','K','V_K','L','err_mag']
			row_table = [JD, Phase, Jmag, Hmag, Kmag, V_K, Lmag, 0.008]
			d['Feast+ 2008'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_GaiaDR2(star): 
	data = {'Gaia DR2':[]}
	s = Simbad.query_object(star)
	width = u.Quantity(1, u.arcsec)
	height = u.Quantity(1, u.arcsec)
	r = Gaia.query_object_async(coordinate=coordinates(star), width=width, height=height)
	r['source_id', 'ra', 'dec', 'phot_g_mean_mag']
	id = r['source_id']
	url = 'http://geadata.esac.esa.int/data-server/data?RETRIEVAL_TYPE=epoch_photometry&ID=%d&VALID_DATA=false&FORMAT=CSV'
	try:
		f = urlopen(url%(int(id)))
	except:
		return
	print('      - downloading Gaia DR2 data from ESA Archive...')
	for l in f.readlines():
		if not 'mag' in list(data.keys()):
			cols =  l.split(',')
			data.update({c:[] for c in cols})
		elif len(l)>10 and l.split(',')[3]!='':
			for i,v in enumerate(l.split(',')):
				try:
					v = float(v)
				except:
					pass
				data[cols[i]].append(v)
			# -- make SPIPS data points:
			f = {'G':'G_GAIA_GAIA2',
				'BP':'Gbp_GAIA_GAIA2',
				'RP':'Grp_GAIA_GAIA2'}
			data['Gaia DR2'].append([data['time'][-1]+55197.0, f[data['band'][-1]], data['mag'][-1], 2.5/np.log(10)*1/data['flux_over_error'][-1]])
		else:
			pass
	d={}
	tmp = [s for s in data['Gaia DR2']]
	obs=[]
	for t in tmp:
		t[-1]=float(t[-1]*5)
	obs.extend(tmp)
	MJD_G, MJD_Gbp, MJD_Grp, G, Gbp, Grp, err_G, err_Gbp, err_Grp= [], [], [], [], [], [], [], [], []
	for i in range(len(obs)):
		if obs[i][1]=='G_GAIA_GAIA2':
			MJD_G.append(obs[i][0])
			G.append(obs[i][2]) 
			err_G.append(obs[i][3])
		if obs[i][1]=='Gbp_GAIA_GAIA2':
			MJD_Gbp.append(obs[i][0])
			Gbp.append(obs[i][2]) 
			err_Gbp.append(obs[i][3])
		if obs[i][1]=='Grp_GAIA_GAIA2':
			MJD_Grp.append(obs[i][0])
			Grp.append(obs[i][2]) 
			err_Grp.append(obs[i][3])
	legend = ['MJD_G', 'MJD_Gbp', 'MJD_Grp', 'G', 'Gbp', 'Grp', 'err_G', 'err_Gbp', 'err_Grp']
	row_table = [MJD_G, MJD_Gbp, MJD_Grp, G, Gbp, Grp, err_G, err_Gbp, err_Grp]
	if len(MJD_G)==0:
		return
	else:
		d['Gaia DR2'] = dict(list(zip(legend,row_table)))
		return d

##############################################################################################################################
def GetData_eGDR3(star):    
	Vizier.ROW_LIMIT = -1
	result = Vizier.query_region(coordinates(star), radius="1s",catalog='I/350/gaiaedr3')
	d = {}
	if str(result)=='Empty TableList':
		plx = 'Not in eGDR3'
		e_plx = 'Not in eGDR3'
	else:
		plx = ma.getdata(result[0][0]['Plx'])
		e_plx = ma.getdata(result[0][0]['e_Plx'])
	legend = ['Plx eGDR3']
	row_table = [[plx,e_plx]]
	d = dict(list(zip(legend,row_table)))
	return d

##############################################################################################################################
def GetData_Gorynya1992(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='III/229')
	d = {}
	if str(star2)=='Empty TableList':
		star2 = v.query_region(coordinates(star), radius="10m",catalog='III/229')
		if str(star2)=='Empty TableList':
			return
		else:
			star2=star2[0][0][0]
			result = Vizier.query_constraints(catalog='III/229/catalog',Name=star2)
			if str(result)=='Empty TableList':
				return
			else:
				print('      - downloading Gorynya+ 1992-1998 data from Vizier...')
				result = Vizier.query_constraints(catalog='III/229/catalog',Name=star2)[0]
				JD = ma.getdata(result['JD'])
				RV = ma.getdata(result['HRV'])
				err_RV = ma.getdata(result['e_HRV'])
				legend = ['JD','RV','err_RV']
				row_table = [JD, RV, err_RV]
				d['Gorynya+ 1992-1998'] = dict(list(zip(legend,row_table)))
				return d
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='III/229/catalog',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Gorynya+ 1992-1998 data from Vizier...')
			result = Vizier.query_constraints(catalog='III/229/catalog',Name=star2)[0]
			JD = ma.getdata(result['JD'])
			RV = ma.getdata(result['HRV'])
			err_RV = ma.getdata(result['e_HRV'])
			legend = ['JD','RV','err_RV']
			row_table = [JD, RV, err_RV]
			d['Gorynya+ 1992-1998'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Groenewegen2013(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/A+A/550/A70')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/A+A/550/A70/table3',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Groenewegen+ 2013 data from Vizier...')
			result = result[0]
			JD = ma.getdata(result['JD'])
			RV = ma.getdata(result['RV'])
			legend = ['JD', 'RV']
			row_table = [JD, RV]
			d['Groenewegen+ 2013'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Hipparcos(star):
	try:
		objet = Simbad.query_objectids(star)
	except:
		objet = Simbad.query_objectids(star)
	obj = list(objet['ID'])
	hip_name = np.nan
	n = 0
	for elt in obj:
		if 'HIP' in elt:
			hip_name = elt.split()[1]
			n = 1
	if hip_name is not np.nan :
		print('      - downloading Hipparcos data from Vizier...')
		html_HIPPARCOS = 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Plot/Vgraph/htm?I/239/' + hip_name
		try:
			hip_html = requests.get(html_HIPPARCOS)
			hip_soup = BeautifulSoup(hip_html.text,"html5lib").text.split('\n')
		except:
			hip_html = requests.get(html_HIPPARCOS)
			hip_soup = BeautifulSoup(hip_html.text,"html5lib").text.split('\n')
		hip_n_init = 0
		for x in hip_soup:
			if 'JD-2440000' in x:
				hip_n_init = hip_soup.index(x)
	
		hip_data = []
		for i in np.arange(hip_n_init+1,len(hip_soup)-2,1):
			hip_data.append(fix_unicode(hip_soup[i].split('|')))
		for i in range(len(hip_data)):
			for j in range(len(hip_data[i])):
				hip_data[i][j]=hip_data[i][j].replace(' ','')
		hip_data = [item for item in hip_data if item[3]=='0']
		hip_MJD = [float(item[0])-0.5+40000 for item in hip_data if item[3]=='0']
		hip_mag = [float(item[1]) for item in hip_data if item[3]=='0']
		hip_err = [float(item[2]) for item in hip_data if item[3]=='0']
	
		hip_legend = ['MJD', 'HP', 'err_HP']
		hip_row_table = [hip_MJD, hip_mag, hip_err]
	d_hip = {}
	if n == 0:
		return
	else: 
		d_hip['Hipparcos'] = dict(list(zip(hip_legend,hip_row_table)))
		return d_hip

##############################################################################################################################
def GetData_Monson2011(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/ApJS/193/12/table4')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/ApJS/193/12/table3',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Monson+ 2011 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/ApJS/193/12/table3',Name=star2)[0]
			JD = ma.getdata(result['HJD'])
			Phase = ma.getdata(result['Phase'])
			Jmag = ma.getdata(result['Jmag'])
			err_J = ma.getdata(result['e_Jmag'])
			Hmag = ma.getdata(result['Hmag'])
			err_H = ma.getdata(result['e_Hmag'])
			Kmag = ma.getdata(result['Kmag'])
			err_K = ma.getdata(result['e_Kmag'])
			legend = ['JD','Phase','J','err_J','H','err_H','K','err_K']
			row_table = [JD, Phase, Jmag, err_J, Hmag, err_H, Kmag, err_K]
			d['Monson+ 2011'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Monson2012(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/ApJ/759/146')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/ApJ/759/146/table3',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Monson+ 2012 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/ApJ/759/146/table3',Name=star2)[0]
			MJD = ma.getdata(result['MJD'])
			mag_3_6 = ma.getdata(result['__3.6_'])
			err_3_6 = ma.getdata(result['e__3.6_'])
			mag_4_5 = ma.getdata(result['__4.5_'])
			err_4_5 = ma.getdata(result['e__4.5_'])
			legend = ['MJD', '3.6um', 'err_3.6', '4.5um', 'err_4.5']
			row_table = [MJD, mag_3_6, err_3_6, mag_4_5, err_4_5]
			d['Monson+ 2012'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Nardetto2009(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/A+A/502/951')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/A+A/502/951/table45',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Nardetto+ 2009 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/A+A/502/951/table45',Name=star2)[0]
			RV = ma.getdata(result['RVcc'])
			err_RV = ma.getdata(result['e_RVcc'])
			Phase = ma.getdata(result['Phase'])
			legend = ['Phase','RV','err_RV']
			row_table = [Phase, RV, err_RV]
			d['Nardetto+ 2009'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Petterson2005(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/MNRAS/362/1167')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/MNRAS/362/1167/tablea',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Petterson+ 2005 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/MNRAS/362/1167/tablea',Name=star2)[0]
			JD = ma.getdata(result['JD'])
			RV = ma.getdata(result['RV'])
			Phase = ma.getdata(result['Phase'])
			RV2 = np.array([x for x in RV if not math.isnan(x)])
			if len(RV2) != 0:
				legend = ['JD','Phase', 'RV']
				row_table = [JD,Phase, RV]
				d['Petterson+ 2005'] = dict(list(zip(legend,row_table)))
				return d
			else: 
				JD = ma.getdata(result['JD'])
				RV = ma.getdata(result['FeIIVel'])
				Phase = ma.getdata(result['Phase'])
				legend = ['JD','Phase', 'RV']
				row_table = [JD,Phase, RV]
				d['Petterson+ 2005'] = dict(list(zip(legend,row_table)))
				return d

##############################################################################################################################
def GetData_Storm2004(star):    
	Vizier.ROW_LIMIT = -1
	d = {}
	if star=='SU Cas':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea2',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='EV Sct':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea3',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='Delta Cep':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea4',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='CV Mon':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea5',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='U Sgr':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea6',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='Eta Aql':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea7',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='X Cyg':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea8',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='T Mon':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea9',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='RS Pup':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea10',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	elif star=='SV Vul':
		print('      - downloading Storm+ 2004 data from Vizier...')
		result = Vizier.query_constraints(catalog='J/A+A/415/531/tablea11',Name=star)[0]
		RV = ma.getdata(result['RV'])
		err_RV = ma.getdata(result['e_RV'])
		JD = ma.getdata(result['HJD'])
		Phase = ma.getdata(result['phase'])
		legend = ['JD','Phase','RV','err_RV']
		row_table = [JD, Phase, RV, err_RV]
		d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
		return d
	else:
		return

##############################################################################################################################
def GetData_Storm2011(star):    
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/A+A/534/A94')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		star2=star2[0][0][0]
		result = Vizier.query_constraints(catalog='J/A+A/534/A94/table3',Name=star2)
		if str(result)=='Empty TableList':
			return
		else:
			print('      - downloading Storm+ 2011 data from Vizier...')
			result = Vizier.query_constraints(catalog='J/A+A/534/A94/table3',Name=star2)[0]
			RV = ma.getdata(result['HRV'])
			err_RV = ma.getdata(result['e_HRV'])
			JD = ma.getdata(result['HJD'])
			Phase = ma.getdata(result['Phase'])
			legend = ['JD','Phase','RV','err_RV']
			row_table = [JD, Phase, RV, err_RV]
			d['Storm+ 2011'] = dict(list(zip(legend,row_table)))
			return d

##############################################################################################################################
def GetData_Tycho(star):
	try:
		objet = Simbad.query_objectids(star)
	except:
		objet = Simbad.query_objectids(star)
	obj = list(objet['ID'])
	tyc_name = np.nan
	n = 0
	for elt in obj:
		if 'TYC' in elt:
			tyc_name = elt.split()[1]
	if tyc_name is not np.nan: 
		html_TYCHO = 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Plot/Vgraph/htm?I/239/0&' + tyc_name
		try:
			tyc_html = requests.get(html_TYCHO)
			tyc_soup = BeautifulSoup(tyc_html.text,"html5lib").text.split('\n')
		except:
			tyc_html = requests.get(html_TYCHO)
			tyc_soup = BeautifulSoup(tyc_html.text,"html5lib").text.split('\n')
		tyc_n_init = 0
		for x in tyc_soup:
			tyc_n_init = tyc_soup.index(x)
		tycdata = []
		for i in np.arange(tyc_n_init+1,len(tyc_soup)-2,1):
			tycdata.append(fix_unicode(tyc_soup[i].split('|')))
		a=0
		for i in range(len(tycdata)):
			for j in range(len(tycdata[i])):
				tycdata[i][j]=tycdata[i][j].replace(' ','')
			if tycdata[i]==['z', 'PAz', 'Dmas', 's.e.', 'Fit', 'Flags']:
				a=i
				n = 1
		if n == 1:
			print('      - downloading Tycho data from Vizier...')
			if tycdata[a+1]==['']:
				del tycdata[a+1]
			tyc_data = []
			for i in range(len(tycdata)):
				if i>a:
					tyc_data.append((tycdata[i]))
			tyc_data = [item for item in tyc_data if item[12]=='0']
			tyc_MJD = [float(item[0])-0.5+40000 for item in tyc_data]
			# only Vmag
			tyc_mag_v = [float(item[4]) for item in tyc_data if item[5]!='']
			tyc_err_v = [float(item[5]) for item in tyc_data if item[5]!='']

			# only Bmag:
			tyc_mag_b = [float(item[1]) for item in tyc_data if item[2]!='']
			tyc_err_b = [float(item[2]) for item in tyc_data if item[2]!='']

			tyc_legend = ['MJD', 'V', 'err_V', 'B', 'err_B']
			tyc_row_table = [tyc_MJD, tyc_mag_v, tyc_err_v, tyc_mag_b, tyc_err_b]
	d_tyc = {}
	if n == 0:
		return
	else:
		d_tyc['Tycho'] = dict(list(zip(tyc_legend,tyc_row_table)))
		return d_tyc

##############################################################################################################################
def GetData_WISE(star):
	v = Vizier(columns=['AllWISE', 'W1mag', 'W2mag', 'W3mag', 'W4mag', 'e_W1mag', 'e_W2mag', 'e_W3mag', 'e_W4mag'])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="5s",catalog='II/311/wise')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		result = v.query_region(coordinates(star), radius="5s",catalog='II/311/wise')[0][0]
		print('      - downloading WISE data from Vizier...')
		d = {}
		mag_W1 = np.float32(ma.getdata(result['W1mag']))
		err_W1 = np.float32(ma.getdata(result['e_W1mag']))
		mag_W2 = np.float32(ma.getdata(result['W2mag']))
		err_W2 = np.float32(ma.getdata(result['e_W2mag']))
		mag_W3 = np.float32(ma.getdata(result['W3mag']))
		err_W3 = np.float32(ma.getdata(result['e_W3mag']))
		mag_W4 = np.float32(ma.getdata(result['W4mag']))
		err_W4 = np.float32(ma.getdata(result['e_W4mag']))
		legend = ['MJD', 'W1mag', 'W2mag', 'W3mag', 'W4mag', 'e_W1mag', 'e_W2mag', 'e_W3mag', 'e_W4mag']
		row_table = [None, mag_W1, mag_W2, mag_W3, mag_W4, err_W1, err_W2, err_W3, err_W4]
		d['WISE'] = dict(list(zip(legend,row_table)))
		return d

##############################################################################################################################
def GetData_2MASS(star):
	v = Vizier(columns=['2MASS', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag', 'JD'])
	v.ROW_LIMIT=-1
	result = v.query_region(coordinates(star), radius="2s",catalog='II/246/out')
	d = {}
	if str(result)=='Empty TableList':
		return
	else:
		result = v.query_region(coordinates(star), radius="5s",catalog='II/246/out')[0][0]
		print('      - downloading 2MASS data from Vizier...')
		MJD = result['_tab1_36']-2400000.5
		J = result['Jmag']
		H = result['Hmag']
		Ks = result['Kmag']
		e_J = result['e_Jmag']
		e_H = result['e_Hmag']
		e_Ks = result['e_Kmag']
		legend = ['MJD', 'J', 'H', 'Ks', 'e_J', 'e_H', 'e_Ks']
		row_table = [MJD, J, H, Ks, e_J, e_H, e_Ks]
		d['2MASS'] = dict(list(zip(legend,row_table)))
		return d

##############################################################################################################################
######################################################      AAVSO      #######################################################
##############################################################################################################################
def GetData_AAVSO(star):
	try:
		objet = Simbad.query_objectids(star)
	except:
		objet = Simbad.query_objectids(star)
	obj = list(objet['ID'])
	for elt in obj:
		if 'V*' in elt:
			star2 = elt.replace('V* ', '')
			star2 = star2.replace(' ','+',1)
	dslr_ccd = 'dslr+ccd'
	page = '1'
	filters_available = 'https://www.aavso.org/filters'
	html_aavso = 'https://www.aavso.org/apps/webobs/results/?star=%s&num_results=200&obs_types=%s&page=%s'%(star2,dslr_ccd,page)
	aavso_html = requests.get(html_aavso)
	aavso_soup = BeautifulSoup(aavso_html.text,"html5lib").text.split('\n')
	rani = list(range(len(aavso_soup)))
	rani.reverse()
	pagemax = 1
	for i in rani:
		aavso_soup[i] = aavso_soup[i].replace(' ','')
		if 'NextLast' in str(aavso_soup[i]):
			pagemax = float(aavso_soup[i][aavso_soup[i].index('N')-1])
			if aavso_soup[i][aavso_soup[i].index('N')-2]==',' or aavso_soup[i][aavso_soup[i].index('N')-2]=='.':
				pass
			else:
				pagemax = float('%s'%(aavso_soup[i][aavso_soup[i].index('N')-2])+'%s'%(pagemax))
	page = np.arange(1,pagemax+1,1)
	for i in range(len(page)):
		page[i]=int(page[i])
	aavso_data = []
	for p in page:
		p=int(p)
		html_aavso = 'https://www.aavso.org/apps/webobs/results/?star=%s&num_results=200&obs_types=%s&page=%s'%(star2,dslr_ccd,p)
		aavso_html = requests.get(html_aavso)
		aavso_soup = BeautifulSoup(aavso_html.text,"html5lib").text.split('\n')
		rani = list(range(len(aavso_soup)))
		rani.reverse()
		for i in rani:
			aavso_soup[i] = aavso_soup[i].replace(' ','')
			try:
				if float(aavso_soup[i]):
					aavso_soup[i] = float(aavso_soup[i])
				if aavso_soup[i]<30000. and aavso_soup[i]>30.:
					del aavso_soup[i]
			except:
				if aavso_soup[i]=='B' or aavso_soup[i]=='V' or aavso_soup[i]=='R' or aavso_soup[i]=='I' or aavso_soup[i]=='J' or aavso_soup[i]=='H' or aavso_soup[i]=='K' or aavso_soup[i]=='TG' or aavso_soup[i]=='CV':
					aavso_soup[i]=str(aavso_soup[i])
				else:
					del(aavso_soup[i])
		rani = list(range(len(aavso_soup)))
		rani.reverse()
		for i in rani:
			if isinstance(aavso_soup[i], str):
				if aavso_soup[i-1]>1. and not isinstance(aavso_soup[i-1], str):
					del(aavso_soup[i-1])
		rani = list(range(len(aavso_soup)))
		rani.reverse() 
		for i in rani:
			if isinstance(aavso_soup[i], str) and not isinstance(aavso_soup[i-1],float):
				del(aavso_soup[i])
		for i in range(len(aavso_soup)):
			aavso_data.append(aavso_soup[i])
	# catalog of all dslr+ccd data
	aavso_MJD = []
	aavso_mag = []
	aavso_err = []
	aavso_band = []
	for i in range(len(aavso_data)):
		if aavso_data[i]>2400000. and not isinstance(aavso_data[i], str) and (aavso_data[i]!='0') and (aavso_data[i]!='.00'):
			aavso_MJD.append(aavso_data[i])
		if aavso_data[i]<1. and not isinstance(aavso_data[i], str):
			aavso_err.append(aavso_data[i])
		if aavso_data[i]>1. and aavso_data[i]<30. and not isinstance(aavso_data[i], str):
			aavso_mag.append(aavso_data[i])
		if isinstance(aavso_data[i], str):
			aavso_band.append(aavso_data[i])
	aavso_MJD_b, aavso_MJD_v, aavso_MJD_r, aavso_MJD_i, aavso_MJD_j, aavso_MJD_h, aavso_MJD_k = [],[],[],[],[],[],[]
	aavso_mag_b, aavso_mag_v, aavso_mag_r, aavso_mag_i, aavso_mag_j, aavso_mag_h, aavso_mag_k = [],[],[],[],[],[],[]
	aavso_err_b, aavso_err_v, aavso_err_r, aavso_err_i, aavso_err_j, aavso_err_h, aavso_err_k = [],[],[],[],[],[],[]
	for i in range(len(aavso_band)):
		if aavso_band[i]=='B':
			aavso_MJD_b.append(aavso_MJD[i]-2400000.5)
			aavso_mag_b.append(aavso_mag[i])
			aavso_err_b.append(aavso_err[i])
		if aavso_band[i]=='V':
			aavso_MJD_v.append(aavso_MJD[i]-2400000.5)
			aavso_mag_v.append(aavso_mag[i])
			aavso_err_v.append(aavso_err[i])
		if aavso_band[i]=='R':
			aavso_MJD_r.append(aavso_MJD[i]-2400000.5)
			aavso_mag_r.append(aavso_mag[i])
			aavso_err_r.append(aavso_err[i])
		if aavso_band[i]=='I':
			aavso_MJD_i.append(aavso_MJD[i]-2400000.5)
			aavso_mag_i.append(aavso_mag[i])
			aavso_err_i.append(aavso_err[i])
	aavso_legend = []
	aavso_row_table = []
	if len(aavso_MJD_b)!=0:
		aavso_legend.extend(('MJD_B','B','err_B'))
		aavso_row_table.extend((aavso_MJD_b,aavso_mag_b,aavso_err_b))
	if len(aavso_MJD_v)!=0:
		aavso_legend.extend(('MJD_V','V','err_V'))
		aavso_row_table.extend((aavso_MJD_v,aavso_mag_v,aavso_err_v))
	if len(aavso_MJD_r)!=0:
		aavso_legend.extend(('MJD_R','R','err_R'))
		aavso_row_table.extend((aavso_MJD_r,aavso_mag_r,aavso_err_r))
	if len(aavso_MJD_i)!=0:
		aavso_legend.extend(('MJD_I','I','err_I'))
		aavso_row_table.extend((aavso_MJD_i,aavso_mag_i,aavso_err_i))
	d_aavso = {}
	if len(aavso_legend)==0:
		return
	else:
		print('      - downloading AAVSO data...')
		d_aavso['AAVSO'] = dict(list(zip(aavso_legend,aavso_row_table)))
		return d_aavso

##############################################################################################################################
######################################################      CROCUS      ######################################################
##############################################################################################################################
def Reference_McMaster(star):
	html = requests.get('https://crocus.physics.mcmaster.ca/Cepheid/URL/MW/'+star+'.html',verify=False)
	if str(html)=='<Response [404]>':
		ref_p = []
		ref_v = []
	else:
		soup = BeautifulSoup(html.text,"html5lib")
		ligne = soup.find_all('li')
		list_ref = []
		list_lign = []
		for link in soup.find_all('a'):
			list_ref.append(link.get('href'))
		for ligne in soup.find_all('li'):
			list_lign.append(ligne.get_text())
		list_ref_p = []
		list_ref_v = []
		for i in range(len(list_ref)):
			adress = list_ref[i].split('/')
			fic = adress[4].split('.')
			if ((fic[-1] == '1') or (fic[-1] == '2')) & (fic[-2] == 'p'):
				list_ref_p.append(list_ref[i])
			if ((fic[-1] == '1') or (fic[-1] == '2')) & (fic[-2] == 'v'):
				list_ref_v.append(list_ref[i])
		ref_p = []
		ref_v = []
		for i in np.arange(0,len(list_ref_p),1):
			a = str(list_lign[i].split()[0].split(',')[0])
			b = list_ref_p[i].split('/')[3].split('_')[0][0:4]
			ref_p.append(a + '+ ' + b)
			print(('      - downloading %s data from McMaster...'%ref_p[i]))
		for i in np.arange(0,len(list_ref_v),1):
			a = str(list_lign[i+len(list_ref_p)].split()[0].split(',')[0])
			b = list_ref_v[i].split('/')[3].split('_')[0][0:4]
			ref_v.append(a + '+ ' + b)
			print(('      - downloading %s data from McMaster...'%ref_v[i]))
		ref_p = np.array(ref_p)
		ref_v = np.array(ref_v)
	return ref_p,ref_v

##############################################################################################################################
def GetData_McMaster(star):
	legend_table = [] 
	legend_table.append('Photometry')
	legend_table.append('Radial Velocity')
	try:
		objet = Simbad.query_objectids(star)
	except:
		objet = Simbad.query_objectids(star)
	obj = list(objet['ID'])
	star2 = star.replace(' ','_',1)
	for elt in obj:
		if elt[0:3]=='V* ':
			star2 = elt.replace('V* ', '')
			star2 = star2.replace(' ','_',2)
	if star=='Delta Cep':
		star2='Delta_Cep'
	elif star=='Zeta Gem':
		star2='Zeta_Gem'
	elif star=='Beta Dor':
		star2='Beta_Dor'
	elif star=='Alpha UMi':
		star2='Alpha_UMi'
	elif star=='Eta Aql':
		star2='Eta_Aql'
	ref_p,ref_v = Reference_McMaster(star2)
	crocus = 'https://www.physics.mcmaster.ca/'
	r = requests.get(crocus + 'Cepheid/URL/MW/'+star2+'.html',verify=False)
	if str(r)=='<Response [404]>':
		row_table = ([],[])
		d = dict(list(zip(legend_table,row_table)))
	else:
		soup = BeautifulSoup(r.text,"html5lib")
		list_ref = []
		list_lign = []
		for link in soup.find_all('a'):
			list_ref.append(link.get('href'))
		for ligne in soup.find_all('li'):
			list_lign.append(ligne.get_text())
		data_phot = []
		data_vit = []
		for i in range(len(list_ref)):
			adress = list_ref[i].split('/')
			fic = adress[4].split('.')
			if ((fic[-1] == '1') or (fic[-1] == '2')) & (fic[-2] == 'p'):
				data = unicode_to_liste(requests.get(crocus + list_ref[i],verify=False))
				data = np.array(object_to_tab(data))
				data_i = requests.get(crocus + list_ref[i+1],verify=False)
				unit,code = unit_and_code(data_i)
				data_phot.append([adress[3],data,unit,code])
			if ((fic[-1] == '1') or (fic[-1] == '2')) & (fic[-2] == 'v'):
				data = unicode_to_liste(requests.get(crocus + list_ref[i],verify=False))
				data = np.array(object_to_tab(data))
				data_i = requests.get(crocus + list_ref[i+1],verify=False)
				unit,code = unit_and_code(data_i)
				data_vit.append([adress[3],data,unit,code])        
		for i in range(len(data_phot)):
			data_phot[i][0] = ref_p[i]
		for i in range(len(data_vit)):
			data_vit[i][0] = ref_v[i]
		row_table = (data_phot,data_vit)
		d = dict(list(zip(legend_table,row_table)))
	return d
  
###################################################################################################################################
######################################################      Other data      #######################################################
###################################################################################################################################
def GetData_OtherData(star):
	legend=[]
	row_table=[]
	phot = {}
	rv = {}
	teff = {}
	diam = {}
	if star not in list(OtherData.keys()):
		return
	else:
		print('      - downloading other data from personal database...')
		if OtherData[star]['Teff']!=[]:
			teff['Other Data'] = OtherData[star]['Teff']
		if OtherData[star]['Diameter']!=[]:
			diam['Other Data'] = OtherData[star]['Diameter']
		if OtherData[star]['Radial Velocity']!=[]:
			rv['Other Data'] = OtherData[star]['Radial Velocity']
		if OtherData[star]['Photometry']!=[]:
			phot['Other Data'] = OtherData[star]['Photometry']
	return teff, diam, rv, phot

##############################################################################################################################
######################################################      INFOS      #######################################################
##############################################################################################################################
def GetData_Groenewegen2018(star):    
	v = Vizier(columns=["**"])
	v.ROW_LIMIT = -1
	star2 = v.query_region(coordinates(star), radius="2s",catalog='J/A+A/619/A8/table1')
	d = {}
	if str(star2)=='Empty TableList':
		return
	else:
		result = v.query_region(coordinates(star), radius="2s",catalog='J/A+A/619/A8/table1')[0][0]
		star_type = result['Type']
		period = result['Period']
		Vmag = result['Vmag']
		refVis = result['r_Vmag']
		Jmag = result['Jmag']
		Hmag = result['Hmag']
		Kmag = result['Kmag']
		refIR = result['r_Jmag']
		ebv = result['E_B-V_']
		err_ebv = result['e_E_B-V_']
		if refIR=='1':
			refIR = 'Monson+ 2011'
		elif refIR=='2':
			refIR = 'Laney+ 1992'
		elif refIR=='3':
			refIR = 'Laney+ PC'
		elif refIR=='4':
			refIR = 'Feast+ 2008'
		elif refIR=='5':
			refIR = 'Barnes+ 1997'
		elif refIR=='6':
			refIR = 'Welch+ 1984'
		elif refIR=='7':
			refIR = 'Genovali+ 2014 (templates)'
		elif refIR=='10':
			refIR = '2MASS' 
		else:
			refIR = 'Not enough data (2MASS and others)'
		Jsmag, Hsmag, Ksmag = 0., 0., 0.
		if (refIR=='Laney+ 1992') or (refIR=='Laney+ PC') or (refIR=='Feast+ 2008'):
			err_Jmag = 0.008
			err_Hmag = 0.008
			err_Kmag = 0.008
			# Original in Carter System
			JmagCIT = Jmag - 0.001 - 0.134*(Jmag-Kmag)
			HmagCIT = Hmag + 0.004 - 0.022*(Jmag-Kmag)
			KmagCIT = Kmag - 0.003 - 0.027*(Jmag-Kmag)
			Ksmag = Kmag - 0.019 + 0.001*(Jmag-Kmag)
			Jsmag = Ksmag - 0.020 + 1.068*(Jmag-Kmag)
			Hsmag = Ksmag + 0.034 + 1.000*(Hmag-Kmag)
		elif (refIR=='Monson+ 2011'):
			err_Jmag = 0.008
			err_Hmag = 0.008
			err_Kmag = 0.008
			Ksmag = Kmag - 0.019 + 0.001*(Jmag-Kmag)
			Jsmag = Ksmag - 0.020 + 1.068*(Jmag-Kmag)
			Hsmag = Ksmag + 0.034 + 1.000*(Hmag-Kmag)
		elif (refIR=='Welch+ 1984') or (refIR=='Barnes+ 1997'):
			err_Jmag = 0.01
			err_Hmag = 0.01
			err_Kmag = 0.01
			Ksmag = Kmag - 0.019 + 0.001*(Jmag-Kmag)
			Jsmag = Ksmag - 0.020 + 1.068*(Jmag-Kmag)
			Hsmag = Ksmag + 0.034 + 1.000*(Hmag-Kmag)
		elif refIR=='Genovali+ 2014':
			err_Jmag = 0.025
			err_Hmag = 0.025
			err_Kmag = 0.025
			Ksmag = Kmag
			Jsmag = Jsmag
			Hsmag = Hsmag
		elif refIR=='2MASS': 
			err_Jmag = 0.25
			err_Hmag = 0.25
			err_Kmag = 0.25
			Ksmag = Kmag
			Jsmag = Jsmag
			Hsmag = Hsmag
		else:
			err_Jmag = 0.025
			err_Hmag = 0.025
			err_Kmag = 0.025
			Ksmag = Kmag
			Jsmag = Jsmag
			Hsmag = Hsmag
		legend = ['Name', 'Type', 'Period', 'Vmag', 'Jmag_original', 'Hmag_original', 'Kmag_original', 'Jmag_2MASS', 'Hmag_2MASS', 'Kmag_2MASS', 'E(B-V)', 'refVis', 'refIR']
		row_table = [star, star_type, period, [Vmag, 0.02], [Jmag, err_Jmag], [Hmag, err_Hmag], [Kmag, err_Kmag], [Jsmag, err_Jmag], [Hsmag, err_Hmag], [Ksmag, err_Kmag], [ebv, err_ebv], refVis, refIR]
		d = dict(list(zip(legend,row_table)))
		return d

##############################################################################################################################
def info_cepheids(star):
	print('      - Getting informations about the star...')
	# Calculate distance from Hipparcos parallax
	# search with simbad icrs coordinates in a catalog: result = Vizier.query_region(coord.SkyCoord("22h29m10.2650s +58d24m54.713s",ICRS),radius=Angle(1, "arcsec"),catalog="I/239/hip_main")
	try:
		coords = coordinates(star)
	except:
		coords = coordinates(star)
	RA_star = coords.ra.value
	DEC_star = coords.dec.value
	# Find Period and epoch of maximum/minimum light in VSX
	v = Vizier(columns=['OID', 'Name', 'Epoch', 'Period'])
	try:
		result = v.query_region(coords, radius="0.1s", catalog='B/vsx/vsx')[0]
		mjd0_star = result[0][2]-2400000.5 # epoch of maximum or minimum in days (HJD-2400000.5)
		period_star = result[0][3] # period in days
	except:
		try:
			result = v.query_region(coords, radius="1s", catalog='B/vsx/vsx')[0]
			mjd0_star = result[0][2]-2400000.5 # epoch of maximum or minimum in days (HJD-2400000.5)
			period_star = result[0][3] # period in days
		except:
			try:
				result = v.query_region(coords, radius="3s", catalog='B/vsx/vsx')[0]
				mjd0_star = result[0][2]-2400000.5 # epoch of maximum or minimum in days (HJD-2400000.5)
				period_star = result[0][3] # period in days
			except:
				result = v.query_region(coords, radius="50s", catalog='B/vsx/vsx')[0]
				mjd0_star = result[0][2]-2400000.5 # epoch of maximum or minimum in days (HJD-2400000.5)
				period_star = result[0][3] # period in days
	eGDR3 = GetData_eGDR3(star)
	try:
		Plx_Gaia = np.float32(eGDR3['Plx eGDR3'][0])
		e_Plx_Gaia = np.float32(eGDR3['Plx eGDR3'][1])
	except:
		Plx_Gaia = eGDR3['Plx eGDR3'][0]
		e_Plx_Gaia = eGDR3['Plx eGDR3'][1]
	# Find Color excess (if available) from Kovtyukh 2008
	if GetData_Groenewegen2018(star) is not None:
		excess_star = GetData_Groenewegen2018(star)['E(B-V)'][0]
	else:
		v = Vizier(columns=['Name', 'E(B-V)'])
		Vizier.ROW_LIMIT = -1
		result = v.query_region(coords, radius="2s", catalog='J/MNRAS/389/1336/table2')
		if str(result)=='Empty TableList':
			excess_star = 0.2
		else:
			excess_star = result[0][0][1]
	# Find temperature (if available) from Luck 2011
	v = Vizier(columns=["Name"])
	Vizier.ROW_LIMIT = -1
	star2 = v.query_region(coords, radius="2s",catalog='J/AJ/142/136/table4')
	if str(star2)=='Empty TableList':
		v = Vizier(columns=["Star","Teff"])
		result = v.query_region(coords, radius="2s",catalog='J/MNRAS/428/3252')
		if str(result)=='Empty TableList':
			teff_star = 5200.
		else:
			teff_star = result[0][0][1]
	else:
		star2=star2[0][0][0]
		v = Vizier(columns=["Name",'Teff'])
		result = v.query_constraints(catalog='J/AJ/142/136/table3',Name=star2)
		teff_star = result[0][0][1]
	dic_info_star = {'MJD0': mjd0_star, 'Period': period_star, 'Plx eGDR3': [Plx_Gaia, e_Plx_Gaia], 'Coordinates': [RA_star, DEC_star], 'E(B-V)': excess_star, 'Teff': teff_star}
	return dic_info_star

##############################################################################################################################