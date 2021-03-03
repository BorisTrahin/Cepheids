import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy as unp
import pickle
import scipy.optimize
import time

### Loads the stars from Trahin database:
print('\n Loading Trahin database: _data_spips.dpy...')
f1 = open('_data_spips.dpy', 'rb')   
database = pickle.load(f1, encoding='latin1')
f1.close()
stars = sorted(list(database.keys())) 


###########################################################################
### Returns the phase of an observation:
def phase(mjd, params):
# mjd = epoch of the observation
# params = dictionnary with at least ephemeris MJD0 and period PERIOD. Can include period changes such as:
# period = PERIOD + (mjd-MJD0) * PERIOD1/(24*3600*365.25) + sum( (mjd-MJD0/1e4)**n*PERIODn )
	mjd_mjd0 = mjd - params['MJD0']
	period = params['PERIOD']*np.ones(len(mjd))
	if 'PERIOD1' in params.keys():
		period += (mjd_mjd0)*params['PERIOD1']/(24*3600*365.25)
	px = []
	for x in params.keys():
		if x.startswith('PERIOD') and x!='PERIOD' and x!='PERIOD1' and not x.startswith('PERIOD '):
			px.append(x)
	if len(px)>0:
		for k in px:
			x = float(k.split('PERIOD')[1])
			period += ((mjd_mjd0)/1e4)**x*params[k]
	return ((mjd_mjd0)/period)%1.0

###########################################################################
### Returns the mean magnitudes from the fit of the light-curves and derive mean magnitude from flux-average:
def mean_mag(star, band, fmode = 4, show_LC = False, save_LC = False):
	### Updates main parameters for specific stars:
	params={'PERIOD': database[star]['Info']['Period'], 'MJD0': database[star]['Info']['MJD0']}
	# Period changes for some long period stars
	if star=='S Vul':
		params={'MJD0': 48333.87377856825,'PERIOD': 68.5243,'PERIOD1': -2814,'PERIOD2': 0.064,'PERIOD3': 1.96,'PERIOD4': -1.83,'PERIOD5': 0.66,'PERIOD6': 0.225}
	elif star=='SV Vul':
		params={'MJD0': 48307.75795377858,'PERIOD': 44.9386, 'PERIOD1': -769,'PERIOD2': 0.432, 'PERIOD3': 0.507, 'PERIOD4': -1.233, 'PERIOD5': 0.388, 'PERIOD6': 0.1360}
	elif star == 'AQ Pup':
		params = {'MJD0': 54587.06559742254, 'PERIOD': 30.16601, 'PERIOD1': 139.1, 'PERIOD2': 0.00347 }
	elif star == 'RS Pup':
		params = {'MJD0': 54215.99341561817, 'PERIOD': 41.45710, 'PERIOD1': -227.3, 'PERIOD2': 0.3865, 'PERIOD3': 1.3013, 'PERIOD4': 1.2644, 'PERIOD5': 0.4013 }
	elif star == 'GY Sge':
		params = {'MJD0': 48311.29679835036, 'PERIOD': 51.5059, 'PERIOD1': 728, 'PERIOD2': 4.48, 'PERIOD3': -2.95, 'PERIOD4': -23.54, 'PERIOD5': 31.96, 'PERIOD6': -9.90}
	elif star == 'T Mon':
		params = {'MJD0': 43783.93472354179,'PERIOD': 27.02530,'PERIOD1': 17.6}
	elif star == 'U Car':
		params = {'MJD0': 48337.59629466054,'PERIOD': 38.8679,'PERIOD1': 504, 'PERIOD2': -0.517, 'PERIOD3': -1.195,'PERIOD4': 2.07, 'PERIOD5': 2.92, 'PERIOD6': -3.55, 'PERIOD7': -2.19, 'PERIOD8': 2.27} 
	elif star == 'l Car':
		params = {'MJD0': 47773.9285738106,'PERIOD': 35.56987, 'PERIOD1': -44.2,'PERIOD2': -0.0372,'PERIOD3': 0.0464}
	### Loads J, H, K data:
	if band in ['J', 'H', 'K']:
		band_name = band + '_CTIO_ANDICAM'
		### Takes all CTIO data points in the requested band with value > 0.5 (removes obvious outliers):
		mjddata = [x[0]      for x in database[star]['Data'] if x[2]==band_name and x[3]>0.5]
		ydata =   [x[3]      for x in database[star]['Data'] if x[2]==band_name and x[3]>0.5]
		e_ydata = [x[4]      for x in database[star]['Data'] if x[2]==band_name and x[3]>0.5]
		refs =    [x[1][5::] for x in database[star]['Data'] if x[2]==band_name and x[3]>0.5]
	### Loads V and I data:
	elif band in ['V', 'I']:
		if band == 'V':
			band_name = band + '_GCPD_Johnson'
		elif band == 'I':
			band_name = band + '_GCPD_Cousins'
		### Takes all GCPD_Johnson data points from Berdnikov in the requested band:
		mjddata = [x[0]      for x in database[star]['Data'] if x[2]==band_name]
		ydata   = [x[3]      for x in database[star]['Data'] if x[2]==band_name]
		e_ydata = [x[4]*3    for x in database[star]['Data'] if x[2]==band_name]
		refs    = [x[1][5::] for x in database[star]['Data'] if x[2]==band_name]

	list_ref = []
	for r in refs:
		if r in list_ref:
			pass
		else:
			list_ref.append(r)

	nb_pt = len(mjddata)
	if nb_pt <= 8:
		mean_mag_from_flux, e_mean_mag_from_flux, list_ref, fmode = 0.0, 0.0, '---', 0.
	else:
		### Phase the data:
		xdata = phase(np.array(mjddata), params)
		xdata = list(xdata)
		A0 = (max(ydata)+min(ydata))/2
		A1 = abs(A0-min(ydata))
		### Sorts all the data:
		Z = list(zip(xdata, ydata))
		Z2 = list(zip(xdata, e_ydata))
		Z.sort()
		Z2.sort()
		### Magnitudes converted into flux:
		### Flux for Vega (zero-point) in erg/cm2/s/A.
		F_vega = {'V': 363.1e-11, 'I': 112.6e-11, 'J': 30.524e-11, 'H': 12.003e-11, 'K': 4.4793e-11}
		flux = [F_vega[band]*10**(-ufloat(ydata[i], e_ydata[i])/2.5) for i in range(len(ydata))]
		ydata_flux   = [flux[i].nominal_value for i in range(len(flux))]
		e_ydata_flux = [flux[i].std_dev       for i in range(len(flux))]
		A0_flux = (max(ydata_flux)+min(ydata_flux))/2
		A1_flux = abs(A0_flux-min(ydata_flux))
		Z_flux = list(zip(xdata, ydata_flux))
		Z2_flux = list(zip(xdata, e_ydata_flux))
		Z_flux.sort()
		Z2_flux.sort()

		xdata        = np.array([X[0] for X in Z])
		ydata        = np.array([X[1] for X in Z])
		e_ydata      = np.array([X[1] for X in Z2])
		ydata_flux   = np.array([X[1] for X in Z_flux])
		e_ydata_flux = np.array([X[1] for X in Z2_flux])

		fmode = fmode

		### Fourier modes:
		fg      = {'A0':A0,      'A1':A1,      'A2':0., 'PHI1':0., 'PHI2':0., 'WAV':1.}
		fg_flux = {'A0':A0_flux, 'A1':A1_flux, 'A2':0., 'PHI1':0., 'PHI2':0., 'WAV':1.}
		for p in range(2, fmode+1):
			fg.update({'A%i'%p: 0., 'PHI%i'%p: 0.})
			fg_flux.update({'A%i'%p: 0., 'PHI%i'%p: 0.})

		### Fits the light curve in magnitudes:
		fit = leastsqFit(fourier, xdata, fg, ydata, err=e_ydata, verbose=2, doNotFit=['WAV'])
		fit = randomParam(fit)
		mean_mag = np.mean(fit['model'])
		e_mean_mag = 2*np.sqrt(np.sum([(ydata[i] - fit['model'][int(xdata[i]*1000)])**2 for i in range(len(xdata))]))/(len(xdata)-(2*fmode+2))
		ampl = abs( max(fit['model']) - min(fit['model']))

		### Fits the light curve in flux:
		fit_flux = leastsqFit(fourier, xdata, fg_flux, ydata_flux, err=e_ydata_flux, verbose=2, doNotFit=['WAV'])
		fit_flux = randomParam(fit_flux)
		mean_flux = np.mean(fit_flux['model'])
		e_mean_flux = 2*np.sqrt(np.sum([(ydata_flux[i] - fit_flux['model'][int(xdata[i]*1000)])**2 for i in range(len(xdata))]))/(len(xdata)-(2*fmode+2))
		mean_M = -2.5*unp.log10(ufloat(mean_flux, e_mean_flux)/F_vega[band])
		mean_mag_from_flux, e_mean_mag_from_flux = mean_M.nominal_value, mean_M.std_dev
		if e_mean_mag_from_flux < 0.006:
			e_mean_mag_from_flux = 0.006

		plt.figure(figsize=(14, 5))
		plt.subplots_adjust(left=0.05, right=0.98, top=0.94, bottom=0.11)

		a1 = plt.subplot(121)
		a1.plot(np.arange(0, 1, 0.001), fit_flux['model'], '-', linewidth=1., color='k', label='Best fit for flux')
		a1.errorbar(xdata, ydata_flux, yerr=e_ydata_flux, color='k', fmt='o', markersize=3, elinewidth=0.5)
		a1.plot([0., 1.], [mean_flux, mean_flux], '--', color='limegreen', linewidth=1., label='Mean flux')
		a1.fill_between([0., 1.], [mean_flux + e_mean_flux, mean_flux + e_mean_flux], [mean_flux - e_mean_flux, mean_flux - e_mean_flux], facecolor='limegreen', alpha=0.25)
		a1.set_xlim(0., 1.)
		a1.set_ylabel('Flux in %s (erg/cm2/s/A) '%band)
		a1.set_xlabel('phase')	
		a1.set_title(' %s  (P=%.3f d)   (Fourier Mode: %i)' %(star, params['PERIOD'], fmode))
		plt.gca().invert_yaxis()
		plt.legend(fontsize=9, loc='best')

		a2 = plt.subplot(122)
		a2.plot(np.arange(0, 1, 0.001), fit['model'], '-',  color='b', linewidth=1., label='Best fit (value p.t.p.: %.2f mag)'%ampl)
		a2.errorbar(xdata, ydata, yerr=e_ydata, color='k', fmt='o', markersize=3, elinewidth=0.5, label=[x for x in list_ref])
		a2.plot([0., 1.], [mean_mag_from_flux, mean_mag_from_flux], '--', color='limegreen', linewidth=1., label='Mean mag from flux: $m_K$ = %.3f +/- %.3f mag' %(mean_mag_from_flux, e_mean_mag_from_flux))
		a2.fill_between([0., 1.], [mean_mag_from_flux + e_mean_mag_from_flux, mean_mag_from_flux + e_mean_mag_from_flux], [mean_mag_from_flux - e_mean_mag_from_flux, mean_mag_from_flux - e_mean_mag_from_flux], facecolor='limegreen', alpha=0.25)
		a2.set_xlim(0., 1.)
		a2.set_ylabel('%s (mag)'%band)
		a2.set_xlabel('phase')	
		a2.set_title(' %s  (P=%.3f d)' %(star, params['PERIOD']))
		plt.gca().invert_yaxis()
		plt.legend(fontsize=9, loc='best')
		if save_LC == True:
			plt.savefig('%s_%s.pdf'%(band, star))
		if show_LC == True:
			plt.show()
		plt.close()
			
	return (mean_mag_from_flux, e_mean_mag_from_flux, nb_pt, list_ref, fmode)

###########################################################################
### Function used to fit the light-curves:
def leastsqFit(func, x, params, y, err=None, fitOnly=None, verbose=False, doNotFit=[], epsfcn=1e-7, ftol=1e-5, fullOutput=True, normalizedUncer=True, follow=None, maxfev=200, showBest=True):
	global Ncalls
	if fitOnly is None:
		if len(doNotFit)>0:
			fitOnly = [x for x in list(params.keys()) if x not in doNotFit]
		else:
			fitOnly = list(params.keys())
		fitOnly.sort() # makes some display nicer
	pfit = [params[k] for k in fitOnly]
	pfix = {}
	for k in list(params.keys()):
		if k not in fitOnly:
			pfix[k]=params[k]
	Ncalls=0
	plsq, cov, info, mesg, ier = \
			  scipy.optimize.leastsq(_fitFunc, pfit, args=(fitOnly,x,y,err,func,pfix,verbose,follow,), full_output=True, epsfcn=epsfcn, ftol=ftol, maxfev=maxfev)
	if cov is None:
		cov = np.zeros((len(fitOnly), len(fitOnly)))
	for i,k in enumerate(fitOnly):
		pfix[k] = plsq[i]
	xmodel = np.arange(0,1,0.001)
	model = func(xmodel,pfix)
	tmp = _fitFunc(plsq, fitOnly, x, y, err, func, pfix)
	try:
		chi2 = (np.array(tmp)**2).sum()
	except:
		chi2=0.0
		for x in tmp:
			chi2+=np.sum(x**2)
	reducedChi2 = chi2/float(np.sum([1 if np.isscalar(i) else
									 len(i) for i in tmp])-len(pfit)+1)
	if not np.isscalar(reducedChi2):
		reducedChi2 = np.mean(reducedChi2)
	uncer = {}
	for k in list(pfix.keys()):
		if not k in fitOnly:
			uncer[k]=0 # not fitted, uncertatinties to 0
		else:
			i = fitOnly.index(k)
			if cov is None:
				uncer[k]= -1
			else:
				uncer[k]= np.sqrt(np.abs(np.diag(cov)[i]))
				if normalizedUncer:
					uncer[k] *= np.sqrt(reducedChi2)
	if fullOutput:
		if normalizedUncer:
			try:
				cov *= reducedChi2
			except:
				pass
		cor = np.sqrt(np.diag(cov))
		cor = cor[:,None]*cor[None,:]
		cor = cov/cor
		pfix={'best':pfix, 'uncer':uncer, 'chi2':reducedChi2, 'model':model, 'cov':cov, 'fitOnly':fitOnly, 'info':info, 'cor':cor, 'x':x, 'y':y, 'err':err, 'func':func}
	return pfix

###########################################################################
### Function used to fit the light-curves:
def randomParam(fit, N=None, x=None):
	if N is None:
		N = len(fit['x'])
	m = np.array([fit['best'][k] for k in fit['fitOnly']])
	res = [] # list of dictionnaries
	for k in range(N):
		p = dict(list(zip(fit['fitOnly'],np.random.multivariate_normal(m, fit['cov']))))
		p.update({k:fit['best'][k] for k in list(fit['best'].keys()) if not k in
				 fit['fitOnly']})
		res.append(p)
	ymin, ymax = None, None
	tmp = []
	if x is None:
		x = fit['x']
	for r in res:
		tmp.append(fit['func'](x, r))
	tmp = np.array(tmp)
	fit['r_param'] = res
	fit['r_ym1s'] = np.percentile(tmp, 16, axis=0)
	fit['r_yp1s'] = np.percentile(tmp, 84, axis=0)
	fit['r_x'] = x
	fit['r_y'] = fit['func'](x, fit['best'])
	fit['all_y'] = tmp
	return fit
randomParam = randomParam

###########################################################################
verboseTime=time.time()
Ncalls=0
### Function used to fit the light-curves:
def _fitFunc(pfit, pfitKeys, x, y, err=None, func=None, pfix=None, verbose=False, follow=None):
	global verboseTime, Ncalls
	Ncalls+=1
	params = {}
	for i,k in enumerate(pfitKeys):
		params[k]=pfit[i]
	for k in pfix:
		params[k]=pfix[k]
	if err is None:
		err = np.ones(np.array(y).shape)
	if type(y)==np.ndarray and type(err)==np.ndarray:
		if len(err.shape)==2:
			tmp = func(x,params)
			res = np.dot(np.dot(tmp-y, err), tmp-y)
			res = np.ones(len(y))*np.sqrt(res/len(y))
		else:
			y = np.array(y)
			res= ((func(x,params)-y)/err).flatten()
	else:
		res = []
		tmp = func(x,params)
		if np.isscalar(err):
			err = 0*y + err
		for k in range(len(y)):
			df = (np.array(tmp[k])-np.array(y[k]))/np.array(err[k])
			try:
				res.extend(list(df))
			except:
				res.append(df)
	if verbose and time.time()>(verboseTime+5):
		verboseTime = time.time()
		try:
			chi2=(res**2).sum/(len(res)-len(pfit)+1.0)
		except:
			chi2 = 0
			N = 0
			res2 = []
			for r in res:
				if np.isscalar(r):
					chi2 += r**2
					N+=1
					res2.append(r)
				else:
					chi2 += np.sum(np.array(r)**2)
					N+=len(r)
					res2.extend(list(r))
			res = res2
	return res

###########################################################################
### Function used to fit the light-curves:
def fourier(x, params):
	phi = [k for k in list(params.keys()) if k[:3]=='PHI']
	res = np.copy(x)
	res *=0
	for f in phi:
		res += params['A'+f[3:]]*np.cos(float(f[3:])*2*np.pi*x/params['WAV']+ params[f])
	if 'A0' in params:
		res += params['A0']
	return res


