# Database of Galactic Cepheids

Create a homogeneous database of radial velocities, photometry, effective temperature and angular diameters of Galactic Cepheids.

Last update: December 27th, 2020

Author: Boris Trahin

Last use: Breuval et al. 2021

The .dpy files contains data for more than 300 Cepheids.

Names are VSX names such as Delta_Cep, L_Car, RS_Pup, etc.

## cepheids_data.py:
Collect data from online database such as Vizier, McMaster, AAVSO, etc..., using SIMBAD names (then converted into coordinates):
* coordinates: Get the coordinates of a star in SkyCoord form.
* GetData_Anderson2014 : Get Anderson+ 2014 (J/A+A/566/L10, RV) data from Vizier for the specified star.
* GetData_Anderson2016 : Get Anderson+ 2016 (J/ApJS/226/18, RV) data from Vizier for the specified star.
* GetData_Barnes2005 : Get Barnes+ 2005 (J/ApJS/156/227, RV) data from Vizier for the specified star.
* GetData_Berdnikov2008 : Get Berdnikov+ 2008 (II/285, V,B-V,V-R,R-I,V-I) data from Vizier for the specified star.
* GetData_Berdnikov2014 : Get Berdnikov+ 2014 (J/PAZh/40/147, V,B,I) data from Vizier for the specified star.
* GetData_Borgniet2019: Get Borgniet+ 2019 (J/A+A/631/A37, RV) data from Vizier for the specified star.
* GetData_Feast2008 : Get Feast+ 2008 (J/MNRAS/386/2115, J,H,K,L,V-K) data from Vizier for the specified star.
* GetData_GaiaDR2 : Get Gaia DR2 (ESA Archive, G, Gbp, Grp) data for the specified star.
* GetData_GaiaeDR3 : Get Gaia eDR3 (I/350/gaiaedr3) parameters from Vizier for the specified star.
* GetData_Gorynya1992 : Get Gorynya+ 1992-1998 (III/229, RV) data from Vizier for the specified star.
* GetData_Groenewegen2013 : Get Groenewegen+ 2013 (J/A+A/550/A70, RV) data from Vizier for the specified star.
* GetData_Hipparcos : Get Hipparcos (ESA Archive, HP) data for the specified star.
* GetData_Monson2011 : Get Monson+ 2011 (J/ApJS/193/12/table4, J,H,K) data from Vizier for the specified star.
* GetData_Monson2012 : Get Monson+ 2012 (J/ApJ/759/146, Spitzer I1 and I2) data from Vizier for the specified star.
* GetData_Nardetto2009 : Get Nardetto+ 2009 (J/A+A/502/951, RV) data from Vizier for the specified star (HARPS data of 8 galactic cepheids with cross-correlation).
* GetData_Petterson2005 : Get Petterson+ 2005 (J/MNRAS/362/1167, RV) data from Vizier for the specified star.
* GetData_Storm2004 : Get Storm+ 2004 (J/A+A/415/531/, RV) data from Vizier for the specified star.
* GetData_Storm2011 : Get Storm+ 2011 (J/A+A/534/A94, RV) data from Vizier for the specified star.
* GetData_Tycho : Get Tycho (ESA Archive, B,V) data for the specified star.
* GetData_WISE : Get WISE (II/311/wise, W1, W2, W3, W4) data from Vizier for the specified star.
* GetData_2MASS : Get 'II/246/out (II/311/wise, J, H, Ks) data from Vizier for the specified star.
* GetData_AAVSO : Get AAVSO WebObs data for the specified star.
* Reference_McMaster : Create 2 arrays containing references of photometry and radial velocity available in McMaster for the specified star.
* GetData_McMaster : Get McMaster data available in the references obtained with Reference_McMaster for the selected star.
* GetData_OtherData : Other data from own observations or not in the electronic form.
* GetData_Groenewegen2018 : Get Groenewegen+ 2018 (J/A+A/619/A8/) mean magnitudes in the 2MASS system from Vizier for the specified star.
* info_cepheids : Get some informations of the selected star (period, distance, ...).


## OtherData.py:
Data not present in an online database. Collected directly from articles or personal observations.


## data_original.py:
Create a dictionary \_data_original.dpy with the original data for each Cepheids, as they were published and collected using cepheids_data.py and OtherData.py.
* readDb_obs : read data from .dpy file for a specified star.
* writeDb_obs : write data in .dpy file for a specified star.
* Dictionnaire_data : Create a dictionnary with data (from McMaster and Vizier) and information about the star (useful in SPIPS).
For each star: 2 dictionaries: Data and Info.
* Info: Informations such as Period (in days), E(B-V), mean Teff, MJD0, eGDR3 distance, etc...
* Data: dictionnaires of data for each reference.


## data_spips.py:
Create a dictionary \_data_spips.dpy with the data for each Cepheids (using \_data_original.dpy file created by data_original.py), in a homogeneous structure, readable by SPIPS (Mérand et al. 2015, Breitfelder et al. 2016, Kervella et al. 2018, Trahin 2019).
For each star: 2 dictionaries: 'Data' and 'Info'.
* Info: Informations such as Period (in days), E(B-V), mean Teff, MJD0, eGDR3 distance, etc...
* Data: a list containing the data in the form: \[MJD, 'type; reference', 'SVOfilter or base', value, error]
* MJD = epoch of the observation in days
* type = vrad, mag, teff, color or diam. Reference is the reference of the data ex.
	* ex: 'mag; Berdnikov+ 2008' for magnitudes published in Berdnikov et al. 2008.
* SVO filter = band and system used for the photometric observations, as it appears in the [SVO Filter Profile Service](http://svo2.cab.inta-csic.es/theory/fps/).
	* ex: 'V_GCPD_Johnson' for V magnitude in GCPD_Johnson system, Ks_2MASS for Ks magnitude in 2MASS system, etc.
* base = Wavelength and base of the interferometric measurements of angular diameter.
	* ex: \[1.675, 130] for PIONIER measurements (H band) using a 130m baseline.
* value = value of the observation.
* error = error of the observation.
