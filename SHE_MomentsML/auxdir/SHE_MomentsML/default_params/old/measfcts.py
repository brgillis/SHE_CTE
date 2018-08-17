import megalut.meas


def default(catalog, stampsize, gain):
	"""
	Default measfct, runs on "img".
	"""	
	
	# HSM adamom
	catalog = megalut.meas.galsim_adamom.measfct(catalog, stampsize=stampsize, variant="wider")
	
	# And skystats
	catalog = megalut.meas.skystats.measfct(catalog, stampsize=stampsize)
	
	# And snr
	catalog = megalut.meas.snr_2hlr.measfct(catalog, gain=gain)
	
	
	return catalog


default_groupcols = [
'adamom_flag',
'adamom_flux',
'adamom_x',
'adamom_y',
'adamom_g1',
'adamom_g2',
'adamom_sigma',
'adamom_rho4',
'adamom_size',
'skystd',
'skymad',
'skymean',
'skymed',
'skystampsum',
'skyflag',
'snr',
]

default_removecols = []

