# ----------------------------------------------------------------------------------------------------------------------------
# python utility functions for QENS fits
# both backscattering and TOF 
# amongst other models of 2 coupled diffusive states, global fit of several q-values
# fit of D2O solvent background or use as fixed parameters in sample fits
# the energy axis is always in micro eV (converted from meV at load level for TOF data)
# by Tilo Seydel, project started in summer 2017
# requires python 3.5.2 or higher
# ----------------------------------------------------------------------------------------------------------------------------
import numpy as np
import h5py

# Define units used
Dunits   =   10 / 6.58212

# ----------------------------------------------------------------------------------------------------------------------------
def gaussian( x, sig, x0 ):
	return ( 1. / ( np.sqrt( 2. * np.pi ) * sig ) ) * np.exp( - np.power( x - x0, 2. ) / ( 2. * np.power( sig, 2. ) ) )

# ----------------------------------------------------------------------------------------------------------------------------
def lorentzian( x, gamma, x0 ):
	return gamma / np.pi / ( np.power( x - x0, 2. ) + np.power( gamma, 2. ) )

# ---------------------------------------------------------------------------------
def V( x, alpha, gamma ): # Voigt line shape, gamma is the Lorentzian HWHM
	from scipy.special import wofz # Faddeeva function exp(-z**2)*erfc(-i*z)
	import pylab
	sigma = alpha # / np.sqrt( 2. * np.log(2.) ) the input is the Gaussian sigma, not HWHM, from the Gaussian fit 
	return np.real( wofz( (x + 1j*gamma) / sigma / np.sqrt(2.) ) ) / sigma / np.sqrt( 2.*np.pi )

# --------------------------------------------------------------------------------
def group_4channels( data, *args ): # group by 4 channels in energy
	class out:
		pass
	try:
		out.q    = data.q
		out.hw   = (   data.hw[ 0::4 ]       +   data.hw[ 1::4 ]       +   data.hw[ 2::4 ]       +   data.hw[ 3::4 ]       )      / float( 4 )
		out.sqw  = (  data.sqw[ 0::4, : ]    +  data.sqw[ 1::4, : ]    +  data.sqw[ 2::4, : ]    +  data.sqw[ 3::4, : ]    )      / float( 4 )
		out.dsqw = ( data.dsqw[ 0::4, : ]**2 + data.dsqw[ 1::4, : ]**2 + data.dsqw[ 2::4, : ]**2 + data.dsqw[ 3::4, : ]**2 )**0.5 / float( 4 )
	except:
		out.y    = (    data[ 0::4 ]    +  data[ 1::4 ]   +   data[ 2::4 ]       +   data[ 3::4 ]   )      / float( 4 )
		try:
			for count, thing in enumerate( args ):
				out.dy   = (   thing[ 0::4 ]**2 + thing[ 1::4 ]**2 + thing[ 2::4 ]**2 + thing[ 3::4 ]**2 )**0.5 / float( 4 )
		except:
			pass
	return out

# --------------------------------------------------------------------------------
def group_2channels( data, *args ): # group by 2 channels in energy
	class out:
		pass
	try:
		out.q    = data.q
		out.hw   = (   data.hw[ 0::2 ]       +   data.hw[ 1::2 ]       )      / float( 2 )
		out.sqw  = (  data.sqw[ 0::2, : ]    +  data.sqw[ 1::2, : ]    )      / float( 2 )
		out.dsqw = ( data.dsqw[ 0::2, : ]**2 + data.dsqw[ 1::2, : ]**2 )**0.5 / float( 2 )
	except:
		out.y    = (    data[ 0::2 ]    +  data[ 1::2 ]    )      / float( 2 )
		try:
			for count, thing in enumerate( args ):
				out.dy   = (   thing[ 0::2 ]**2 + thing[ 1::2 ]**2 )**0.5 / float( 2 )
		except:
			pass
	return out

# ---------------------------------------------------------------------------------
# model for 2 coupled diffusive modes (switching diffusion) according to Roosen-Runge, Bicout, and Barrat:
# cf. Grimaldo et al., PCCP 2015
def two_state( D1, D2, t1in, t2in, q, qn ):
	class out:
		pass
	hbar = 658.212    # reduced Planck constant in micro eV ps
	t1 = t1in #/ hbar  # the linewidths in the fits are in units of hbar*omega, therefore, the residence times have to be taken in the same units
	t2 = t2in #/ hbar
	G1 = D1 * q[ qn ]**2 # assume 2 coupled diffusion processes
	G2 = D2 * q[ qn ]**2
	L  = np.sqrt( np.power( G1 - G2 + 1. / t1 - 1. / t2, 2 ) + 4. / ( t1 * t2 ) ) 
	g1 = ( G1 + 1. / t1 + G2 + 1. / t2 + L ) / 2.
	g2 = ( G1 + 1. / t1 + G2 + 1. / t2 - L ) / 2.
	a  = 1. / ( g2 - g1 ) * (  t1 / ( t1 + t2 ) * ( G2 + 1. / t1 + 1. / t2 - g1 ) \
			 	 + t2 / ( t1 + t2 ) * ( G1 + 1. / t1 + 1. / t2 - g1 ) ) 
	out.g1 = g1
	out.g2 = g2
	out.a  = a
	return( out )

# ---------------------------------------------------------------------------------
def array2vector( x ):
	try:
		len(x[0])
	except:
		out = x
	else:
		out = x[ :, 0 ]
		for k in range( 1, len( x[1] ) ):
			out = np.concatenate( ( out, x[ :, k ] ) )
	return out

# ---------------------------------------------------------------------------------
def vector2array( x, q ):
	try:
		len(q)
	except:
		out = x
	else:
		m   = int( np.round( len(x)/len(q) ) )
		n   = len( q )
		out = np.zeros(( m, n ))
		for k in range( 0, n ):
			out[:,k] = x[ (0 + m*k) : (m + m*k) ] 
	return out

# --------------------------------------------------------------------------------------------------------
# ------------- WRAPPER FUNCTION FOR SCIPY.CURVE_FIT ----------------------------------------------------
# ---- build matrices from individual fit parameters to conveniently interface to the model functions ---- 
# --------------------------------------------------------------------------------------------------------
#
# "wrap" a variable number of fit parameters into the fit function - wrapper for S(q,omega) parser
def wrapper_fit_func_parser( hw, q, qn, r, N, *args, **kwargs ):
	a  = list( args[0][:N] )
	fp = np.asarray( a )  # fit parameters
	return model_sqw_interface_parser( hw, fp, q, qn, r, **kwargs )

# ---------------------------------------------------------------------------------------------------
# "wrap" a variable number of fit parameters into the resolution function - the fit here is per one q
def wrapper_resolution_func( hw, N, *args ):
	a   = list( args[0][:N] )
	fp  = np.asarray( a )
	b   = fp[ 0 ]    # global constant background
	a0  = fp[ 1::3 ] # every 3rd entry: amplitude
	sig = fp[ 2::3 ] # every 3rd entry: width
	x0  = fp[ 3::3 ] # every 3rd entry: x-offset
	return resolution_model_IN5( hw, b, a0, sig, x0 )

# ---------------------------------------------------------------------------------------------------
# "wrap" a variable number of fit parameters into the resolution function - the fit here is per one q
def wrapper_1Gauss_fix_resolution_func( hw, b0, sig1, x01, N, *args ):
	a   = list( args[0][:N] )
	fp  = np.asarray( a )
	b   = b0         # fixed global constant background
	a0  = fp[ 0 ]    # first Gauss free amplitude
	sig = sig1	 # first Gauss fixed width
	x0  = x01        # first Gauss fixed center
	a0  = np.hstack( [ a0,  fp[ 1::3 ] ] ) # free amplitude for other Gauss
	sig = np.hstack( [ sig, fp[ 2::3 ] ] ) # free width for other Gauss
	x0  = np.hstack( [ x0,  fp[ 3::3 ] ] ) # free center for other Gauss
	return resolution_model_IN5( hw, b, a0, sig, x0 )

# ---------------------------------------------------------------------------------
# ------------------------ MODEL FUNCTIONS FOR FITTING ----------------------------
# ---------------------------------------------------------------------------------
#
# requirement for "curve_fit": use hw and sqw as long vectors with all q "in a line"
def model_sqw_interface_parser( hw, fp, q, qn, r, **kwargs ):
	hwm  = vector2array( hw, q )
	out  = model_sqw_parser( hwm, fp, q, qn, r, **kwargs ) # pass the keyword arguments unchanged
	outm = array2vector( np.transpose( out.y ) )
	return( outm )

# ---------------------------------------------------------------------------------
def model_sqw_subplots_parser( hw, fp, q, qn, r, **kwargs ):
	hwm  = vector2array( hw, q )
	out  = model_sqw_parser( hwm, fp, q, qn, r, **kwargs ) # pass the keyword arguments unchanged
	return( out )

# ---------------------------------------------------------------------------------
#
# build the model S(q,omega) by parsing keywords:
def model_sqw_parser( hw, fp, q, qn, r, **kwargs ): # for a flexible configuration of the model according to keyword arguments
#
	class out:      # object of returned values
		pass
	ID2O = np.zeros( np.shape( hw ) )
	GD2O = np.zeros( np.shape( hw ) )
	for key, value in kwargs.items():
		if key == 'SolventIntensities':
			ID2O = value
		if key == 'SolventWidths':
			GD2O = value
		if key == 'SolventDirac':
			DD2O = value
		if key == 'BGSlope':
			bgs = value
		if key == 'BGFlat':
			bgf = value
		if key == 'intD1':
			D1i = value
			intflag = 1
		else: # If no intD1 argument is passed
			intflag = 0
		if key == 'intTau1':
			tau1i = value
	for key, value in kwargs.items():
		if key == 'SQWModel':
			if value == 'IN5SolventPerQ':
				k      = qn[ 0 ]
				p      = np.zeros( 7 )
				p[ 0 ] = fp[ 0 ]    # intensity of first solvent Lorentzian
				p[ 1 ] = fp[ 1 ]    # width of first solvent Lorentzian
				p[ 2 ] = fp[ 2 ]    # intensity of second solvent Lorentzian
				p[ 3 ] = fp[ 3 ]    # width of second solvent Lorentzian
				p[ 4 ] = fp[ 4 ]    # amplitude of guest molecule Dirac (elastic line e.g. from empty can)
				p[ 5 ] = fp[ 5 ]    # background slope
				p[ 6 ] = fp[ 6 ]    # background offset
			if value == 'BATSSolventFixedGammaPerQ':
				k      = qn[ 0 ]
				p      = np.zeros( 5 )
				p[ 0 ] = fp[ 0 ]    # intensity of first solvent Lorentzian
				p[ 1 ] = GD2O       # width of first solvent Lorentzian
				p[ 2 ] = fp[ 1 ]    # amplitude of guest molecule Dirac (elastic line e.g. from empty can)
				p[ 3 ] = fp[ 2 ]    # background slope
				p[ 4 ] = fp[ 3 ]    # background offset
			if value == 'BATSSolventPerQ':
				k      = qn[ 0 ]
				p      = np.zeros( 5 )
				p[ 0 ] = fp[ 0 ]    # intensity of first solvent Lorentzian
				p[ 1 ] = fp[ 1 ]    # width of first solvent Lorentzian
				p[ 2 ] = fp[ 2 ]    # amplitude of guest molecule Dirac (elastic line e.g. from empty can)
				p[ 3 ] = fp[ 3 ]    # background slope
				p[ 4 ] = fp[ 4 ]    # background offset

			elif value == 'BATSGlobal2State':
				ng = 5                                              # number of global parameters: D, D1, tau1, D2, tau2 
				p  = np.zeros( [ len(qn), 11 ] )                    # parameter matrix to be returned	
				D   = fp[0]
				D1  = fp[1]
				D2  = fp[2]
				t1  = fp[3]
				t2  = fp[4]
				tsp = two_state( D1, D2, t1, t2, q, qn )
				for k in range( 0, len( qn ) ):
					p[ k, 0 ] = fp[ ng + k ]       * fp[ ng + k + len(qn) ]                     # intensity center-of-mass Lorentzian: beta(q) * A_0(q)
					p[ k, 1 ] = D * q[ k ]**2                                                   # width center-of-mass Lorentzian
					p[ k, 2 ] = fp[ ng + k ] * ( 1 - fp[ ng + k + len(qn) ] ) * tsp.a[k]        # intensities 1st coupled Lorentzian: beta(q)*(1-A_0(q))*tsp.a(q)
					p[ k, 3 ] = D * q[ k ]**2 + tsp.g1[ k ]		   		            # 1st coupled Lorentzian width convoluted with center-of-mass
					p[ k, 4 ] = fp[ ng + k ] * ( 1 - fp[ ng + k + len(qn) ] ) * ( 1- tsp.a[k] ) # intensities 2nd coupled Lorentzian
					p[ k, 5 ] = D * q[ k ]**2 + tsp.g2[ k ]				            # 2nd coupled Lorentzian width convoluted with center-of-mass
					p[ k, 6 ] = ID2O[ k ]                                                       # intensity solvent Lorentzian
					p[ k, 7 ] = GD2O[ k ]                                                       # width solvent Lorentzian
					p[ k, 8 ] = DD2O[ k ]                                                       # amplitude of empty can Dirac
					p[ k, 9 ] = bgs[ k ]                                                        # background slope
					p[ k,10 ] = bgf[ k ]                                                        # background offset

			elif value == 'BATSGlobalFickianInternalPerQ':
				ng = 1                                              # number of global parameters: D, D1, tau1, D2, tau2 
				p  = np.zeros( [ len(qn), 9 ] )                    # parameter matrix to be returned	
				D   = fp[0]
				for k in range( 0, len( qn ) ):
					p[ k, 0 ] = fp[ ng + k ]       * fp[ ng + k + len(qn) ]                     # intensity center-of-mass Lorentzian: beta(q) * A_0(q)
					p[ k, 1 ] = D * q[ k ]**2                                                   # width center-of-mass Lorentzian
					p[ k, 2 ] = fp[ ng + k ] * ( 1 - fp[ ng + k + len(qn) ] )                   # intensities internal Lorentzian: beta(q)*(1-A_0(q))
					p[ k, 3 ] = D * q[ k ]**2 +      fp[ ng + k + len(qn) * 2 ]		    # internal Lorentzian width convoluted with center-of-mass
					p[ k, 4 ] = ID2O[ k ]                                                       # intensity solvent Lorentzian
					p[ k, 5 ] = GD2O[ k ]                                                       # width solvent Lorentzian
					p[ k, 6 ] = DD2O[ k ]                                                       # amplitude of empty can Dirac
					p[ k, 7 ] = bgs[ k ]                                                        # background slope
					p[ k, 8 ] = bgf[ k ]                                                        # background offset

			elif value == 'IN16BSolventPerQ':
				p      = np.zeros( 2 )
				p[ 0 ] = fp[ 0 ]    # intensity of first solvent Lorentzian
				p[ 1 ] = GD2O       # width of first solvent Lorentzian
			elif value == 'IN16BSolventPerQCanSubtracted':
				k      = qn[ 0 ]
				p      = np.zeros( 3 )
				p[ 0 ] = fp[ 0 ]    # intensity of first solvent Lorentzian
				p[ 1 ] = GD2O       # width of first solvent Lorentzian fixed from
			elif value == 'IN16BSolventAndDiracPerQBackground':
				p      = np.zeros( 5 )
				p[ 0 ] = fp[ 0 ]    # intensity of first solvent Lorentzian
				p[ 1 ] = GD2O       # width of first solvent Lorentzian
				p[ 2 ] = fp[ 1 ]    # Dirac intensity
				p[ 3 ] = fp[ 2 ]    # background slope
				p[ 4 ] = fp[ 3 ]    # background slope
			elif value == 'IN16BSolventAndDiracPerQ':
				p      = np.zeros( 3 )
				p[ 0 ] = fp[ 0 ]    # intensity of first solvent Lorentzian
				p[ 1 ] = GD2O       # width of first solvent Lorentzian
				p[ 2 ] = fp[ 1 ]    # Dirac intensity
			elif value == 'IN16BGlobalFickianInternalFree':
				p = np.zeros( [ len(qn), 6 ] )                                 # parameter matrix to be returned	
				for k in range( 0, len( qn ) ):
					p[ k, 0 ]  = fp[ 1+k ] * fp[ 1+k+len(qn) ]             # beta(q) * A_0(q)
					p[ k, 1 ]  = fp[ 0 ] * q[ k ]**2                       # Fickian center-of-mass diffusion
					p[ k, 2 ]  = fp[ 1+k ] * ( 1 - fp[ 1+k+len(qn) ] )     # beta(q) * ( 1 - A_0(q) )
					p[ k, 3 ]  = fp[ 0 ] * q[ k ]**2 + fp[ 1+k+2*len(qn) ] # internal diffusion convoluted with center-of-mass diffusion
					p[ k, 4 ]  = ID2O[ k ]                                 # intensity solvent Lorentzian
					p[ k, 5 ]  = GD2O[ k ]                                 # width solvent Lorentzian
			elif value == 'IN16BGlobalFickianInternalFreeFixedDirac':
				p = np.zeros( [ len(qn), 7 ] )                                 # parameter matrix to be returned	
				for k in range( 0, len( qn ) ):
					p[ k, 0 ]  = fp[ 1+k ] * fp[ 1+k+len(qn) ]             # beta(q) * A_0(q)
					p[ k, 1 ]  = fp[ 0 ] * q[ k ]**2                       # Fickian center-of-mass diffusion
					p[ k, 2 ]  = fp[ 1+k ] * ( 1 - fp[ 1+k+len(qn) ] )     # beta(q) * ( 1 - A_0(q) )
					p[ k, 3 ]  = fp[ 0 ] * q[ k ]**2 + fp[ 1+k+2*len(qn) ] # internal diffusion convoluted with center-of-mass diffusion
					p[ k, 4 ]  = ID2O[ k ]                                 # intensity solvent Lorentzian
					p[ k, 5 ]  = GD2O[ k ]                                 # width solvent Lorentzian
					p[ k, 6 ]  = DD2O[ k ]                                 # Dirac intensity
			elif value == 'IN16BGlobalFickianInternalFreeDiracFree':
				p = np.zeros( [ len(qn), 9 ] )                                 # parameter matrix to be returned	
				for k in range( 0, len( qn ) ):
					p[ k, 0 ]  = fp[ 1+k ] * fp[ 1+k+len(qn) ]             # beta(q) * A_0(q)
					p[ k, 1 ]  = fp[ 0 ] * q[ k ]**2                       # Fickian center-of-mass diffusion
					p[ k, 2 ]  = fp[ 1+k ] * ( 1 - fp[ 1+k+len(qn) ] )     # beta(q) * ( 1 - A_0(q) )
					p[ k, 3 ]  = fp[ 0 ] * q[ k ]**2 + fp[ 1+k+2*len(qn) ] # internal diffusion convoluted with center-of-mass diffusion
					p[ k, 4 ]  = ID2O[ k ]                                 # intensity solvent Lorentzian
					p[ k, 5 ]  = GD2O[ k ]                                 # width solvent Lorentzian
					p[ k, 6 ]  = fp[ 1+k+3*len(qn) ]                       # Dirac intensity as free fit parameter
					p[ k, 7 ]  = bgs[ k ]                                  # background slope
					p[ k, 8 ]  = bgf[ k ]                                  # background offset
			elif value == 'IN16BGlobalFickianInternalJumpBackground':
				nf = 3                                                           # number of fixed parameters
				p  = np.zeros( [ len(qn), 8 ] )                                  # parameter matrix to be returned	
				for k in range( 0, len( qn ) ):
					gamma = fp[ 0 ] * q[ k ]**2
					Gamma = fp[ 1 ] * q[ k ]**2 / ( 1  + fp[ 1 ] * q[ k ]**2 * fp[ 2 ] )
					p[ k, 0 ]  = fp[ nf+k ] * fp[ nf+k+len(qn) ]             # beta(q) * A_0(q)
					p[ k, 1 ]  = gamma                                       # Fickian center-of-mass diffusion
					p[ k, 2 ]  = fp[ nf+k ] * ( 1 - fp[ nf+k+len(qn) ] )     # beta(q) * ( 1 - A_0(q) )
					p[ k, 3 ]  = gamma + Gamma                               # internal diffusion convoluted with center-of-mass diffusion
					p[ k, 4 ]  = ID2O[ k ]                                   # intensity solvent Lorentzian
					p[ k, 5 ]  = GD2O[ k ]                                   # width solvent Lorentzian
					p[ k, 6 ]  = bgs[ k ]                                    # background slope
					p[ k, 7 ]  = bgf[ k ]                                    # background offset
			elif value == 'IN16BGlobalFickianInternalJump':
				if intflag:
					nf = 1
				else:
					nf = 3                                                           # number of fixed parameters	
				p  = np.zeros( [ len(qn), 6 ] )                                  # parameter matrix to be returned	
				for k in range( 0, len( qn ) ):
					gamma = fp[ 0 ] * q[ k ]**2
					if intflag:
						Gamma = D1i / Dunits * q[ k ]**2 / ( 1  + D1i / Dunits * q[ k ]**2 * tau1i )
					else:
						Gamma = fp[ 1 ] * q[ k ]**2 / ( 1  + fp[ 1 ] * q[ k ]**2 * fp[ 2 ] )
					p[ k, 0 ]  = fp[ nf+k ] * fp[ nf+k+len(qn) ]             # beta(q) * A_0(q)
					p[ k, 1 ]  = gamma                                       # Fickian center-of-mass diffusion
					p[ k, 2 ]  = fp[ nf+k ] * ( 1 - fp[ nf+k+len(qn) ] )     # beta(q) * ( 1 - A_0(q) )
					p[ k, 3 ]  = gamma + Gamma                               # internal diffusion convoluted with center-of-mass diffusion
					p[ k, 4 ]  = ID2O[ k ]                                   # intensity solvent Lorentzian
					p[ k, 5 ]  = GD2O[ k ]                                   # width solvent Lorentzian
			elif value == 'IN16BGlobalFickianInternalJumpDiracFree':
				nf = 3                                                           # number of fixed parameters
				p  = np.zeros( [ len(qn), 9 ] )                                  # parameter matrix to be returned	
				for k in range( 0, len( qn ) ):
					p[ k, 0 ]  = fp[ nf+k ] * fp[ nf+k+len(qn) ]             # beta(q) * A_0(q)
					p[ k, 1 ]  = fp[ 0 ] * q[ k ]**2                         # Fickian center-of-mass diffusion
					p[ k, 2 ]  = fp[ nf+k ] * ( 1 - fp[ nf+k+len(qn) ] )     # beta(q) * ( 1 - A_0(q) )
					p[ k, 3 ]  = ( fp[ 0 ] + fp[ 1 ] ) * q[ k ]**2 / ( 1  + ( fp[ 0 ] + fp[ 1 ] ) * q[ k ]**2 * fp[ 2 ] ) # internal diffusion convoluted with center-of-mass diffusion
					p[ k, 4 ]  = ID2O[ k ]                           # intensity solvent Lorentzian
					p[ k, 5 ]  = GD2O[ k ]                                   # width solvent Lorentzian
					p[ k, 6 ]  = fp[ nf+k+2*len(qn) ]                        # Dirac intensity as free fit parameter
					p[ k, 7 ]  = bgs[ k ]                                    # background slope
					p[ k, 8 ]  = bgf[ k ]                                    # background offset
			elif value == 'IN16BProteinPerQ':
				p = np.zeros( 6 )                      # parameter matrix to be returned
				p[ 0 ]  = fp[ 0 ] * fp[ 1 ]            # beta(q) * A_0(q)
				p[ 1 ]  = fp[ 2 ]                      # center-of-mass diffusion at this q
				p[ 2 ]  = fp[ 0 ] * ( 1 - fp[ 1 ] )    # beta(q) * ( 1 - A_0(q) )
				p[ 3 ]  = fp[ 2 ]  + fp[ 3 ]           # internal diffusion width convoluted with center-of-mass diffusion at this q
				p[ 4 ]  = ID2O                         # intensity solvent Lorentzian at this q
				p[ 5 ]  = GD2O                         # width solvent Lorentzian atthis q

			elif value == 'IN5GuestPerQFixedSolvent':
				k = qn[0]                                          # fit parameters: Lorentzian amplitudes, width of guest molecule Lorentzian, amplitude of Dirac
				p = np.zeros( 7 )                                  # parameter matrix to be returned	
				p[ 0 ] = fp[ 0 ]                                   # intensity of the guest molecule Lorentzian
				p[ 1 ] = fp[ 2 ]                                   # width guest molecule Lorentzian
				p[ 2 ] = fp[ 1 ] * ID2O[ 0 ]                    # intensity 1st solvent Lorentzian
				p[ 3 ] = GD2O[ 0 ]                              # width 1st solvent Lorentzian
				p[ 4 ] = fp[ 1 ] * ID2O[ 1 ]                    # intensity 2nd solvent Lorentzian
				p[ 5 ] = GD2O[ 1 ]                              # width 2nd solvent Lorentzian
				p[ 6 ] = fp[ 3 ]                                   # amplitude of guest molecule Dirac
			elif value == 'BATSGuest1LorGlobalSolventAmpFixed':
				k = qn[0]                                          # fit parameters: Lorentzian amplitudes, width of guest molecule Lorentzian, amplitude of Dirac
				p = np.zeros( [ len(qn), 7 ] )                     # parameter matrix to be returned	
				for k in range( 0, len( qn ) ):
					p[ k, 0 ]  = fp[ 2+k ]                    # intensity of the guest molecule Lorentzian
					p[ k, 1 ]  = fp[ 0 ] * q[ k ]**2 / ( 1 + fp[ 0 ] * q[ k ]**2 * fp[ 1 ] )                 # width guest molecule Lorentzian
					p[ k, 2 ]  = ID2O[ k ]                            # intensity 1st solvent Lorentzian
					p[ k, 3 ]  = GD2O[ k ]                            # width 1st solvent Lorentzian
					p[ k, 4 ]  = fp[ 2+k+len(qn) ]                    # amplitude of guest molecule Dirac
					p[ k, 5 ]  = bgs[ k ]                             # background slope
					p[ k, 6 ]  = bgf[ k ]                             # background offset
			elif value == 'IN5Guest2LorGlobalSolventAmpFixed':
				k = qn[0]                                          # fit parameters: Lorentzian amplitudes, width of guest molecule Lorentzian, amplitude of Dirac
				p = np.zeros( [ len(qn), 11 ] )  # parameter matrix to be returned	
				for k in range( 0, len( qn ) ):
					p[ k, 0 ]  = fp[ 0 ] * fp[ 9+k ]                   # intensity of the guest molecule Lorentzian
					p[ k, 1 ]  = fp[ 1 ] * q[ k ]**2 / ( 1 + fp[ 1 ] * q[ k ]**2 * fp[ 2 ] )                       # width guest molecule Lorentzian
					p[ k, 2 ]  = fp[ 3 ] * fp[ 9+k ]                  # intensity of the guest molecule Lorentzian
					p[ k, 3 ]  = fp[ 4 ] * q[ k ]**2 / ( 1 + fp[ 4 ] * q[ k ]**2 * fp[ 5 ] )  + fp[ 6 ]                     # width guest molecule Lorentzian
					p[ k, 4 ]  = ID2O[ k, 0 ]                    # intensity 1st solvent Lorentzian
					p[ k, 5 ]  = GD2O[ k, 0 ]                              # width 1st solvent Lorentzian
					p[ k, 6 ]  = ID2O[ k, 1 ]                    # intensity 2nd solvent Lorentzian
					p[ k, 7 ]  = GD2O[ k, 1 ]                              # width 2nd solvent Lorentzian
					p[ k, 8 ]  = fp[ 7 ] * np.exp( - q[ k ]**2 * fp[ 8 ] )  * fp[ 9+k ]                      # amplitude of guest molecule Dirac
					p[ k, 9 ]  = bgs[ k ]                                  # background slope
					p[ k, 10 ] = bgf[ k ]                             # background offset
			elif value == 'IN5GuestPerQFixedSolventBG':
				k = qn[0]                                          # fit parameters: Lorentzian amplitudes, width of guest molecule Lorentzian, amplitude of Dirac
				p = np.zeros( 9 )                                  # parameter matrix to be returned	
				p[ 0 ] = fp[ 0 ]                                   # intensity of the guest molecule Lorentzian
				p[ 1 ] = fp[ 2 ]                                   # width guest molecule Lorentzian
				p[ 2 ] = fp[ 1 ] * ID2O[ 0 ]                    # intensity 1st solvent Lorentzian
				p[ 3 ] = GD2O[ 0 ]                              # width 1st solvent Lorentzian
				p[ 4 ] = fp[ 1 ] * ID2O[ 1 ]                    # intensity 2nd solvent Lorentzian
				p[ 5 ] = GD2O[ 1 ]                              # width 2nd solvent Lorentzian
				p[ 6 ] = fp[ 3 ]                              # Dirac intentsity
				p[ 7 ] = bgs                                  # background slope
				p[ 8 ] = bgf                                  # background offset
			elif value == 'IN5GuestPerQFixedPhiFixedSolventBG':
				k = qn[0]                                          # fit parameters: Lorentzian amplitudes, width of guest molecule Lorentzian, amplitude of Dirac
				p = np.zeros( 9 )                                  # parameter matrix to be returned	
				p[ 0 ] = fp[ 0 ]                                   # intensity of the guest molecule Lorentzian
				p[ 1 ] = fp[ 1 ]                                   # width guest molecule Lorentzian
				p[ 2 ] = ID2O[ 0 ]                    # intensity 1st solvent Lorentzian
				p[ 3 ] = GD2O[ 0 ]                              # width 1st solvent Lorentzian
				p[ 4 ] = ID2O[ 1 ]                    # intensity 2nd solvent Lorentzian
				p[ 5 ] = GD2O[ 1 ]                              # width 2nd solvent Lorentzian
				p[ 6 ] = fp[ 2 ]                              # Dirac intentsity
				p[ 7 ] = bgs                                  # background slope
				p[ 8 ] = bgf                                  # background offset
			elif value == 'IN5Guest2LorPerQFixedPhiFixedSolventBG':
				k = qn[0]                                          # fit parameters: Lorentzian amplitudes, width of guest molecule Lorentzian, amplitude of Dirac
				p = np.zeros( 11 )                                  # parameter matrix to be returned	
				p[ 0 ]  = fp[ 0 ]                                   # intensity of the guest molecule Lorentzian
				p[ 1 ]  = fp[ 1 ]                                   # width guest molecule Lorentzian
				p[ 2 ]  = fp[ 2 ]                                   # intensity of the guest molecule Lorentzian
				p[ 3 ]  = fp[ 3 ]                                   # width guest molecule Lorentzian
				p[ 4 ]  = ID2O[ 0 ]                    # intensity 1st solvent Lorentzian
				p[ 5 ]  = GD2O[ 0 ]                              # width 1st solvent Lorentzian
				p[ 6 ]  = ID2O[ 1 ]                    # intensity 2nd solvent Lorentzian
				p[ 7 ]  = GD2O[ 1 ]                              # width 2nd solvent Lorentzian
				p[ 8 ]  = fp[ 4 ]                              # Dirac intentsity
				p[ 9 ]  = bgs                                  # background slope
				p[ 10 ] = bgf                                  # background offset
			elif value == 'BATSGuest1LorPerQFixedPhiFixedSolventBG':
				k = qn[0]                                          # fit parameters: Lorentzian amplitudes, width of guest molecule Lorentzian, amplitude of Dirac
				p = np.zeros( 7 )                                  # parameter matrix to be returned	
				p[ 0 ]  = fp[ 0 ]                                   # intensity of the guest molecule Lorentzian
				p[ 1 ]  = fp[ 1 ]                                   # width guest molecule Lorentzian
				p[ 2 ]  = ID2O                                 # intensity 1st solvent Lorentzian
				p[ 3 ]  = GD2O                                 # width 1st solvent Lorentzian
				p[ 4 ]  = fp[ 2 ]                              # Dirac intentsity
				p[ 5 ]  = bgs                                  # background slope
				p[ 6 ]  = bgf                                  # background offset
	return( convoluted_model_parser( hw, q, qn, p, r, **kwargs ) )

# ---------------------------------------------------------------------------------
# ----------- CONVOLUTED MODELS WITH RESOLUTION FUNCTION -------------------------
# ---------------------------------------------------------------------------------
def sumNVoigts( x, l, r, p ):
	# build sum of Voigt profiles according to number of Gaussians in resolution model
	# take resolution matrix r at q-index l, use Lorenztian width p
	# attention if used with IN16B resolution: reshaping of r is necessary in calling macro/function
	y = np.zeros( np.size( x ) )
	for k in range( 0, len( r[0, 1::3] ) ): # loop over Gaussians in resolution function
		y = y + r[ l, 3*k+1 ] * V( x - r[ l, 3*k+3 ], r[ l, 3*k+2 ], p )  # r[:,0] is the constant background b0 not emplpoyed in convolution
	return( y )

# ---------------------------------------------------------------------------------
def convoluted_model_empty_can( hw, p, qn, r ): 
	# input energy axis, used q-index, matrices of model and resolution parameters
	# asuming resolution function consisting of several Gaussians
	y = p[0]   * sumNVoigts( hw,    qn,    r, 0 ) + p[1] 
	return( y )

# ------------------------------------------------------------------------------------------------------------------------------------
def convoluted_model_parser( hwq, q, qn, pq, r, **kwargs ): 
	# input matrix of energy axes, vector of q-values, vector of used q-indices, matrices of model and resolution parameters
	# list of keywords
	# assuming resolution function consisting of a sum of Gaussians
	# The model is always a sum of 
	# a given number >=0 of Lorentzians plus optionally an elastic contribution ("Dirac") plus optionally sloped and flat backgrounds
	#
	class out:
		pass
	#
	# initialize internal variables:
	y   = None # function values to be returned
	nL  = 0    # number of Lorentzians
	uD  = 0    # use of an elastic contribution ("Dirac")
	sB  = 0    # use of a sloped background
	fB  = 0    # use of a constant background
	#
	# loop over the keywords relevant for this subroutine:
	for key, value in kwargs.items():
		if key == 'NumberLorentzians': # total number of Lorentzians including e.g. solvent contributions
			if value:
				nL = value
		if key == 'UseDirac':
			if value:
				uD = 1
		if key == 'UseSlopedBackground':
			if value:
				sB = 1
		if key == 'UseFlatBackground':
			if value:
				fB = 1
	if ( np.size(qn) == 1 ):
		out.L =[ [ [ 0 for x in range(nL) ] for y in range( len(qn) ) ] for z in range( len(hwq) ) ] # [] inititialize list
	else:
		out.L =[ [ [ 0 for x in range(nL) ] for y in range( len(qn) ) ] for z in range( len(hwq[:,0]) ) ] # [] inititialize list
	#
	# loop over all q:
	for k in range( 0, len( qn ) ):
		idx = 0    # fit parameter index
		if ( np.size(qn) == 1 ): # in case of fit per q
			p  = pq
			hw = hwq
		else: # one-dinensional slice of p and q matrices
			p  =  pq[ k, : ]
			hw = hwq[ :, qn[k] ]
		y = np.zeros( np.shape( hw ) )
		# loop over Lorentzians:
		for l in range( 0, nL ):
			if ( np.size( qn ) == 1 ):
				yy  = p[ idx ] * sumNVoigts( hw, qn,    r, p[ idx + 1 ] )  # add a Lorentzian
			else:
				yy  = p[ idx ] * sumNVoigts( hw, qn[k], r, p[ idx + 1 ] )  # add a Lorentzian
			y   = y + yy # add result to total fit
			idx = idx + 2 # increment parameter index
			# store component for later plotting:
			out.L[l][k][:] = np.transpose( yy ) #.append( np.transpose( yy ) )
		if uD:
			if ( np.size( qn ) == 1 ):
				yy  = p[ idx ] * sumNVoigts( hw, qn,    r, 0 )            # add elastic Line (Dirac)
			else:
				yy  = p[ idx ] * sumNVoigts( hw, qn[k], r, 0 )            # add elastic Line (Dirac)
			y   = y + yy # add result to total fit
			idx = idx + 1
			# store component for later plotting:
			if k == 0:
				out.D = np.transpose( yy )
			else:
				out.D = np.vstack( [ out.D, np.transpose( yy ) ] )
		if sB:
			yy  = p[ idx ] * hw                                     # add sloped background
			y   = y + yy # add result to total fit
			idx = idx + 1
			if k == 0:
				out.S = np.transpose( yy )
			else:
				out.S = np.vstack( [ out.S, np.transpose( yy ) ] )
		if fB:
			yy  = np.zeros( np.shape( y ) ) + p[ idx ]              # add flat background
			y   = y + yy # add result to total fit
			idx = idx + 1
			if k == 0:
				out.F = np.transpose( yy )
			else:
				out.F = np.vstack( [ out.F, np.transpose( yy ) ] )
		if k == 0:
			out.y = np.transpose( y )
		else:
			out.y = np.vstack( [ out.y, np.transpose( y ) ] )
	return( out )

# ---------------------------------------------------------------------------------
# ----------- LOAD ROUTINES -------------------------------------------------------
# ---------------------------------------------------------------------------------
def loadINX( filename ):   # read the old ILL "inx" format
	# initialize some variables:
	f     = open( filename, 'r' )
	k     = 1	# line counter for entire block
	l     = 0       # line counter for (hw, sqw, dsqw)-data lines
	q     = None
	hw    = None
	sqw   = None
	dsqw  = None
	hwm   = None
	sqwm  = None
	dsqwm = None
	class out:      # object of returned values
		pass
	for line in f:   # loop over all lines in the ASCII file
		columns = line.split()  # this is the text line from the file
		if k == 1:         # read number of data lines to folow in block
			n = float( columns[ -1 ] )
		if k == 3:          # read q-value
			if q is None:
				q = float( columns[ 0 ] )
			else:
				q = np.append( q,  float( columns[ 0 ] ) )
		if (k > 4) & (k <= n + 4):  # between block headers
			l += 1              # increase (hw, sqw, dsqw)-data line counter
			if l == 1:          # first line of new q-block
				hw   = float( columns[ 0 ] )
				sqw  = float( columns[ 1 ] )
				dsqw = float( columns[ 2 ] )
			else:               # all other lines of q-block
				hw   = np.vstack( [hw,   float( columns[ 0 ] ) ] )
				sqw  = np.vstack( [sqw,  float( columns[ 1 ] ) ] )
				dsqw = np.vstack( [dsqw, float( columns[ 2 ] ) ] )
		k += 1
		if k == (n + 5):    # this is where a new q-block starts
			k = 1       # reset block line counter
			l = 0       # reset data line counter
			if hwm is None:
				hwm   = hw
				sqwm  = sqw
				dsqwm = dsqw
			else:
				hwm   = np.hstack( [   hwm,   hw ] )
				sqwm  = np.hstack( [  sqwm,  sqw ] )
				dsqwm = np.hstack( [ dsqwm, dsqw ] )
	indices  = np.argsort( q )    # indices for sorting by ascending q
	out.q    = np.array(     q[    indices ] )
	out.hw   = np.array(   hwm[ :, indices ] ) * 1e3 # in units of micro-eV
	out.sqw  = np.array(  sqwm[ :, indices ] )
	out.dsqw = np.array( dsqwm[ :, indices ] )
	return( out )

# ---------------------------------------------------------------------------------
def align_SD( ws, vanad ):    # align single detectors (for BATS data) based on Vanadium
	indSD  = np.argmax( vanad.sqw[:,0] ) # peak position in single detector channels
	indPSD = np.argmax( vanad.sqw[:,5] ) # peak position in PSD channels
	print( 'SD peak channel: ', indSD )
	print( 'PSD peak channel: ', indPSD )
	sqwsd  = np.zeros( np.shape(  ws.sqw[:,0:2] ) )
	dsqwsd = np.zeros( np.shape( ws.dsqw[:,0:2] ) )
	sqwsd[  (indPSD-indSD):, : ]  = ws.sqw[  :-(indPSD-indSD) , 0:2 ]
	dsqwsd[ (indPSD-indSD):, : ]  = ws.dsqw[ :-(indPSD-indSD) , 0:2 ]
	ws.sqw[  :,0:2] = sqwsd 
	ws.dsqw[ :,0:2] = dsqwsd 
	return( ws )

# ---------------------------------------------------------------------------------
def align_all( ws, **kwargs ):    # align all detectors (for BATS data) based on Vanadium
	import copy
	vana_ws = copy.deepcopy( ws )
	print( 'Aligning peaks in workspace.' )
	for key, value in kwargs.items():
		if key == 'Vana':
			if value: # align based on vanadium
				vana_ws = copy.deepcopy( value )
				print( '... based on provided alignment workspace' )
	sqwsd   = np.zeros( np.shape( ws.sqw ) )
	dsqwsd  = np.zeros( np.shape( ws.dsqw ) )
	for k in range( 0, len( ws.q ) ):                         # loop over all spectra in workspace
		indpeak  = np.argmax( vana_ws.sqw[:,k] )          # peak position in single detector channels
		indzero  = np.argmin( np.absolute( vana_ws.hw ) ) # index of zero energy
		print( 'Peak channel: ', indpeak )
		print( 'Zero energy channel: ', indzero )
		sqwsd[  (indpeak-indzero):, k ] = ws.sqw[  :-(indpeak-indzero), k ]
		dsqwsd[ (indpeak-indzero):, k ] = ws.dsqw[ :-(indpeak-indzero), k ]
		ws.sqw  = sqwsd 
		ws.dsqw = dsqwsd 
	return( ws )

# ---------------------------------------------------------------------------------
def loadIN16B( filename ):   # read Mantid-reduced IN16B QENS data
	class out:
		pass
	f = h5py.File( filename, 'r' )
	out.x      = np.transpose( np.array( f['/mantid_workspace_1/workspace/axis1/']  ) )
	out.y      = np.transpose( np.array( f['/mantid_workspace_1/workspace/axis2/']  ) )
	out.sqw    = np.transpose( np.array( f['/mantid_workspace_1/workspace/values/'] ) )
	out.dsqw   = np.transpose( np.array( f['/mantid_workspace_1/workspace/errors/'] ) )
	out.hw     = ( out.x[0:-1] + out.x[1:] ) / 2.0 * 1e3; # from Mantid bins to channel center
	wavelength = np.array( f['/mantid_workspace_1/logs/wavelength/value/'] )
	out.q      = 4. * np.pi / wavelength * np.sin( out.y * np.pi / 360.0 );
	return( out )

# ---------------------------------------------------------------------------------
def loadIN16B_Q( filename ):   # read Mantid-reduced IN16B QENS data, alreay converted to Q
	class out:
		pass
	f = h5py.File( filename, 'r' )
	out.x    = np.transpose( np.array( f['/mantid_workspace_1/workspace/axis1/']  ) )
	out.q    = np.transpose( np.array( f['/mantid_workspace_1/workspace/axis2/']  ) )
	out.sqw  = np.transpose( np.array( f['/mantid_workspace_1/workspace/values/'] ) )
	out.dsqw = np.transpose( np.array( f['/mantid_workspace_1/workspace/errors/'] ) )
	out.hw   = ( out.x[0:-1] + out.x[1:] ) / 2.0 * 1e3; # from Mantid bins to channel center
	return( out )

# in-plane q for tests:
#	out.q = 4. * np.pi / 6.271 * np.sin( np.concatenate( (np.array([ 10.94, 16.8 ]), np.array( np.linspace(33.1,150.1,16) ) ) ) * np.pi / 360.0 )


# ---------------------------------------------------------------------------------
def loadIN5( filename ):   # read lamp-reduced IN5 QENS data
	class out:
		pass
	f = h5py.File( filename, 'r' )
	out.hw   = np.transpose( np.array( f['/entry1/data1/X/']  ) ) * 1e3 # convert to micro-eV
	out.q    = np.transpose( np.array( f['/entry1/data1/Y/']  ) )
	out.sqw  = np.transpose( np.array( f['/entry1/data1/DATA/'] ) )
	out.dsqw = np.transpose( np.array( f['/entry1/data1/errors/'] ) )
	return( out )

# ---------------------------------------------------------------------------------
def saveIN16B( filename, d ):    # save data from the script to HDF5
	f    = h5py.File( filename, 'w' )
	x    = f.create_dataset( '/mantid_workspace_1/workspace/axis1/',  data = d.hw * 1e-3  ) # convert to meV
	y    = f.create_dataset( '/mantid_workspace_1/workspace/axis2/',  data = d.q    )
	sqw  = f.create_dataset( '/mantid_workspace_1/workspace/values/', data = np.transpose( d.sqw )  ) # aspect ratio as in Mantid
	dsqw = f.create_dataset( '/mantid_workspace_1/workspace/errors/', data = np.transpose( d.dsqw ) )
	f.close()

# ---------------------------------------------------------------------------------
# ----------- RESOLUTION FUNCTIONS ------------------------------------------------
# ---------------------------------------------------------------------------------
def resolution_model_IN16B( x, b, a0, sig, x0 ):
	return a0 * gaussian( x, sig, x0 ) + b

# ---------------------------------------------------------------------------------
def resolution_model_IN5( x, b, a0, sig, *args ): 
# resolution a single Gaussian or a sum of Gaussians ( number determined from length of sig)
	import sys
	if np.size( sig ) == 1: # if single Gaussian, test if an energy-offset was provided
		x00 = list( args )
		x0  = np.asarray( x00 ).flatten()
		if not x0:
			x0 = 0 # if the center value has not been given by the calling function
		return a0 * gaussian( x, sig, x0 ) + b # single Gaussian centered at 0 + constant background
	else:
		x0 = np.asarray( args ).flatten()
		y  = np.zeros( np.shape( x ) )
		for l in range( 0, len( sig ) ):
			y = y + a0[l] * gaussian( x, sig[l], x0[l] )  
		return( y + b ) # return sum of Gaussians + constant background

# ----------------------------------------------------------------------------------------------------------------------------
def resolutionIN16B( q, hw, sqw, dsqw ):  # fit resolution function of IN16B with a single Gaussian
	from scipy.optimize import curve_fit
	qn = range( 0, len(q) )
	class out:
		pass	
	# remove NaN values and other unreasonable points:
	hwd   = hw[  (hw>-4) & (hw<4) ]
	sqwd  = sqw[ (hw>-4) & (hw<4), : ]
	dsqwd = dsqw[(hw>-4) & (hw<4), : ]
	sqwc  = sqwd
	dsqwc = dsqwd
	sqwc[  np.isnan( sqwd )  ] = 0
	sqwc[  np.isinf( sqwd )  ] = 0
	sqwc[  sqwd < 5e-4*sqwd.max() ] = 0 # reasonable signal-to-noise range
	dsqwc[ np.isnan( sqwd )  ] = np.inf
	dsqwc[ np.isinf( sqwd )  ] = np.inf
	dsqwc[ np.isnan( dsqwd ) ] = np.inf
	dsqwc[  sqwd==0 ] = np.inf # no signal at all should have infinite error
	dsqwc[ dsqwd==0 ] = np.inf # "a measurement without error is nonsense"
	f0 = [  0.0, 1.0, 0.4,  0.0 ]
	l  = [  0.0, 0.0, 0.1, -0.5 ]
	u  = [  1.0, 1e3, 0.6,  0.5 ]
	for k in range( 0, len( qn ) ):
		y = sqwc[ :, k ]
		dy = dsqwc[ :, k ]
		popt, pcov = curve_fit( resolution_model_IN16B, hwd, y, p0=f0, bounds=(l,u), sigma=dy**2, ftol=1e-12, xtol=1e-12  ) 
		if k == 0:
			out.popt = popt
		else:
			out.popt = np.vstack( [ out.popt, popt ] )
	return( out )

# ----------------------------------------------------------------------------------------------------------------------------
def resolution_function( q, hw, sqw, dsqw, *args, **kwargs ):  # fit resolution function - IN5: 1 Gaussian, or with keyword 'LET' several Gaussians
# for IN5, the fit does not seem to work with bounds
# the length of the vector of start values defined the number of Gaussians for the resolution function
# lowest entry in this vector: constant background, then subsequently: (a_i, sigma_i, x0_i)
# default is 5 Gaussians
	from scipy.optimize import curve_fit
	import matplotlib # for plotting if "verbose" is on
	matplotlib.rc('text', usetex=True)
	matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
	matplotlib.use('ps') # used for saving figures to png file
	import matplotlib.pyplot as plt
	qmin = 0      # default assumption: use full q-range
	qmax = len(q) # default assumption: use full q-range
	class out:
		pass
	datapath = ''
	filenameextension = ''
	cryofit  = None
	rawvanad = None
	verbose  = 0
	useLET   = 0 # use several Gaussians for LET spectrometer resolution
	useIN5   = 0 # use several Gaussians for IN5 spectrometer resolution
	useIN16B = 0 # use two Gaussians for IN16B
	useBATS  = 0
	for key, value in kwargs.items():
		if key == 'BATS':
			if value:
				useBATS = 1
		if key == 'BATSPATH':
			datapath = value
		if key == 'FILENAME':
			filenameextension = value
		if key == 'CRYO':
			cryofit = value
		if key == 'RAWVANA':
			rawvanad = value
	erange  = 4 # energy range for fit
	f0 = [  0.0,    2e2,    35.0    ] # start values for single Gaussian centered at 0
	l  = [  0.0,    0.0,     0.0    ]
	u  = [  np.inf, np.inf, np.inf  ]
	if useBATS:
		f0 = [     0,  8e0,  2.5,    0, \
		              9e-3,   25,    0, \
		              4e-4,    6,   -1, \
		              4e-5,    6,    1, \
		              5e-5, 20.0,   -2 ]
		l  = [ -1e-5,    0,    1, -12, \
		                 0,    5, -12, \
		                 0,    5, -15, \
		                 0,    5, -15, \
		                 0,    5, -15 ]
		u  = [   1e2,  1e3,  30,  12, \
		               5e2,  60,  12, \
		               5e1,  60,  15, \
		               5e1,  60,  15, \
		               5e1,  40,  15 ]
		verbose = 1
		erange = 30
		print( 'Fitting BATS resolution ...' )
		print( 'Saving figures to: '+datapath )
	else:
		for count, thing in enumerate( args ):
			if (count == 0) & (thing == 'verbose'):
				verbose = 1
				print('Fitting resolution, verbose mode ...')
			if (count == 1 ) & (thing == 'IN16B'):
				useIN16B = 1
	# 			        b   a01  sig1  x01  a02  sig2  x02
				f0 = [  0.0, 1.0,    0.4,  0.0,    0.1,    0.8,  0.0,    0.05,    2.5,  0.0 ]
				l  = [  0.0, 0.0,    0.1, -0.5,    0.0,    0.1, -0.5,    0.0,     0.1, -0.5 ]
				u  = [  1.0, np.inf, 0.6,  0.5, np.inf, np.inf,  0.5, np.inf,  np.inf,  0.5 ]
			elif count == 1:
				erange = thing
				if verbose:
					print('Fit-range in energy:',-erange,erange)
			if (count == 2) & useIN16B:
				datapath = thing
				print( 'figures saving path: ', datapath )
			elif count == 2:
				qmin = thing
				if verbose:
					print('q_min=:',qmin)
			if count == 3:
				qmax = thing
				if verbose:
					print('q_max=:',qmax)
			if (count == 4 ) & (thing == 'LET'): # start values and boundaries for 5 Gaussians for LET resolution function
				useLET = 1
				f0 = [     0, 8e4,   25,    1, \
					      9e2,  250,   28, \
					      4e4,   52,  -57, \
					      4e2,   23,  -17, \
					      5e1,  200,  -50 ]
				l  = [ -1e-5,    0,   10, -100, \
					         0,   10, -100, \
					         0,   10, -150, \
					         0,   10, -150, \
					         0,   50, -150 ]
				u  = [   1e2,  1e9,  300,  100, \
					       5e7,  300,  100, \
					       5e6,  300,  150, \
					       5e6,  300,  150, \
					       5e6,  400,  150 ]
			if (count == 4 ) & (thing == 'IN5'): # start values and boundaries for several Gaussians for IN5 resolution function
				useIN5 = 1
			if count == 5:
				f0 = thing
				if verbose:
					print('using user-provided start values:',f0)
			if count == 6:
				l = thing
				if verbose:
					print('using user-provided lower limits:',l)
			if count == 7:
				u = thing
				if verbose:
					print('using user-provided upper limits:',u)
			if count == 8:
				datapath = thing
				print( 'figures saving path: ', datapath )
	qn = range( qmin, qmax )
	out.sqw   = np.zeros( np.shape( sqw ) )
	out.dsqw  = np.zeros( np.shape( dsqw ) )
	out.hw    = hw
	for k in qn:
		# remove NaN values and other unreasonable points:
		try:
			hw0   = hw[ :, k ]
		except:
			hw0   = hw
		hwd   = hw0[  (hw0 > -erange) & (hw0 < erange) ]
		sqwd  = sqw[  (hw0 > -erange) & (hw0 < erange), : ] * 1e-2
		dsqwd = dsqw[ (hw0 > -erange) & (hw0 < erange), : ] * 1e-2
		sqwc  = sqwd.copy()
		dsqwc = dsqwd.copy()
		sqwc[  np.isnan( sqwd )  ] = 0
		sqwc[  np.isinf( sqwd )  ] = 0
		sqwc[  sqwd < 1e-6*sqwd.max() ] = 0 # reasonable signal-to-noise range
		dsqwc[ np.isnan( sqwd )  ] = np.inf
		dsqwc[ np.isinf( sqwd )  ] = np.inf
		dsqwc[ np.isnan( dsqwd ) ] = np.inf
		dsqwc[  sqwd==0 ] = np.inf # no signal at all should have infinite error
		dsqwc[ dsqwd==0 ] = np.inf # "a measurement without error is nonsense"
		y  =  sqwc[ :, k ]
		dy = dsqwc[ :, k ]
		if useLET:
			print( "Fitting LET resolution, q=", q[k], "/A" )
			if k > qmin:
				f0 = popt # iteratively use previous q fit result as start values
			popt, pcov = curve_fit( lambda x, *p: wrapper_resolution_func( x, len( f0 ), p ), hwd, y, p0=f0, sigma=dy, bounds=(l,u), ftol=1e-12, xtol=1e-12 )
		elif useBATS:
			print( "Fitting BATS resolution, q=", q[k], "/A" )
			if k > qmin:
				f0 = popt # iteratively use previous q fit result as start values
#			print( 'x_min', x.min(), 'x_max', x.max()  )
			popt, pcov = curve_fit( lambda x, *p: wrapper_resolution_func( x, len( f0 ), p ), hwd, y, p0=f0, sigma=dy, bounds=(l,u), ftol=1e-12, xtol=1e-12, max_nfev=1e6 )
		elif useIN5: # first fit single Gaussian, use background, width, and offset of first Gaussian to fix for next iteration
			print( "Fitting IN5 resolution, q=", q[k], "/A" )
			f01 = [ 0,     8e1,        34,     0 ] # start values for single Gauss
			l01 = [ 0,      0,          0,  -100 ] # bounds for single Gauss
			u01 = [ np.inf, np.inf, np.inf,  100 ]
			popt0, pcov0 = curve_fit( lambda x, *p: wrapper_resolution_func( x, len( f01 ), p ), hwd, y, p0=f01, sigma=dy, bounds=(l01,u01), ftol=1e-15, xtol=1e-15 )
			print( 'Fit results 1st iteration:', popt0 )
			f0[0] = popt0[1] *0.6 # main Gauss amplitude start value for 2nd iterative fit
			l[0]  = popt0[1] *0.1 # main Gauss amplitude bounds for 2nd iterative fit
			u[0]  = popt0[1] *1.0
			print( 'Start values 2nd iteration:', f0 )
			print( 'Lower boundaries 2nd iteration:', l )
			print( 'Upper boundaries 2nd iteration:', u )
			print( 'Fixed values 2nd iteration:', popt0[0], popt0[2], popt0[3] )
			popt1, pcov1 = curve_fit( lambda x, *p: wrapper_1Gauss_fix_resolution_func( x, popt0[0], popt0[2], popt0[3], len( f0 ) , p ), hwd, y, p0=f0, sigma=dy, bounds=(l,u), ftol=1e-15, xtol=1e-15, gtol=1e-12, method='dogbox', jac='3-point', max_nfev=1e6 )
			print( '2nd iteration results: ', popt1 )
			print( '2nd iteration amplitudes 2nd Gauss and higher: ', popt1[1::3] )
			print( '2nd iteration widths     2nd Gauss and higher: ', popt1[2::3] )
			print( '2nd iteration offsets    2nd Gauss and higher: ', popt1[3::3] )
			popt = np.array( [ popt0[0], popt1[0], popt0[2], popt0[3] ] )
			print( 'Self-found number of Gaussians:', len(popt1[1::3]) )
			for si in range( 0, len(popt1[1::3]) ):
				popt = np.hstack( [ popt, popt1[ si*3+1 ], popt1[ si*3+2 ], popt1[ si*3+3 ] ] )
			pcov = np.zeros( [ len(popt), len(popt) ] ) # preliminary
		elif useIN16B:
			print( "Fitting IN16B resolution, q=", q[k], "/A" )
			popt, pcov = curve_fit( lambda x, *p: wrapper_resolution_func( x, len( f0 ), p ), hwd, y, p0=f0, sigma=dy, bounds=(l,u), ftol=1e-12, xtol=1e-12, max_nfev=1e8 )
		else:   # one Gaussian centered at 0 (e.g. IN5):
			popt, pcov = curve_fit( resolution_model_IN5, hwd, y, p0=f0, sigma=dy ) # , bounds=(l,u), ftol=1e-14, xtol=1e-14
			popt = np.hstack( [popt, 0 ] ) # add the center energy entry for IN5 ( not fitted, Gauss fixed at energy center ) 
		if k == qmin:
			out.popt = popt
		else:
			out.popt = np.vstack( [ out.popt, popt ] )
		if verbose:			
			print( "Resolution fit results: ", popt, pcov[0][0] )
			print( "building plot ..." )
			plt.errorbar( hwd, y, yerr=dy, marker='o', color='b', linestyle='None' )
			if ( useLET | useIN16B | useIN5 | useBATS ):
				colors = [ 'r', 'c', 'b', 'y', 'm', 'k', 'g', 'r', 'b' ]
				if cryofit is not None:
					for nbg in range( 0, len(popt[1::3]) ):
						print( "Plotting Gaussian no."+str(nbg)+" with color "+colors[nbg] )
						plt.plot( np.linspace( -erange, erange, 100 ), resolution_model_IN5( np.linspace( -erange, erange, 100 ), popt[0], popt[nbg*3+1], popt[nbg*3+2], popt[nbg*3+3] ), color=colors[nbg] )
				else:
					for nbg in range( 0, 1 ):
						print( "Plotting Gaussian no."+str(nbg)+" with color "+colors[nbg] )
						plt.plot( np.linspace( -erange, erange, 100 ), resolution_model_IN5( np.linspace( -erange, erange, 100 ), popt[0], popt[nbg*3+1], popt[nbg*3+2], popt[nbg*3+3] ), color=colors[nbg] )
				plt.plot( np.linspace( -erange, erange, 100 ), resolution_model_IN5( np.linspace( -erange, erange, 100 ), popt[0], popt[1::3], popt[2::3], popt[3::3] ), color='g' )
				out.sqw[ :, k ] = resolution_model_IN5( hw, popt[0], popt[1::3], popt[2::3], popt[3::3] )
			else:
				plt.plot( np.linspace( -erange, erange, 100 ), resolution_model_IN5( np.linspace( -erange, erange, 100 ), popt[0], popt[1], popt[2] ), color='g' )
			if cryofit is not None:
				plt.plot( cryofit.hw, cryofit.sqw[:,k], color='k' )
			if rawvanad is not None:
				plt.errorbar( rawvanad.hw, rawvanad.sqw[:,k], yerr=rawvanad.dsqw[:,k], marker='o', color='lightblue', linestyle='None' )
			plt.xlim( [ -erange, erange ] )
#			try:
			if cryofit is not None:
				plt.ylim( [ 0.5 * np.array( [ resolution_model_IN5( np.linspace( -erange, erange, 100 ), popt[0], popt[1::3], popt[2::3], popt[3::3] ).min(), 1e-9 ] ).max(),
					    1.5 * np.array( [ cryofit.sqw[:,k].max(), resolution_model_IN5( np.linspace( -erange, erange, 100 ), popt[0], popt[1::3], popt[2::3], popt[3::3] ).max() ] ).max() ] )
			else:
				plt.ylim( [ 0.5 * np.array( [ resolution_model_IN5( np.linspace( -20, 20, 100 ), popt[0], popt[1::3], popt[2::3], popt[3::3] ).min(), 1e-9 ] ).max(),
					    1.1 * resolution_model_IN5( np.linspace( -erange, erange, 100 ), popt[0], popt[1::3], popt[2::3], popt[3::3] ).max() ] )
				plt.xlim( [ -20, 20 ] )
#			except:
#				pass
			plt.yscale('log')
			plt.xlabel(r'$\hbar\omega$ [$\mu$eV]',fontweight='bold')
			plt.ylabel(r'Intensity [arb. units]',fontweight='bold')
			plt.tight_layout()
			try:		
				plt.savefig(datapath+'Fig_python_resolution'+filenameextension+'_fit_q'+str(k)+'.eps', format='eps' )
				plt.clf()
			except:
				print( "Could not display the plot" )
				print( datapath+'Fig_python_resolution'+filenameextension+'_fit_q'+str(k)+'.png' )
	return( out )


# --------------------------------------------------------------------------------
# ----------- OTHER HELPER FUNCTIONS ---------------------------------------------
# ---------------------------------------------------------------------------------
def subtract_EC( sample, can, *arg ): # subtract the empty can, apply an optional scaling factor
	import copy
	class out:
		pass 
	if len( arg ) == 1:
		sf = arg # scaling factor for EC subtraction
	else:
		sf = 1   # in case no scaling factor is provided
	sample_sqw  = copy.deepcopy( sample.sqw )
	can_sqw     = copy.deepcopy( can.sqw )
	sample_dsqw = copy.deepcopy( sample.dsqw )
	can_dsqw    = copy.deepcopy( can.dsqw )
	Delta_sqw   = sample_sqw - sf * can_sqw
	dDelta_sqw  = np.sqrt( sample_dsqw**2 + ( sf * can_dsqw )**2 )
	out.hw      = sample.hw
	out.q       = sample.q
	out.sqw     = Delta_sqw
	out.dsqw    = dDelta_sqw 
	return( out )

# -------------------------------------------------------------------------------------------------
def goodness_of_fit( x, y, p, q, qn, r, sigma, **kwargs ):
	yf  = model_sqw_interface_parser( x, p, q, qn, r, **kwargs )
	r2  = np.sum( ( y - yf )**2 / sigma **2 ) # weighted residuals
	return r2 / ( len( y ) - len( p ) )

def L1Err( x, y, p, q, qn, r, sigma, **kwargs ):
	""" The L1 loss function as error estimate """
	yf = model_sqw_interface_parser( x, p, q, qn, r, **kwargs )
	L1 = np.sum( abs(y - yf)) # L1 error is the difference between true value, estimated value in absolute numbers
	return L1/len(y) # Division by number of data points, as per definition (NOT degrees of freedom)

def L2Err( x, y, p, q, qn, r, sigma, **kwargs ):
        """ The L2 loss function as error estimate """
        yf = model_sqw_interface_parser( x, p, q, qn, r, **kwargs )
        L2 = np.sum( (y - yf)**2 ) # L2 error is the squared difference
        return L2/len(y) # Division by number of data points, as per definition (NOT degrees of freedom)

def BayesSivia(x, y, p, q, qn, r, sigma, **kwargs):
	""" Implementation of the Bayesian Analysis procedure suggested by Sivia et al. (1992) \
	for model selection."""
	yf = model_sqw_interface_parser( x, p, q, qn, r, **kwargs )
	r2  = np.sum( ( y - yf )**2 / sigma **2 ) # weighted residuals
	chi_square = r2 / ( len( y ) - len( p ) ) # Chi-square
	# Compute maximum amplitude and width of any Lorentzian
	
	# Compute the number of Lorentzians used

	# Compute determinant of the Hessian matrix, with the optimal parameters found

	# Compute the probability
	
