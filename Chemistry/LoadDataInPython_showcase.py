# exec(open("./qens_script_example_from_other_experiment.py").read())
import utilities_qens as utilq
import numpy as np
import copy
import sys # To exit whenever


Dunits   =   10 / 6.58212 # assuming microeV Energy axis, transforms Lorentzian HWHM to A^2/ns

# fill in only these two parameters: datapath and sample number to use:
datapath    = '/Users/ericfagerberg/NeutronData/exp_8-04-868_in16b/processed_PaalmanPingsCorrections/CorrectionsFiles/'
datapath2   = '/Users/ericfagerberg/NeutronData/exp_8-04-813_182_in16b/processed_PaalmanPingsCorrections/CorrectedData/' # For the 200 mg/ml data.




# --------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------


# use temperatures of 280 K, looks better at low temp.
protocol = ([ \
[ 'Hist_050_280K', 'corr_268916:268947', '280', '050', '0.03', 'Hist_050_280' ],
[ 'Hist_100_280K', 'corr_268846:268863', '280', '100', '0.06', 'Hist_100_280' ],
[ 'Hist_150_280K', 'corr_268950:268966', '280', '150', '0.09', 'Hist_150_280' ],
[ 'Hist_200_150mM_280K', 'corr_232428:232436', '280', '200', '0.12', 'Hist_200_150mM_280' ],
])

# Just load the data
sample50data = utilq.loadIN16B( datapath+protocol[0][1] ) # 50 mg/ml
sample100data = utilq.loadIN16B( datapath+protocol[1][1] ) # 100 mg/ml
sample150data = utilq.loadIN16B( datapath+protocol[2][1] ) # 150 mg/ml
sample200data = utilq.loadIN16B( datapath2+protocol[3][1] ) # 200 mg/ml

AllData = [sample50data, sample100data, sample150data, sample200data]

# These should be the same for all samples
hw  = sample50data.hw.copy()
q   = 57.1597*sample50data.q.copy()
qn  = range(0,len(q))
ql  = len(qn)


# Skip solvent, just plot the protein+D2O signal, without solvent subtraction

SQW = [] # Add the sqw data here
DSQW = [] # Add error here

# Clean the data
# Make a loop


for sample in AllData:
	dsqw  = copy.deepcopy( sample.dsqw )
	sqwc  = copy.deepcopy( sample.sqw )
	dsqwc = copy.deepcopy( sample.dsqw )
	sqwc[  np.isnan( sample.sqw )  ] = 0
	sqwc[  np.isinf( sample.sqw )  ] = 0
	sqwc[  sample.sqw < 5e-4*sample.sqw.max() ] = 0 # reasonable signal-to-noise range
	dsqwc[ np.isnan( sample.sqw )  ] = np.inf
	dsqwc[ np.isinf( sample.sqw )  ] = np.inf
	dsqwc[ np.isnan( sample.dsqw ) ] = np.inf
	dsqwc[  sample.sqw==0 ] = np.inf # no signal at all should have infinite error
	dsqwc[ sample.dsqw==0 ] = np.inf # "a measurement without error is nonsense"
	SQW.append(sqwc)
	DSQW.append(dsqwc)

arr_SQW = np.array(SQW)
arr_DSQW = np.array(DSQW)


### Do whatever. ###

# Just to show how the data data is built:



print('Slice of 2nd concentration (i.e. 100 mg/ml), 3rd q-value, the first half of the energies (0-512).')
slice100gml_q3_hw = arr_SQW[1][:512, 2]
print(slice100gml_q3_hw)




