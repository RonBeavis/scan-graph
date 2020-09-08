#
# Copyright © 2019 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
from scan_graph import GetSpectrum,GetBs,GetCharge,GetYs,GetValues
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl

def main():
# provide information for retrieving the spectrum
	scan = 9101 + 1
	peptide = {	'tol':0.2,	#fragment ion tolerance (in Da)
			'z':3,		#maximum fragment ion charge to consider
			'seq': 'NGKITSIVKDSSAARNG', #peptide sequence
			'mods':{1:0.984} #peptide sequence modifications using position:Da pairs
		}
	path = 'PXD018998\\01_001815W_KLH_2.raw'	#path to the spectrum file

# retrieve the spectrum and some text information
	(expt,info) = GetSpectrum(path,scan)
# rescale the spectrum to run from 0 to 100
	iscale = max(expt[1])
	expt[1] = [100.0*x/iscale for x in expt[1]]
	peaks = {'expt':expt}

# process and plot b ions
	bvals = GetBs(peptide)
	ps = [[],[]]
	for z in range(1,peptide['z']+1):
		zvals = GetCharge(bvals,z)
		pvals = GetValues(zvals,expt,peptide['tol'])
		ps[0] += pvals[0]
		ps[1] += pvals[1]
		for i,p in enumerate(pvals[0]):
			print('B\t%i\t%.3f\t%.0f' % (z,p,pvals[1][i]))
	peaks['b-ions'] = ps
# process and plot y ions
	yvals = GetYs(peptide)
	ps = [[],[]]
	for z in range(1,peptide['z']+1):
		zvals = GetCharge(yvals,z)
		pvals = GetValues(zvals,expt,peptide['tol'])
		ps[0] += pvals[0]
		ps[1] += pvals[1]
		for i,p in enumerate(pvals[0]):
			print('Y\t%i\t%.3f\t%.0f' % (z,p,pvals[1][i]))
	peaks['y-ions'] = ps

# display the plot
	fig = plt.figure(figsize=(10, 5), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel('m/z')
	ax.set_ylabel('intensity')
	ax.set_title('#%i, %s' % (scan,peptide['seq']))
	ax.bar(peaks['expt'][0],peaks['expt'][1],color=(0.3,0.3,0.3,.5),width=2,label='unmatched')
	ax.bar(peaks['b-ions'][0],peaks['b-ions'][1],color=(0.1,0.1,0.9,1.0),width=4,label='b-ion')
	ax.bar(peaks['y-ions'][0],peaks['y-ions'][1],color=(0.9,0.1,0.1,1.0),width=4,label='y-ion')
	ax.legend()
	plt.show()

if __name__ == "__main__":
    main()

