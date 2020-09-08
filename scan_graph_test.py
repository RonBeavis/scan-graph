#
# Copyright © 2020 Ronald C. Beavis
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
	if False:
		scan = 10234 + 1
		peptide = {	'tol':0.2,	#fragment ion tolerance (in Da)
				'z':3,		#maximum fragment ion charge to consider
				'seq': 'SAADEVDGLGVARPHYGSVLDNER', #peptide sequence
				'mods':{1:42.011} #peptide sequence modifications using position:Da pairs
			}
		path = 'PXD000865\\00576_E01_P004283_B0E_A00_R1.raw'
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
	
	bvals = GetBs(peptide,'-H2O')
	ps = [[],[]]
	for z in range(1,peptide['z']+1):
		zvals = GetCharge(bvals,z)
		pvals = GetValues(zvals,expt,peptide['tol'])
		ps[0] += pvals[0]
		ps[1] += pvals[1]
		for i,p in enumerate(pvals[0]):
			print('B\t%i\t%.3f\t%.0f' % (z,p,pvals[1][i]))
	peaks['b-H20'] = ps
	
	bvals = GetBs(peptide,'-NH3')
	ps = [[],[]]
	for z in range(1,peptide['z']+1):
		zvals = GetCharge(bvals,z)
		pvals = GetValues(zvals,expt,peptide['tol'])
		ps[0] += pvals[0]
		ps[1] += pvals[1]
		for i,p in enumerate(pvals[0]):
			print('B\t%i\t%.3f\t%.0f' % (z,p,pvals[1][i]))
	peaks['b-NH3'] = ps

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

	yvals = GetYs(peptide,'-H2O')
	ps = [[],[]]
	for z in range(1,peptide['z']+1):
		zvals = GetCharge(yvals,z)
		pvals = GetValues(zvals,expt,peptide['tol'])
		ps[0] += pvals[0]
		ps[1] += pvals[1]
		for i,p in enumerate(pvals[0]):
			print('Y\t%i\t%.3f\t%.0f' % (z,p,pvals[1][i]))
	peaks['y-H20'] = ps

	yvals = GetYs(peptide,'-NH3')
	ps = [[],[]]
	for z in range(1,peptide['z']+1):
		zvals = GetCharge(yvals,z)
		pvals = GetValues(zvals,expt,peptide['tol'])
		ps[0] += pvals[0]
		ps[1] += pvals[1]
		for i,p in enumerate(pvals[0]):
			print('Y\t%i\t%.3f\t%.0f' % (z,p,pvals[1][i]))
	peaks['y-NH3'] = ps

# display the plot
	fig = plt.figure(figsize=(10, 5), dpi=100)
	ax = fig.add_subplot(111)
	ax.set_xlabel('m/z')
	ax.set_ylabel('intensity')
	title = '#%i, %s' % (scan,peptide['seq'])
	if len(peptide['mods']) > 0:
		title += '\n'
		for m in peptide['mods']:
			title += '%s%i+%.3f,' % (peptide['seq'][m-1],m,peptide['mods'][m])
		title = title[0:-1]
	ax.set_title(title)
	ax.set_xlim(0,1.05*max(peaks['expt'][0]))
	ax.bar(peaks['expt'][0],peaks['expt'][1],color=(0.6,0.6,0.6,.5),width=2,label='unmatched')
	ax.bar(peaks['b-NH3'][0],peaks['b-NH3'][1],color=(0.3,0.6,0.9,1.0),width=4,label='b-NH3')
	ax.bar(peaks['b-H20'][0],peaks['b-H20'][1],color=(0.6,0.6,0.9,1.0),width=4,label='b-H20')
	ax.bar(peaks['b-ions'][0],peaks['b-ions'][1],color=(0.1,0.1,0.9,1.0),width=4,label='b-ion')
	ax.bar(peaks['y-NH3'][0],peaks['y-NH3'][1],color=(0.9,0.6,0.3,1.0),width=4,label='y-NH3')
	ax.bar(peaks['y-H20'][0],peaks['y-H20'][1],color=(0.9,0.6,0.6,1.0),width=4,label='y-H20')
	ax.bar(peaks['y-ions'][0],peaks['y-ions'][1],color=(0.9,0.1,0.1,1.0),width=4,label='y-ion')
	ax.legend()
	plt.savefig('%i-%s.png' % (scan,peptide['seq']))
	plt.show()

if __name__ == "__main__":
    main()

