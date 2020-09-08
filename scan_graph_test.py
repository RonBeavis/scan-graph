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

# specify the ion types (and plot colors) to use
	itypes = {'b-NH3':(0.3,0.6,0.9,1.0),
			'b-H2O':(0.6,0.6,0.9,1.0),
			'b-ion':(0.1,0.1,0.9,1.0),
			'y-NH3':(0.9,0.6,0.3,1.0),
			'y-H2O':(0.9,0.6,0.6,1.0),
			'y-ion':(0.9,0.1,0.1,1.0)
			}
# get the ion information
	print('type\tz\tm/z\tintensity')
	for it in itypes:
		if it[0] == 'b':
			bvals = GetBs(peptide,it[1:])
		else:
			bvals = GetYs(peptide,it[1:])
		# reset the peak lists
		ps = [[],[]]
		# iterate through allowed charge states
		for z in range(1,peptide['z']+1):
			# add in charge information
			zvals = GetCharge(bvals,z)
			# get matched peaks
			pvals = GetValues(zvals,expt,peptide['tol'])
			# update the peak lists
			ps[0] += pvals[0]
			ps[1] += pvals[1]
			# create a table entry (mainly for debugging purposes)
			for i,p in enumerate(pvals[0]):
				print('%s\t%i\t%.3f\t%.0f' % (it,z,p,pvals[1][i]))
		peaks[it] = ps

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
	ax.bar(peaks['expt'][0],peaks['expt'][1],color=(0.6,0.6,0.6,0.6),width=2,label='unmatched')
	for it in itypes:
		ax.bar(peaks[it][0],peaks[it][1],color=itypes[it],width=4,label=it)
	ax.legend()
	plt.savefig('%i-%s.png' % (scan,peptide['seq']))
	plt.show()

if __name__ == "__main__":
    main()

