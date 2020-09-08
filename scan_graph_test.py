import scan_graph
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
	path = 'GPM64220031057\\01_001815W_KLH_2.raw'	#path to the spectrum file

# retrieve the spectrum and some text information
	(expt,info) = GetSpectrum(path,scan)

# start constructing the graph
	plt.xlabel('m/z')
	plt.ylabel('intensity')
	plt.title('#%i, %s' % (scan,peptide['seq']))
	plt.gcf().set_size_inches(10,5)
	plt.bar(expt[0],expt[1],color=(0.3,0.3,0.3,.5),zorder=1)

# rescale the spectrum to run from 0 to 100
	iscale = max(expt[1])
	expt[1] = [100.0*x/iscale for x in expt[1]]

# process and plot b ions
	bvals = GetBs(peptide)
	for z in range(1,peptide['z']+1):
	    zvals = GetCharge(bvals,z)
	    pvals = GetValues(zvals,expt,peptide['tol'])
	    plt.bar(pvals[0],pvals[1],color=(0.1,0.1,0.9,1.0),width=4,zorder=10)
	    for i,p in enumerate(pvals[0]):
		print('B\t%i\t%.3f\t%.0f' % (z,p,pvals[1][i]))

# process and plot y ions
	yvals = GetYs(peptide)
	for z in range(1,peptide['z']+1):
	    zvals = GetCharge(yvals,z)
	    pvals = GetValues(zvals,expt,peptide['tol'])
	    plt.bar(pvals[0],pvals[1],color=(0.9,0.1,0.1,1.0),width=4,zorder=10)
	    for i,p in enumerate(pvals[0]):
		print('Y\t%i\t%.3f\t%.0f' % (z,p,pvals[1][i]))

# display the plot
	plt.show()

if __name__ ==  == "__main__":
    main()

