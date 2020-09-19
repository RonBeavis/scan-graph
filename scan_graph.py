#
# Copyright © 2020 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#
# you must install pymsfilereader using the instructions at
# https://github.com/frallain/pymsfilereader
#
# the most recent version of this file is available at
# https://github.com/RonBeavis/scan-graph/blob/master/scan_graph.py
#
# this file contains the class ScanGraph which is designed to be useful for extracting MS/MS data
# from a .raw file and creating lists of markup peaks illustrating the signals assigned to a particular
# peptide sequences
# 

from pymsfilereader import MSFileReader
import re

class ScanGraph:
# some useful constant values
	isotopes = {	
		'p' : 1.007276,
		'H' : 1.007825,
		'C' : 12.0,
		'N' : 14.003074,
		'O' : 15.994915,
		'S' : 31.972071,
		'Se' : 79.91652,
		'P' : 30.973762
	}

	a_to_m = {
		'A' : 71.037114,
		'R' : 156.101111,
		'B' : 114.042927,# interpret B as N
		'N' : 114.042927,
		'D' : 115.026943,
		'C' : 103.009185,
		'E' : 129.042593,
		'Q' : 128.058578,
		'Z' : 128.058578,# interpret Z as Q
		'G' : 57.021464,
		'O' : 237.147727,
		'H' : 137.058912,
		'I' : 113.084064,
		'J' : 113.084064,
		'L' : 113.084064,
		'K' : 128.094963,
		'M' : 131.040485,
		'F' : 147.068414,
		'P' : 97.052764,
		'S' : 87.032028,
		'T' : 101.047679,
		'U' : 150.95363,
		'W' : 186.079313,
		'Y' : 163.06332,
		'V' : 99.068414
	}
	# __init function for class
	def __init__(self):
		return
		
	# interpret a fragment neutral loss mass shift string
	def GetDelta(self,_delta):
		delta = 0;
		if _delta:
			# check for defined strings
			if _delta == '-NH3':
				delta = -1.0*(self.isotopes['N'] + 3*self.isotopes['H'])
			elif _delta == '-H2O':
				delta = -1.0*(2*self.isotopes['H'] + self.isotopes['O'])
			elif _delta == '-CH4SO':
				delta = -1.0*(self.isotopes['C'] + 4*self.isotopes['H'] + self.isotopes['S'] + self.isotopes['O'])
			elif _delta == '-H3PO4':
				delta = -1.0*(3*self.isotopes['H'] + self.isotopes['P'] + 4*self.isotopes['O'])
			elif _delta == '-H5PO5':
				delta = -1.0*(5*self.isotopes['H'] + self.isotopes['P'] + 5*self.isotopes['O'])
			else:
				# check to see if _delta is a number
				try:
					delta = float(_delta)
				except:
					delta = 0.0
		return delta
	
	# calculate the neutral masses of Y-type fragments
	def GetYs(self,_pep,_delta = None):
		delta = self.GetDelta(_delta)
		dY = 2*self.isotopes['H'] + self.isotopes['O']
		# reverse the peptide sequence
		ss = list(_pep['seq'][::-1])
		ms = []
		l = len(_pep['seq'])
		# create an array of residue masses
		for i,aa in enumerate(ss):
			if l-i in _pep['mods']:
				ms.append(self.a_to_m[aa]+_pep['mods'][l-i])
			else:
				ms.append(self.a_to_m[aa])
		rvs= []
		cvalue = dY + delta
		# create an array of fragment masses
		for m in ms:
			cvalue += m
			rvs.append(cvalue)
		return rvs[:-1]

	# calculate the neutral masses of parent ions and fragments
	def GetParents(self,_pep):
		difs = set([self.GetDelta('-H2O'),2.0*self.GetDelta('-H2O')])
		phospho = self.isotopes['H'] + 3*self.isotopes['O'] + self.isotopes['P']
		for m in _pep['mods']:
			if _pep['seq'][m-1] == 'M' and abs(_pep['mods'][m] - self.isotopes['O']) < 0.01:
				difs.add(self.GetDelta('-CH4SO'))
			elif _pep['seq'][m-1] == 'S' and abs(_pep['mods'][m] - phospho) < 0.01:
				difs.add(self.GetDelta('-H3PO4'))
				difs.append(self.GetDelta('-H5PO5'))
			elif _pep['seq'][m-1] == 'T' and abs(_pep['mods'][m] - phospho) < 0.01:
				difs.add(self.GetDelta('-H3PO4'))
				difs.add(self.GetDelta('-H5PO5'))
			elif _pep['seq'][m-1] == 'Y' and abs(_pep['mods'][m] - phospho) < 0.01:
				difs.add(self.GetDelta('-H3PO4'))
				difs.add(self.GetDelta('-H5PO5'))
		dY = 2*self.isotopes['H'] + self.isotopes['O']
		# reverse the peptide sequence
		ss = list(_pep['seq'][::-1])
		ms = []
		l = len(_pep['seq'])
		# create an array of residue masses
		for i,aa in enumerate(ss):
			if l-i in _pep['mods']:
				ms.append(self.a_to_m[aa]+_pep['mods'][l-i])
			else:
				ms.append(self.a_to_m[aa])
		rvs= []
		cvalue = dY
		# create an array of fragment masses
		for m in ms:
			cvalue += m
		rvs.append(cvalue)
		for d in difs:
			rvs.append(cvalue+d)
		return rvs

	# calculate the neutral masses of immonium ions
	def GetImmonium(self,_pep):
		delta = self.isotopes['C'] + self.isotopes['O']
		rvs = set()
		# create an array of fragment masses
		for i,aa in enumerate(_pep['seq']):
			if i+1 in _pep['mods']:
				rvs.add(self.a_to_m[aa]-delta+_pep['mods'][i+1])
			else:
				rvs.add(self.a_to_m[aa]-delta)
		return list(rvs)

	# calculate the neutral masses of B-type fragments
	def GetBs(self,_pep,_delta = None):
		delta = self.GetDelta(_delta)
		dB = 0.0
		ss = list(_pep['seq'])
		ms = []
		# create an array of residue masses
		for i,aa in enumerate(ss):
			if i+1 in _pep['mods']:
				ms.append(self.a_to_m[aa]+_pep['mods'][i+1])
			else:
				ms.append(self.a_to_m[aa])
		rvs= []
		cvalue = dB + delta
		# create an array of fragment masses
		for m in ms:
			cvalue += m
			rvs.append(cvalue)
		return rvs[:-1]

	# calculate the neutral masses of A-type fragments
	def GetAs(self,_pep,_delta = None):
		delta = self.GetDelta(_delta)
		dA = -1.0*(self.isotopes['C'] + self.isotopes['O'])
		ss = list(_pep['seq'])
		ms = []
		# create an array of residue masses
		for i,aa in enumerate(ss):
			if i+1 in _pep['mods']:
				ms.append(self.a_to_m[aa]+_pep['mods'][i+1])
			else:
				ms.append(self.a_to_m[aa])
		rvs= []
		cvalue = dA + delta
		# create an array of fragment masses
		for m in ms:
			cvalue += m
			rvs.append(cvalue)
		return rvs[:-1]

	# take a set of neutral masses (_vs) and convert them into protonated ions for a given charge _z		
	def GetCharge(self,_vs,_z):
		rvs = []
		proton = self.isotopes['p']
		z = float(_z)
		for v in _vs:
			rvs.append((v+z*proton)/z)
		return rvs

	# compare a set of ions (_vals) to a set of peaks (_expts) and return the intersection
	def GetValues(self,_vals,_expts,_tol):
		rvs = [[],[]]
		for v in _vals:
			for i,x in enumerate(_expts[0]):
				if abs(v-x) <= _tol:
					rvs[0].append(x)
					rvs[1].append(_expts[1][i])
		return rvs

	# retrieve an MS/MS spectrum, based on a RAW file path name and a scan number
	def GetSpectrum(self,_path,_scan):
		info = {}
		rawfile = MSFileReader(_path)
		info['Version'] = rawfile.Version()
		info['GetFileName'] = rawfile.GetFileName()
		info['RTFromScanNum'] = rawfile.RTFromScanNum(_scan)
		info['ScanNumFromRT'] = rawfile.ScanNumFromRT(rawfile.RTFromScanNum(_scan))
		info['GetFilterForScanNum'] = rawfile.GetFilterForScanNum(_scan)
		info['GetMSOrderForScanNum'] = rawfile.GetMSOrderForScanNum(_scan)
		spectrum = rawfile.GetSummedMassSpectrum([_scan],centroidResult=True)[0]
		rawfile.Close()
		return ([spectrum[0],spectrum[1]],info)

