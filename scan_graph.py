#
# Copyright © 2020 Ronald C. Beavis
# Licensed under Apache License, Version 2.0, January 2004
#

# you must install pymsfilereader using the instructions at
# https://github.com/frallain/pymsfilereader

from pymsfilereader import MSFileReader
import re

# some useful constant values
isotopes = {	'p' : 1.007276,
		'H' : 1.007825,
		'C' : 12.0,
		'N' : 14.003074,
		'O': 15.994915 }

a_to_m = {	"A":71.037114,
		"R":156.101111,
		"B":114.042927,
		"N":114.042927,
		"D":115.026943,
		"C":103.009185,
		"E":129.042593,
		"Q":128.058578,
		"Z":128.058578,
		"G":57.021464,
		"O":237.147727,
		"H":137.058912,
		"I":113.084064,
		"J":113.084064,
		"L":113.084064,
		"K":128.094963,
		"M":131.040485,
		"F":147.068414,
		"P":97.052764,
		"S":87.032028,
		"T":101.047679,
		"U":150.95363,
		"W":186.079313,
		"Y":163.06332,
		"V":99.068414 }

# calculate the neutral masses of Y-type fragments
def GetYs(_pep):
	dY = 2*isotopes['H'] + isotopes['O']
	# reverse the peptide sequence
	ss = list(_pep['seq'][::-1])
	ms = []
	l = len(_pep['seq'])
	# create an array of residue masses
	for i,aa in enumerate(ss):
		if l-i in _pep['mods']:
			ms.append(a_to_m[aa]+_pep['mods'][l-i])
		else:
			ms.append(a_to_m[aa])
	rvs= []
	cvalue = dY
	# create an array of fragment masses
	for m in ms:
		cvalue += m
		rvs.append(cvalue)
	return rvs

# calculate the neutral masses of B-type fragments
def GetBs(_pep):
	dB = 0.0
	ss = list(_pep['seq'])
	ms = []
	# create an array of residue masses
	for i,aa in enumerate(ss):
		if i+1 in _pep['mods']:
			ms.append(a_to_m[aa]+_pep['mods'][i+1])
		else:
			ms.append(a_to_m[aa])
	rvs= []
	cvalue = dB
	# create an array of fragment masses
	for m in ms:
		cvalue += m
		rvs.append(cvalue)
	return rvs

# take a set of neutral masses (_vs) and convert them into protonated ions for a given charge _z		
def GetCharge(_vs,_z):
	rvs = []
	proton = isotopes['p']
	z = float(_z)
	for v in _vs:
		rvs.append((v+z*proton)/z)
	return rvs

# compare a set of ions (_vals) to a set of peaks (_expts) and return the union
def GetValues(_vals,_expts,_tol):
    rvs = [[],[]]
    for v in _vals:
        for i,x in enumerate(_expts[0]):
            if abs(v-x) <= _tol:
                rvs[0].append(x)
                rvs[1].append(_expts[1][i])
    return rvs

# retrieve an MS/MS spectrum, based on a RAW file path name and a scan number
def GetSpectrum(_path,_scan):
	info = {}
	rawfile = MSFileReader(path)
	info['Version'] = rawfile.Version()
	info['GetFileName'] = rawfile.GetFileName()
	info['RTFromScanNum'] = rawfile.RTFromScanNum(_scan)
	info['ScanNumFromRT'] = rawfile.ScanNumFromRT(rawfile.RTFromScanNum(_scan))
	info['GetFilterForScanNum'] = rawfile.GetFilterForScanNum(_scan)
	info['GetMSOrderForScanNum'] = rawfile.GetMSOrderForScanNum(_scan)
	spectrum = rawfile.GetSummedMassSpectrum([_scan],centroidResult=True)[0]
	rawfile.Close()
	return ([spectrum[0],spectrum[1]],info)
