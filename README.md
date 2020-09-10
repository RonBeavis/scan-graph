# scan-graph

*scan_graph* is a simple Python 3.0 project that shows how to extract a particular MS/MS scan out of a Thermo .raw file and generate an annotated graph for that scan using an assigned peptide sequence. This action reads the .raw file directly and it does not require performing an intermediate file conversion to some open text format, e.g. mzML or MGF.

There are three files:

1. scan_graph.py - contains methods used to extract a scan from a .raw file and to generate lists of annotation peaks

2. scan_graph_test.py - a demonstration of how to use the methods in scan_graph.py

3. 9102-NGKITSIVKDSSAARNG.png - a sample PNG generated using scan_graph_test.py

Access to .raw files uses *pymsfilereader*, so you must install this package using the instructions at https://github.com/frallain/pymsfilereader. This package is only available for Windows, because accessing the files requires using a COM interface supplied by Thermo.

To run *scan_graph_test.py*, you will need an accessible .raw data file. In the code, a data file is specified in the path variable:

path = 'PXD018998\\01_001815W_KLH_2.raw'

1. Place *scan_graph.py* and *scan_graph_test.py* in a new directory. 
2. Create the subdirectory PXD018998. 
3. Grab the appropriate .raw from from ftp://massive.ucsd.edu/MSV000085375/raw/01_001815W_KLH_2.raw
4. Place 01_001815W_KLH_2.raw into PXD018998
5. From the command line, change into the directory containing the scripts and run *> python scan_graph_test.py*
