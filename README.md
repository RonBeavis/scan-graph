# scan-graph

*scan_graph* is a simple Python 3.0 project that shows how to pull a particular MS/MS scan out of a Thermo .raw file and generate an annotated graph for that scan, given an assigned peptide sequence.

There are three files:

1. scan_graph.py - contains methods used to extract a scan from a .raw file and to generate lists of annotation peaks

2. scan_graph_test.py - a demonstration of how to use the methods in scan_graph.py

3. 9102-NGKITSIVKDSSAARNG.png - a sample PNG generated using scan_graph_test.py

Access to .raw files uses *pymsfilereader*, so you must install this package using the instructions at https://github.com/frallain/pymsfilereader. This package is only available for Windows, because accessing the files requires using a COM interface supplied by Thermo.
