This repo will be used to host methods to convert a set of protein structure coordinates from Cartesian to Bond-Angle-Torsion (BAT) 


The converter are in Python (stored in _./python_) and CPP (stored in _./cpp_)


Python version was tested under Python 3.10. The conda environment used to build the Python version can be found as bat.yml under _./python_. Numba accelerated version is available for efficient processing of large number of input files (i.e. a MD trajectory). Validation was done against MDAnalysis functions. 


CPP version is still under development, and will be updated 