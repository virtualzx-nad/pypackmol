#!/usr/bin/env python
"""
PyPackmol is a Python wrapper for packmol program that provide simple packing of molecules. 
Note that packmol needs to be installed separately.  
Molecules can be given in terms of input files (xyz, etc), or SMILES strings.  
Open Babel is required to do the conversion if the input format is not xyz.
"""

__title__     = "pypackmol"
__author__    = "Xiaolei Zhu"
__copyright__ = "Martinez Group, Stanford University, CA, USA, Planet Earth"
__maintainer__= "Xiaolei Zhu"
__email__     = "virtualzx@gmail.com"
__version__   = "0.0.1"
__status__    = "Prototype"

import tempfile
import os
import numpy as np

# If openbabel is imported, it is possible to convert different file types,
# construct molecules from SMILES and/or perform force field optimizations.
# Note that the default is to use force field optimization for SMILES but 
# keep the geometries for file input.
try :
    import pybel as pb
    _has_openbabel=True
except :
    _has_openbabel=False

    
def convert_to_xyz(input_data,input_format='auto') :
    """Converts an input file or SMILES string into xyz format
    
    When input_format is 'auto' which is the default, the function will try 
    to determine the type of input file.
    In such event, if the input string does not correspond to an existing 
    file that can be opened, the string is instead taken to be a SMILES.
    The return is a temp file which contains xyz format molecule definition.
    
    """
    if input_format == 'auto' :
        try :
            with open(input_data) as input_file :
                input_format = input_data.split('.')[-1]
        except :
            input_format ='smi'
            
    # Create a temp file which saves the output
    output_file = tempfile.NamedTemporaryFile(mode='w+')
    
    # If the file is already an xyz, there is no need to convert
    if input_format == 'xyz' :
        with open(input_data) as in_file :
            output_file.write(in_file.read())
    
    # The file is not xyz format.  Conversion is needed. 
    # If do not have openbabel, then we are pretty screwed here.
    if not _has_openbabel :
        raise Exception["Cannot import Open Babel for format conversion. " +
                       "Please use xyz format instead."]
    
    if input_format=='smi' or imput_format=='smiles' :
        mol=pb.readstring(input_format,input_data)
        mol.make3D()
        output_file.write(mol.write('xyz'))
    output_file.seek(0)
    return output_file



class Packmol :
    """Wrapper class for the packmol molecule packer program
    """
    
    
    class PmMolecule :
        """Packmol molecule class which contains info of molecule file name, 
        count, packing region, etc"""
        def __init__(self,  xyzfile,  count=1):
            pass
        
    def __init__(self, region_type='sphere',  dimensions=10.0 , tolerance=3.0 ,  **kwargs ):
        """Construct packmol run by providing simulation parameters.  """
        _dimensions  = dimensions
        _tolerance   = tolerance
        _region_type = region_type
        _packmol_arguments = kwargs        
        


