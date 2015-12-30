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

PACKMOL_FORMATS = ['xyz','pdb','tinker','moldy']


def convert_to_xyz(input_data,input_format='auto') :
    """Convert input format unrecognized by packmol into xyz files.
    
    If the input molecular structure data can be recognized by packmol, the 
    function simply transfer the data to a temp file.  When the input is not
    recognized by packmol, the function use OpenBabel to convert into xyz 
    format, then save this to a temp file.
    
    Args:
        input_data: Filename or SMILES string to be converted
        input_format: Format of the input.  When the value is 'auto' which is 
            the default, the function will try to determine the type of input 
            file. In such event, if the input string does not correspond to an 
            existing file that can be opened, the string is instead taken to 
            be a SMILES string. 
    
    Returns:
        A named temporary file containing the molecular structure in a packmol 
        recognized format (xyz, pdb, tinker or moldy).  When conversion is 
        needed, the function always convert to xyz.  Note that since this is a 
        temp file, it will automatically be deleted when the file object is 
        unreferenced or destroyed.
    """
    input_format=input_format.lower()
    if input_format == 'auto' :
        try :
            with open(input_data) as input_file :
                input_format = input_data.split('.')[-1]
        except :
            input_format ='smi'
            
    # If the file is packmol recognizable, there is no need to convert
    if input_format in PACKMOL_FORMATS :
        with open(input_data) as in_file :
            output_file = tempfile.NamedTemporaryFile(mode='w+',suffix=input_format)
            output_file.write(in_file.read())
    else :
        # The file is packmol compatible format.  Conversion to xyz is needed. 
        # If do not have openbabel, then we are pretty screwed here.
        if not _has_openbabel :
            raise ImportWarning["Can't import OpenBabel for format conversion. " +
                       "Please use xyz format instead."]
        # Create an xyz temp file and perform conversion
        output_file = tempfile.NamedTemporaryFile(mode='w+',suffix='xyz')
        if input_format=='smi' or imput_format=='smiles' :
            mol=pb.readstring(input_format,input_data)
            mol.make3D()
        else :
            mol=pb.readfile(input_format,input_data)
        output_file.write(mol.write('xyz'))
    
    output_file.flush()
    output_file.seek(0)
    return output_file



class Packmol (object):
    """Wrapper class for the packmol molecule packer program
    """
        
    def __init__(self, region_type='sphere',  dimensions=10.0 , 
                 tolerance=3.0 ,  **kwargs ):
        """Construct packmol run by providing simulation parameters.  
       
        """
        self._dimensions  = dimensions
        self._tolerance   = tolerance
        self._region_type = region_type
        self._packmol_arguments = kwargs
        self._structure_list = []

    def add_structure(self, input_data, input_format='auto', count=1) :
        """Add one or more structure group to be packed.
        
        Args:
            input_data: molecular structure of the new group(s).  A 
                string that contains filename that contains this 
                structure, or a SMILES string.  
        """
        self._structure_list.append({ "structure": convert_to_xyz(input_data) ,
                         "number" : count })

    def pack(self) :
        pass
        #Not implemented yet

    def autosize(self) :
        pass
        #Not implemented yet
        
# When called as standalone, perform a doctest run
if __name__ == __main__ :
    import doctest
    doctest.testmod()