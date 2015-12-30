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
import subprocess
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


def convert_to_xyz(input_data,input_format='auto',force_field='uff') :
    """Convert input format unrecognized by packmol into xyz files.
    
    If the input molecular structure data is in a xyz file, the function 
    simply transfer the data to a temp file. Otherwise, the function use 
    OpenBabel to convert into xyz format, then save this to a temp file.
    
    Args:
        input_data: Filename or SMILES string to be converted
        input_format: Format of the input.  When the value is 'auto' which is 
            the default, the function will try to determine the type of input 
            file. In such event, if the input string does not correspond to an 
            existing file that can be opened, the string is instead taken to 
            be a SMILES string. 
    
    Returns:
        A named temporary file containing the molecular structure in xyz format.  
        Note that since this is a temp file, it will automatically be deleted 
        when the file object is unreferenced or destroyed.
    """
    input_format=input_format.lower()
    if input_format == 'auto' :
        try :
            with open(input_data) as input_file :
                input_format = input_data.split('.')[-1]
        except :
            input_format ='smi'
            
    # Create an xyz temp file and perform conversion
    output_file = tempfile.NamedTemporaryFile(mode='w+',suffix='xyz')            
    # If the file is xyz, there is no need to convert
    if input_format == 'xyz' :
        with open(input_data) as in_file :
            output_file.write(in_file.read())
    else :
        # Conversion to xyz is needed. 
        # If do not have openbabel, then we are pretty screwed here.
        if not _has_openbabel :
            raise ImportWarning["Can't import OpenBabel for format conversion. " +
                       "Please use xyz format instead."]
        # Construct a openbabel `pybel.Molecule` object then export to xyz
        if input_format=='smi' or imput_format=='smiles' :
            mol=pb.readstring(input_format,input_data)
            mol.make3D(forcefield=force_field,steps=200)
        else :
            mol=pb.readfile(input_format,input_data)
        output_file.write(mol.write('xyz'))
    
    output_file.flush()
    output_file.seek(0)
    return output_file



class Packmol (object):
    """Wrapper class for the packmol molecule packer program
    """
        
    def __init__(self, command_line='packmol', region_type='sphere', 
                 dimensions=10.0, tolerance=3.0 , output='packed.xyz', 
                 nloop=200, opt_force_field='uff',  **kwargs ):
        """Construct packmol run by providing simulation parameters.  
       
        """
        self._dimensions  = dimensions
        self._tolerance   = tolerance
        self._nloop       = nloop
        self._region_type = region_type.lower()
        self._opt_force_field   = opt_force_field
        self._packmol_arguments = kwargs
        self._command_line = command_line
        self._structure_list = []

    def _prepare_input_file (self) :
        """Create temp file to be used as packmol input file"""
        self._packmol_input = tempfile.NamedTemporaryFile(mode="w+",suffix="inp")
        self._packmol_input.write('tolerance {}\n'.format(self._tolerance))
        self._packmol_input.write('nloop {}\n'.format(self._nloop))
        self._packmol_input.write('filetype xyz\n')
        self._packmol_input.write('output {}\n'.format(self._output))
        # Write down other arguments 
        for key,value in self._packmol_arguments.items():
            self._packmol_input.write('{0} {1}\n'.format(key,value))
        # 
        for mol_type in self._structure_list :
            self._packmol_input.write('structure {}\n'.format(mol_type['structure'].name))
            self._packmol_input.write('  number {}\n'.format(mol_type['number']))
            if self._region_type == 'sphere' :
                self._packmol_input.write(
                    '  inside sphere 0.0 0.0 0.0 {}\n'.format(self._dimensions))
            elif self._region_type :
                raise NotImplementedError('Currently only sphere regions are supported.')
            self._packmol_input.write('  number {}\n'.format(mol_type['number']))
            self._packmol_input.write('end structure\n')
        self._packmol_input.flush()
        self._packmol_input.seek(0)
        return self._packmol_input
    
    def add_structure(self, input_data, input_format='auto', count=1) :
        """Add one or more structure group to be packed.
        
        Args:
            input_data: molecular structure of the new group(s).  A 
                string that contains filename that contains this 
                structure, or a SMILES string.  
        """
        self._structure_list.append({ 
                "structure": convert_to_xyz(input_data,input_format,self._opt_force_field) ,
                "number"   : count })

    def pack(self, print_output=False) :
        outfile = tempfile.NamedTemporaryFile(mode='w+')
        
        code = subprocess.call(self._command_line, stdin=_prepare_input_file(),
                        stdout=outfile )
        outfile.seek(0)
        packmol_output=outfile.read()
        if print_output :
            sys.stdout.write(packmol_output)
        if code or 'ERROR' in packmol_output or 'STOP' in packmol_output :
            raise Exception('Packing failed. exit code={}'.format(code))

    def autosize(self) :
        pass
        #Not implemented yet
        
# When called as standalone, perform a doctest run
if __name__ == __main__ :
    import doctest
    doctest.testmod()