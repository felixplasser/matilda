"""
version 1.0.1
author: Felix Plasser, University of Vienna, Institute for Theoretical Chemistry
Waehringerstr. 17, 1090, Vienna, Austria
usage: OB_repl is intended as a replacement that supplies some of the functionalities of the python-openbabel
   package in case this is not available.
"""

import struc_linalg

Z_mass_dict = {1:1.008,6:12.011,7:14.007,8:16.,15:30.97,16:32.066,17:35.453}
Z_exact_mass_dict = {1:1.00782504,6:12.,7:14.00307401,8:15.99491464,17:34.96885273}

class OBConversion:
    """
    Emulate the openbabel OBConversion class with support for a reduced number of file formats.
    """
    def __init__(self):
        self.in_format = None
        self.out_format = None
        
    def SetInFormat(self,format):
        self.in_format = format
        
    def SetOutFormat(self,format):
        self.out_format = format
        
    def ReadFile(self,mol,file):
        """
        Read information from a file into mol.
        """
        if self.in_format == 'xyz':
            self.read_xyz(mol,file)
        elif self.in_format == None:
            print "Input format for %s not set!"%file
            exit(14)
        else:
            print "File format %s not supported for input!"%self.in_format
            print "Install python-openbabel for complete support."
            exit(15)
            
    def read_xyz(self,mol,file):
        """
        Read file in xyz format (Angstrom).
        """
        infile = open(file,'r')
        
        # the first two lines do not contain any coordinate information
        line = infile.readline()
        line = infile.readline()
        line = infile.readline()
        
        
        while(line!=''):
              words = line.split()
              obatom = OBAtom()
              obatom.SetAtomicNum(struc_linalg.symbol_Z_dict[words[0]])
              coords = [float(word) for word in words[1:4]] 
              obatom.SetVector(*coords)

              mol.AddAtom(obatom)
              line = infile.readline()
              
    def WriteFile(self,mol,file):
        """
        Print structure to a file.
        """
        if self.out_format == 'xyz':
            self.write_xyz(mol,file)
        elif self.out_format == None:
            print "Output format for %s not set!"%file
            exit(14)
        else:
            print "File format %s not supported for output!"%self.out_format
            print "Install python-openbabel for complete support."
            exit(15)
            
    def write_xyz(self,mol,file):
        outfile = open(file,'w')
        num_at = mol.NumAtoms()
        
        outfile.write('%i\n\n'%num_at)
        
        for ind in xrange(1, num_at+1):
            obatom  = mol.GetAtom(ind)
            outstr  = '%2s'%struc_linalg.Z_symbol_dict[obatom.GetAtomicNum()]
            outstr += '% 14.8f'%(obatom.x())
            outstr += '% 14.8f'%(obatom.y())
            outstr += '% 14.8f'%(obatom.z())
            outfile.write(outstr+'\n')

        outfile.close()
        
class OBMol:
    def __init__(self):
        self.atoms = []
        
    def AddAtom(self,Atom):
        self.atoms.append(Atom)
        
    def GetAtom(self,ind):
        """
        Return an atom, indices start with 1.
        """
        return self.atoms[ind-1]
        
    def NumAtoms(self):
        return len(self.atoms)
    
class OBAtom:
    def __init__(self):
        self.AtomicNum = None
        self.Vector = [None,None,None]
        
    def SetAtomicNum(self,AtomicNum):
        self.AtomicNum = AtomicNum
        
    def GetAtomicNum(self):
        return self.AtomicNum
        
    def SetVector(self,*Vector):
        self.Vector=Vector
        
    def GetVector(self):
        return self.Vector
        
    def x(self):
        return self.Vector[0]
    
    def y(self):
        return self.Vector[1]
    
    def z(self):
        return self.Vector[2]
    
    def GetAtomicMass(self):
        return Z_mass_dict[self.AtomicNum]
    
    def GetExactMass(self):
        return Z_exact_mass_dict[self.AtomicNum]