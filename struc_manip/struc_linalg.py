"""
version 1.0.1
author: Felix Plasser, University of Vienna, Institute for Theoretical Chemistry
Waehringerstr. 17, 1090, Vienna, Austria
usage: wrapper to python-openbabel for molecular structure format conversion
    and a package to perform linear algebra operations on structures.
    
    In addition to the formats supported by openbabel this package provides an interface
    to Columbus "col" and a Tinker interface "txyz2".
"""

# improvment: use return to pass structures by value (?).

import os, shutil, locale
import numpy
try:
    import openbabel
except ImportError:
    print " *** Warning: python-openbabel not found! ***"
    print " Using emulation program with limited capabilities ..."
    import OB_repl as openbabel
import superposition
import file_handler, units

Z_symbol_dict = {1:'H',6:'C',7:'N',8:'O',15:'P',16:'S',17:'Cl'}
symbol_Z_dict = {}
for key,val in Z_symbol_dict.iteritems():
    symbol_Z_dict[val] = key
Z_symbol_dict[99] = 'X'

class structure:
    """
    Class to manipulate a structure.
    """
    def __init__(self, name=''):
        self.name = name
        self.new_types = ['txyz2','col'] # these are defined here

    def read_file(self, file_path, file_type='tmol'):
        """
        Read in the structure from a file.
        """
        self.file_path = file_path
        self.file_type = file_type
        self.mol = openbabel.OBMol()
        
        if self.file_type in self.new_types:
            self.read_new_type()
        else:
            obconversion = openbabel.OBConversion()
            obconversion.SetInFormat(file_type)
            obconversion.ReadFile(self.mol, file_path)
            
    def read_new_type(self):
        """
        Routine for reading a file type not contained in open babel.
        """
        infile = open(self.file_path,'r')
        line = infile.readline()

        if self.file_type == 'txyz2': # Tinker format preserving all information
           self.tinker_symbs = []
           self.tinker_extra = []
           num_at_chk = int(line.split()[0])
           line = infile.readline()
           while(line!=''):
              words = line.split()
              if len(words) == 0:
                  break
              obatom = openbabel.OBAtom()
              self.tinker_symbs.append(words[1])
              lett1 = words[1][0]
              if words[1] in symbol_Z_dict:
                 at_num = symbol_Z_dict[words[1]]
              elif words[1][0] in symbol_Z_dict:
                 at_num = symbol_Z_dict[words[1][0]]
              else:
                 at_num = 99
              obatom.SetAtomicNum(at_num)

              coords = [float(word) for word in words[2:5]] 
              obatom.SetVector(*coords)
              self.tinker_extra.append(words[5:])
              # maybe one could use SetSymbol here

              self.mol.AddAtom(obatom)
              line = infile.readline()

           #print self.mol.NumAtoms()
           assert(num_at_chk==self.mol.NumAtoms())
        elif self.file_type == 'col':
           while(line!=''):
              words = line.split()
              obatom = openbabel.OBAtom()
              obatom.SetAtomicNum(int(float(words[1])))
              coords = [float(word)*units.length['A'] for word in words[2:5]] 
              obatom.SetVector(*coords)

              self.mol.AddAtom(obatom)
              line = infile.readline()
        else:
           print 'type %s not supported'%self.file_type
           exit(1)

        infile.close()

    def get_mol(self, mol, file_path, file_type='xyz'):
        """
        Read in an openbabel mol that is passed from a different routine.
        Can be used for accessing multiple structure xyz files.
        """
        self.file_path = file_path
        self.file_type = file_type
        self.mol = mol

    def read_file_vector(self, def_file_path, file_type, vector):
        """
        Initialise the structure with a default .mol from a file and a vector.
        """
        self.read_file(def_file_path, file_type) # like this because copying objects doesn't work

        for i in xrange(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            atom.SetVector(vector[3*i], vector[3*i+1], vector[3*i+2])
        #print 'read_file_vector done'

    def read_file_3xN_matrix(self, def_file_path, file_type, coor_mat):
        """
        Initialise the structure with a default .mol and a 3xN matrix.
        """
        self.read_file(def_file_path, file_type)

        #print self.mol.NumAtoms(), file_type

        #for i, atom in enumerate(openbabel.OBMolAtomIter(self.mol)): ## this only works with the new version

        self.read_3xN_matrix(coor_mat)

    def read_3xN_matrix(self, coor_mat, at_list=None):
        """
        Read in a 3xN matrix.
        If an optional <at_list> is specified, it is assumed that the coordinates
          specified correspond to these atoms and only these atoms are changed.
        """
        if at_list == None:
            at_list = [i+1 for i in range(self.mol.NumAtoms())]

        for imat,iat in enumerate(at_list):
            atom = self.mol.GetAtom(iat)
            atom.SetVector(coor_mat[imat][0], coor_mat[imat][1], coor_mat[imat][2])

    def ret_vector(self):
        " All the coordinates in one vector "
        vec_list = []
        for i in xrange(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            vec_list += [atom.x(), atom.y(), atom.z()]
           
        return numpy.array(vec_list)

    def ret_3xN_matrix(self, at_list=None):
        """
        Return coordinates in a 3 x N matrix.
        If <at_list> is specified only the atoms with those indices are considered.
        """
        mat_list = []

        if at_list == None:
            at_list = [i+1 for i in range(self.mol.NumAtoms())]

        for i in at_list:
            atom = self.mol.GetAtom(i)
            mat_list += [[atom.x(), atom.y(), atom.z()]]
            
        return numpy.array(mat_list)
        
    def ret_rotated_structure(self, theta, vec, name=''):
        """
        Return the structure rotated by <theta> around <vec>.
        """
        if name == '': name = self.name
        
        coor_mat = self.ret_3xN_matrix()
        
        sup = superposition.superposition()
        sup.set_rotation(theta, vec)
        
        #print coor_mat
        #print sup.ret_rotation_matrix()
        coor_mat = numpy.dot(coor_mat, sup.ret_rotation_matrix().transpose())
        
        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, coor_mat)
        return ret_struc
        
    def ret_moved_structure(self, add_vec, name=''):
        """
        Move the structure by <add_vec>.
        """
        if name == '': name = self.name
        
        coor_mat = self.ret_3xN_matrix()
        coor_mat += add_vec
        
        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, coor_mat)
        return ret_struc

    def ret_superimposed_structure(self, struc, mass_wt_pw=1, name='', manual_wts=[]):
        """
        Return this structure superimposed on another structure.
        <self> is fitted onto <struc>.
        <mass_wt_pw>=1 means regular mass weighting.
        <manual_wts> is a nested list [[ind1,wt1],...] with weights that can be put in by hand. Indices start with 1.
        """
        if name == '': name = self.name
        
        coor_mat = self.ret_3xN_matrix()
        
        mass_vect = self.ret_mass_vector(power=mass_wt_pw)
        for ind, weight in manual_wts:
            mass_vect[ind-1] = weight
        
        sup = superposition.superposition()
        sup.superimpose(ref_points=struc.ret_3xN_matrix(), mv_points=coor_mat, weights=mass_vect)
        sup.print_all_info()

        coor_mat = coor_mat - sup.ret_mv_av()
        coor_mat = numpy.dot(coor_mat, sup.ret_rotation_matrix().transpose())
        coor_mat = coor_mat + sup.ret_ref_av()

        #print coor_mat

        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, coor_mat)
        
#         print 'rmsd', sup.ret_rmsd()
#        print 'angle', sup.ret_rotation_angle()
#        print 'axis', sup.ret_rotation_axis()
#         print 'centers', sup.ret_ref_av(), sup.ret_mv_av()
#        print 'matrix'
#        print sup.ret_rotation_matrix()
#         print 'orthogonal matrix?'
#         print numpy.dot(sup.ret_rotation_matrix(), sup.ret_rotation_matrix().transpose())

        return ret_struc
        
    def ret_renumbered_structure_file(self, ren_file, name=''):
        """
        Renumber the structure with input from <ren_file>.
        """
        ren_lines = open(ren_file, 'r').readlines()
        ren_list = [[eval(word) for word in file_handler.line_to_words(line)] for line in ren_lines]
        return self.ret_renumbered_structure(ren_list, name)

    def ret_renumbered_structure(self, perm_list, name=''):
        """
        Do permutations of atoms according to cyclic permutations in <perm_list>.
        This is important for molecules with symmetry or sigma bond rotations. A superposition can be done afterwards.
        example <perm_list>: [[1,2],[3],[6,7,8]...]
            switches 1 and 2, leaves 3 unchanged, and cyclically exchanges 6,7,8
        All atoms must be included exactly once.
        """
        
        if name == '': name = self.name
        
        start_mat = self.ret_3xN_matrix()

        num_at = self.mol.NumAtoms()
        
        trans_mat = numpy.zeros((num_at, num_at), float) # matrix multiplication with this matrix leads to the structure with permutated atom numbering
        for cycle in perm_list:
            for i in xrange(len(cycle)-1):
                trans_mat[cycle[i]-1, cycle[i+1]-1] = 1
            trans_mat[cycle[len(cycle)-1]-1, cycle[0]-1] = 1

        new_vec = numpy.dot(trans_mat, start_mat)

        ret_struc = structure(name=name)
        ret_struc.read_file_3xN_matrix(self.file_path, self.file_type, new_vec)

        return ret_struc

    def ret_bond_length(self, i, j):
        """
        Return the distance between atoms indexed i and j.
        """
        
        OBAtom_i = self.mol.GetAtom(i)
        OBAtom_j = self.mol.GetAtom(j)
        
        pos_i = numpy.array([OBAtom_i.x(), OBAtom_i.y(), OBAtom_i.z()])
        pos_j = numpy.array([OBAtom_j.x(), OBAtom_j.y(), OBAtom_j.z()])

        return numpy.dot(pos_i - pos_j, pos_i - pos_j)**.5
    
    def ret_bend(self, i, j, k):
        """
        Return the bending angle between atoms indexed i, j, k.
        """
        
        OBAtom_i = self.mol.GetAtom(i)
        OBAtom_j = self.mol.GetAtom(j)
        OBAtom_k = self.mol.GetAtom(k)
        
        pos_i = numpy.array([OBAtom_i.x(), OBAtom_i.y(), OBAtom_i.z()])
        pos_j = numpy.array([OBAtom_j.x(), OBAtom_j.y(), OBAtom_j.z()])
        pos_k = numpy.array([OBAtom_k.x(), OBAtom_k.y(), OBAtom_k.z()])

        vec1 = pos_i - pos_j
        vec2 = pos_k - pos_j
        
        len_1 = numpy.sqrt(numpy.dot(vec1,vec1))
        len_2 = numpy.sqrt(numpy.dot(vec2,vec2))

        return numpy.arccos(numpy.dot(vec1, vec2) / (len_1*len_2)) / numpy.pi * 180
    
    def ret_tors(self, i, j, k, l):
        """
        Return the torsion angle between atoms indexed i, j, k, l.
        """
        # Dihedral angle computed according to (http://en.wikipedia.org/wiki/Dihedral_angle) to get the full 360 deg range.
        
        OBAtom_i = self.mol.GetAtom(i)
        OBAtom_j = self.mol.GetAtom(j)
        OBAtom_k = self.mol.GetAtom(k)
        OBAtom_l = self.mol.GetAtom(l)
        
        pos_i = numpy.array([OBAtom_i.x(), OBAtom_i.y(), OBAtom_i.z()])
        pos_j = numpy.array([OBAtom_j.x(), OBAtom_j.y(), OBAtom_j.z()])
        pos_k = numpy.array([OBAtom_k.x(), OBAtom_k.y(), OBAtom_k.z()])
        pos_l = numpy.array([OBAtom_l.x(), OBAtom_l.y(), OBAtom_l.z()])

        vec1 = pos_j - pos_i
        vec2 = pos_k - pos_j
        vec3 = pos_l - pos_k
        
        cross1 = numpy.cross(vec1, vec2)
        cross2 = numpy.cross(vec2, vec3)
        
        norm2 = numpy.sqrt(numpy.dot(vec2, vec2))
        dot1 = numpy.dot(vec1, cross2)
        dot2 = numpy.dot(cross1, cross2)
        
        return numpy.arctan2(norm2 * dot1, dot2) / numpy.pi * 180
    
    def ret_symbol(self, i):
        """
        Returns the symbol of atom i.
        """
        return Z_symbol_dict[self.mol.GetAtom(i).GetAtomicNum()]

    def ret_num_at(self):
        """
        Returns the number of atoms in the file.
        """
        return self.mol.NumAtoms()
    
    def ret_mass_vector(self, power):
        """
        Returns a vector with the masses of the atoms (each 1 time) taken to the <power> power.
        """
        mass_list = []
        for i in xrange(self.mol.NumAtoms()):
            atom = self.mol.GetAtom(i+1)
            mass_list += [atom.GetAtomicMass()**power]

        return numpy.array(mass_list, float)
        
    def ret_partition(self,cutBonds=[]):
        """
        Return a partition according to different non-bonded molecules.
        cutBonds is a list of tuples for bonds to cut. The order does not matter
        e.g. cutBonds=[(3,4),(7,8)]
        """
        at_lists = []
        chk_list = []
        remaining_atoms = [self.mol.NumAtoms()-i for i in xrange(self.mol.NumAtoms())]
        
        while(len(remaining_atoms)>0): # do the loop as long as there are still atoms which are not in at_lists
            if len(chk_list) > 0:
                curr_at = chk_list.pop()
            else:
                #print 'new at_list'
                at_lists.append([])
                curr_at = remaining_atoms.pop()
                at_lists[-1].append(curr_at)
            #print 'curr_at:',curr_at
            
            atom = self.mol.GetAtom(curr_at)
            for bonded in openbabel.OBAtomAtomIter(atom):
                bind = bonded.GetIdx()
                if bind in remaining_atoms:
                    if  (curr_at,bind) in cutBonds or (bind,curr_at) in cutBonds:
                        print 'cutting bond %i-%i'%(curr_at,bind)
                    else:
                        del remaining_atoms[remaining_atoms.index(bind)]
                        chk_list.append(bind)
                        at_lists[-1].append(bind)
                        #print bind
            
            
        
        return at_lists
    
    def make_coord_file(self, file_path, file_type='tmol'):
        """
        Write the structure file.
        """
        if file_type in self.new_types:
            self.make_coord_new(file_path, file_type)
        else:
            obconversion = openbabel.OBConversion()
            obconversion.SetOutFormat(file_type)
            obconversion.WriteFile(self.mol, file_path)

    def make_coord_new(self, file_path, file_type):
        """
        Routine for reading a file type not contained in open babel.
        """
        outfile = open(file_path,'w')
        num_at = self.mol.NumAtoms()

        if file_type == 'txyz2':
          outfile.write('%i from MSMT\n'%num_at) 
          for ind in xrange(1, num_at+1):
            obatom = self.mol.GetAtom(ind)
            outstr  = ' %6i'%ind
            outstr += ' %4s'%self.tinker_symbs[ind-1]
            outstr += ' % 10.6f'%obatom.x()
            outstr += ' % 10.6f'%obatom.y()
            outstr += ' % 10.6f'%obatom.z()
            for extr in self.tinker_extra[ind-1]:
                outstr += ' %6s'%extr
            outfile.write(outstr+'\n')
        elif file_type == 'col':
          for ind in xrange(1, num_at+1):
            obatom  = self.mol.GetAtom(ind)
            outstr  = '%2s'%Z_symbol_dict[obatom.GetAtomicNum()]
            outstr += ' %7.1f'%obatom.GetAtomicNum()
            outstr += '% 14.8f'%(obatom.x() / units.length['A'])
            outstr += '% 14.8f'%(obatom.y() / units.length['A'])
            outstr += '% 14.8f'%(obatom.z() / units.length['A'])
            #outstr += '   %12.8f'%obatom.GetAtomicMass()
            # take the mass of the most abundant isotope rather than the average
            outstr += '% 14.8f'%obatom.GetExactMass()
            outfile.write(outstr+'\n')

        outfile.close()

class mol_calc:
    """
    Class for doing calculations.
    A default file path is needed.
    For larger molecules it is not really efficient to create the whole structure instance everytime.
    It would be more useful to do it only with the coordinates.
    """
    def __init__(self, def_file_path, file_type='tmol'):
        # def_file_path is separately read in because with passing arguments the original objects would also be changed
            # maybe this would not occur if "return" would be used
        self.def_file_path = def_file_path
        self.file_type = file_type

    # return a set of structures
    def gram_schmidt(self, structures):
        """
        Returns an orthogonal system that spans the same space as *structures*.
        """
        b_list = [] # list of the orthogonal vectors

        for struc in structures:
            cn1 = struc.ret_vector()

            bn1 = cn1
            for bk in b_list:
                bn1 -= numpy.dot(cn1, bk) / numpy.dot(bk,bk) * bk

            bn1 = bn1 / numpy.dot(bn1,bn1)**.5
            b_list += [bn1]

        #print numpy.array([[numpy.dot(bi, bj) for bi in b_list] for bj in b_list])

        return [self.make_structure(bi) for bi in b_list]

        
    # return a structure
    def make_structure(self, vector, name=''):
        """
        Return the structure that corresponds to a vector.
        """
        ret_struc = structure(name)        
        ret_struc.read_file_vector(self.def_file_path, self.file_type, vector)
        return ret_struc
    
    def add(self, struc1, struc2):        
        return self.make_structure(struc1.ret_vector() + struc2.ret_vector())
        
    def subtract(self, struc1, struc2):
        #print struc1.ret_vector()
        #print struc2.ret_vector()
        return self.make_structure(struc1.ret_vector() - struc2.ret_vector())
        #print "subtract done"

    def scalar_mult(self, scalar, struc):    
        return self.make_structure(scalar * struc.ret_vector())

    def mean_structure(self, struc_list):
        """
        Return the mean structure of <struc_list>.
        """
        factor = 1. / len(struc_list)
        sum_struc = reduce(lambda struc1, struc2: self.add(struc1, struc2), struc_list)

        return self.scalar_mult(factor, sum_struc)

    def projected_structure(self, struc, ref_struc, diff_strucs, mass_wt=False):
        " Project the structure <struc> into the (hyper)plane defined by structures <ref_struc> and <diff_strucs> "
        shifted_struc = self.subtract(struc, ref_struc)

        o_diff_strucs = self.gram_schmidt(diff_strucs)

        ret_struc = ref_struc
        for o_struc in o_diff_strucs:
            factor = self.inner_product(o_struc, shifted_struc)
            plus_struc = self.scalar_mult(factor, o_struc)
            ret_struc = self.add(ret_struc, plus_struc)

        return ret_struc

    # return a scalar        
    def inner_product(self, struc1, struc2, mass_wt_pw=1):
        if mass_wt_pw == 0:
            return numpy.dot(struc1.ret_vector(), struc2.ret_vector())
        else:
            ## divide by the trace ???
            vec1 = numpy.dot(struc1.ret_vector(), self.ret_mass_matrix(power=mass_wt_pw))
            return numpy.dot(vec1, struc2.ret_vector())

    def norm(self, struc, mass_wt_pw=1):
        return numpy.sqrt(self.inner_product(struc, struc, mass_wt_pw))
        
    def distance(self, struc1, struc2, mass_wt_pw=1):
        """
        Return the distance between to structures in amu**(-1/2)*A. The distance as defined here is the RMSD times the squareroot of the molecular mass.
        """
        return self.norm(self.subtract(struc1, struc2), mass_wt_pw)
        #print 'distance done'

    def angle(self, struc1, struc2, mass_wt_pw=False):
        """
        Return the angle between two structure vectors or normal modes in degrees.
        """
        return numpy.arccos(self.inner_product(struc1, struc2, mass_wt_pw=mass_wt_pw) / self.norm(struc1, mass_wt_pw=mass_wt_pw) / self.norm(struc2, mass_wt_pw=mass_wt_pw)) / numpy.pi * 180

    # return different output
    def ret_mass_vector(self, power=1):
        """
        Returns a vector with the masses of the atoms (each 1 time) taken to the <power> power.
        """
        def_struc = structure()        
        def_struc.read_file(self.def_file_path, self.file_type) # mass weighing in this case not needed
        return def_struc.ret_mass_vector(power=power)
        
    def ret_mass_matrix(self, power=1):
        """
        Returns a diagonal matrix that contains the masses of the atoms (each 3 times).
        <power=.5> gives the squareroots which is typical mass weighting.
        """
        def_struc = structure()        
        def_struc.read_file(self.def_file_path, self.file_type) # mass weighing in this case not needed
        mass_list = []
        #for atom in openbabel.OBMolAtomIter(def_struc.mol):
        for i in xrange(def_struc.mol.NumAtoms()):
            atom = def_struc.mol.GetAtom(i+1)
            mass_list += 3*[atom.GetAtomicMass()**power]

        ret_mat = numpy.zeros((len(mass_list), len(mass_list)), dtype=float)
        for i in xrange(len(mass_list)):
            ret_mat[i,i]=mass_list[i]

        return ret_mat
    
    def distance_table(self, struc_list, mass_wt_pw, digits=4):
        """
        Return a table with the RMSDs between the structures in struc_list.
        <digits> specifies how many digits after the decimal point are printed out.
        """
        #print 'in distance_table'
        tm = file_handler.table_maker([5] + (len(struc_list)-1) * [digits+4])
        tm.write_line([''] + [struc.name for struc in struc_list[1:]])
        for i,struc1 in enumerate(struc_list[:-1]):
            tm.write_line([struc1.name] + [''] * i + [locale.format("%.*f", (digits, self.distance(struc1, struc2, mass_wt_pw=mass_wt_pw))) for struc2 in struc_list[i+1:]])
        return tm.return_table()
                 
if __name__ == '__main__':
    test_struc = structure()
    test_struc.read_file(file_path = '/home2/plasserf/BIP/opt/gs/b3lyp/SVP/coord', file_type='tmol')
    print test_struc.ret_3xN_matrix()
