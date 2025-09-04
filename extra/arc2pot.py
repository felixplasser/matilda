#!/usr/bin/python

print """
==================================================================
arc2pot.py

version: 1.0
author: Felix Plasser (felix.plasser@unive.ac.at)
 University of Vienna, Institute for Theoretical Chemistry
 Waehringerstr. 17, 1090, Vienna, Austria

usage: convert a tinker archive into potential.x file for Columbus

syntax: "arc2pot.py <tinker-arc> <prm-file>"
input file: arc2pot.in

==================================================================
"""

import os, sys

au2ang = 0.529177249

class converter:
   def __init__(self,arc_file,prm_file,at_list,frames,mirror=None,scale=1.,sep=False,new_charge={},link_atom = {},zero_scatter=[]):
      """
      <sep> specifies that data is read from separate data files rather than an archive.
           <arc_file> refers then to the stem name
      """
      self.arc_file = arc_file
      self.prm_file = prm_file
      self.at_list = at_list
      self.frames = frames
      self.mirror = mirror
      self.scale = scale
      self.sep = sep
      self.new_charge = new_charge
      self.link_atom = link_atom # structure: {QM_ind:[MM_ind, ratio, symbol]}
      self.zero_scatter = zero_scatter
      
      if mirror == None:
         self.num_frames = len(self.frames)
      else:
          self.num_frames = 2*len(self.frames)
          
      self.link_atom_aux = {} # structure: {MM_ind:coors}
      for lind,lat in self.link_atom.iteritems():
          MMind = lat[0]
          self.link_atom_aux[MMind] = []
          #if MMind in self.new_charge:
          #    print 'Keeping user defined charge for MM link atom'
          #    print '  new_charge[%i]=%.5f'%(MMind,self.new_charge[MMind])
          #else:
          #    print 'Charge of MM link atom set to zero'
          #    print '  new_charge[%i]=0.0'%(MMind)
          #    new_charge[MMind] = 0.0

      
   def get_charges(self):
      """
      Get the charges from the parameter file.
      """
      if self.sep:
          ofil = '%s.%03i'%(self.arc_file,1)
      else:
          ofil = self.arc_file
      try:
         arc = open(ofil,'r')
      except:
         print " Could not open %s"%ofil
         print " Please use 'sep = True' if you have separate files rather than one archive!"
         print
         raise
      
      words = arc.readline().split()
      self.num_at = int(words[0])
      print 'Number of atoms: ', self.num_at
      
      # first find out what kinds of atoms there are
      self.ch_dict = {}
      for i in xrange(self.num_at):
         words = arc.readline().split()
         at_typ = words[5]
         if not at_typ in self.ch_dict:
            self.ch_dict[at_typ] = None
            
      arc.close()            
      
      # now read the charges for the atoms
      prm = open(self.prm_file,'r')
      for line in prm:
         if 'charge' in line:
            words = line.split()
            at_typ = words[1]
            if at_typ in self.ch_dict:
               charg = float(words[2])
               self.ch_dict[at_typ] = charg
               
      prm.close()
      print "Charges read in:"
      for at_typ,charg in self.ch_dict.iteritems():
         print " %s -> %9.5f"%(at_typ,charg)
      print
      
   def prep_run(self):
       """
       This is a preparation pass through the file.
       """
       print "Computing charges to scatter ..."
       
       all_zero = []
       all_scatter = []
       tar_charges = []
       orig_scatter = {}
       
       for zero,scatter in self.zero_scatter:
           tar_charges.append(0.)
           all_zero+= zero
           all_scatter+= scatter
        # charge of the MM system including definitions from the .prm file and "new_charge" definitions
       
       if self.sep:
          ofil = '%s.%03i'%(self.arc_file,1)
       else:
          ofil = self.arc_file
       arc = open(ofil,'r')
       
       arc.readline()
       for j in xrange(self.num_at):
          line = arc.readline()
          words = line.split()
          at_ind = int(words[0])
          at_typ = words[5]
          
          if (at_ind in all_zero) or (at_ind in all_scatter):
             if at_ind in self.new_charge:
                 charg_full = self.new_charge[at_ind]
                 print "Taking user defined charge of % .4f for atom %i"%(charg_full,at_ind)
             else:
                try:
                    charg_full = self.ch_dict[at_typ]
                except KeyError:
                    print "No charge found for:"
                    print line
                    print "Assuming 0"
                    charg_full = 0.
                    
             if at_ind in all_zero:
                for zsind,zs in enumerate(self.zero_scatter):
                    zero,scatter = zs
                    if at_ind in zero:
                        tar_charges[zsind]+= charg_full
                        if lvprt>=1: print "  at_ind: %i, charge: % .5f, zsind: %i, new tar charge: % .5f"%(at_ind,charg_full,zsind,tar_charges[zsind])
             if at_ind in all_scatter:
                orig_scatter[at_ind] = charg_full
             
       arc.close()
       
       for zsind,zs in enumerate(self.zero_scatter):
           zero,scatter = zs
           if len(scatter) > 0:
              print " Distributing charge of % .5f over %i point charges"%(tar_charges[zsind],len(scatter))
              dist_charg = tar_charges[zsind] / float(len(scatter))
              for scind in scatter:
                  new_charge[scind] = orig_scatter[scind] + dist_charg
                  if lvprt>=1: print "  Index: %i, original: % .5f, new: % .5f"%(scind,orig_scatter[scind],new_charge[scind])
           else:
              print " scatter not specified, discarding charge of % .5f"%tar_charges[zsind]
           for zind in zero:
               new_charge[zind] = 0.
      
   def write_pot(self,out_dir=None):
      """
      Read the Tinker archive and write the potential.x file as input for Columbus.
      This is performed at the same time to save memory.
      If <out_dir> is given, separate tinker files are taken from <out_dir>.
      """
      tot_charg = 0.
      if not self.sep: arc = open(self.arc_file,'r')
      pot = open('potential.xyz','w')
      if not out_dir == None:
          try:
            os.mkdir(out_dir)
          except OSError:
            print 'out_dir %s not created'%out_dir
      
      num_pot = self.num_frames * len(self.at_list)
#      if not self.mirror==None:
#         num_pot = num_pot * 2
      pot.write(str(num_pot)+'\n')
      print "Writing out %i point charges"%(num_pot)
      
      last_frame = 0
      for frame in self.frames:
         step_str = ''
         ENV_str = '%s\n'%len(self.at_list)
         
         print "Analyzing frame", frame
         frame_skip = frame - last_frame - 1
         if frame_skip < 0:
            print 'Please enter frames in ascending order!'
            sys.exit(1)
            
         if self.sep:
            ofil = '%s.%03i'%(self.arc_file,frame)
            print ' Opening file %s ...'%ofil
            arc = open(ofil,'r')
         else:
            # skip the unused frames
            for i in xrange(frame_skip):
                for j in xrange(self.num_at + 1):
                    line = arc.readline()
                  
         coor=[0.,0.,0.]
         line = arc.readline() # header with the number of atoms
         step_str += line
         for j in xrange(self.num_at):
            line = arc.readline()
            step_str += line
            
            words = line.split()
            at_ind = int(words[0])
            if at_ind in self.at_list:
               ENV_str += line
               name = words[1]
               coor = [float(word)/au2ang for word in words[2:5]]
               
               at_typ = words[5]
               if at_ind in self.new_charge:
                    charg_full = self.new_charge[at_ind]
               else:
                    try:
                        charg_full = self.ch_dict[at_typ]
                    except:
                        print "No charge found for:"
                        print line
                        raise
                   
               charg = charg_full / float(self.num_frames) * self.scale
               tot_charg+= charg
               outstr = '%5s %14.8f %14.8f %14.8f %14.8f\n'%(name,charg,coor[0],coor[1],coor[2])
#               for icoor in coor:
 #                  if icoor == 0.:
  #                     print 'Warning coor=0. for'
   #                    print outstr
               pot.write(outstr)
               # perform additional permuations to symmetrize
               if not self.mirror==None:
                  coor[mirror] = -coor[mirror]
                  outstr = '%5s %14.8f %14.8f %14.8f %14.8f\n'%(name,charg,coor[0],coor[1],coor[2])
                  pot.write(outstr)
                  
            if at_ind in self.link_atom_aux:
               coor = [float(word) for word in words[2:5]]
               self.link_atom_aux[at_ind] = coor
               
         last_frame = frame
         if not out_dir == None:
             step_file = open(out_dir+'/tinker_step.%5.5i'%frame,'w')
             step_file.write(step_str)
             step_file.close
             
             step_file = open(out_dir+'/tinker_ENV.%5.5i'%frame,'w')
             step_file.write(ENV_str)
             step_file.close
      
      
      arc.close()
      pot.close()
      
      print "potential.xyz written"
      print "Total charge: % .6f"%tot_charg
      
   def write_geom(self):
      """
      Write out the geometry that is not in <at_list> for the first frame. Typically this could be the QM region.
      """
      if self.sep:
          ofil = '%s.%03i'%(self.arc_file,1)
      else:
          ofil = self.arc_file
      arc = open(ofil,'r')
      geo = open('geom.xyz','w')
      coor=[0.,0.,0.]
      
      num_geo = self.num_at - len(self.at_list) + len(link_atom)
      geo.write('%i\n\n'%num_geo)
      
      arc.readline() # header with the number of atoms
      for j in xrange(self.num_at):
        line = arc.readline()
        words = line.split()
        at_ind = int(words[0])
        if at_ind not in self.at_list:
            name = words[1][0] # only take the first letter. this may have to be manually corrected in some cases
            coor = [float(word) for word in words[2:5]]
            outstr = '%5s %14.8f %14.8f %14.8f\n'%(name,coor[0],coor[1],coor[2])
            geo.write(outstr)
        
        if at_ind in link_atom:
            [mm_ind,ratio,symb] = link_atom[at_ind]
            coor_qm = [float(word) for word in words[2:5]] 
            
            # the position of the MM atom in the last frame is taken here. Another possibility would be to take an average.
            coor_mm = self.link_atom_aux[mm_ind]
            
            coor_new = [coor_qm[i] + ratio*(coor_mm[i]-coor_qm[i]) for i in xrange(3)]
            outstr = '%5s %14.8f %14.8f %14.8f\n'%(symb,coor_new[0],coor_new[1],coor_new[2])
            geo.write(outstr)
        
      arc.close()
      geo.close()

if __name__=='__main__':
   from file_handler import to
   if len(sys.argv) < 1+2:
      print 'At least two arguments required!'
      sys.exit()

   arc_file = sys.argv[1]
   prm_file = sys.argv[2]
   #at_list = range(1,2989)
   #frames = [2,4,6]
   mirror = None
   scale = 1.
   sep = False
   zero_scatter = []
   new_charge = {}
   link_atom = {}
   lvprt = 0
   
   try:
    execfile('arc2pot.in')
   except IOError:
    print "Error: Please create input file arc2pot.in, example:"
    print "# --- arc2pot.in ---"
    print "at_list = to(1,10) + to(20,100)\nframes = [2,4,6]\nscale = 1."
    #print "new_charge[5] = 0.122\nnew_charge[100] = -0.345"
    exit(1)
   
   conv = converter(arc_file, prm_file,at_list,frames,mirror,scale,sep,new_charge,link_atom,zero_scatter)
        
   conv.get_charges()
   if not zero_scatter == []:
       conv.prep_run()
   
   #conv.write_pot(out_dir='TINKER_steps')
   conv.write_pot()
   conv.write_geom()