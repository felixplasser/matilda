#!/usr/bin/python

print """
==================================================================
pot_reduce.py

version: 1.0
author: Felix Plasser (felix.plasser@unive.ac.at)
 University of Vienna, Institute for Theoretical Chemistry
 Waehringerstr. 17, 1090, Vienna, Austria

usage: reduce a potential.xyz file by considering close lying point charges

==================================================================
"""

import math,time
#import os, sys

au2ang = 0.529177249

def reduce_pot():
    """
    Reduce point charges.
    """
    conv = converter(orig=orig,scale=scale,read=True)
    conv.summ_chrg()
     
    conv.reduce_both(C=C,r0=r0)
    conv.summ_chrg()
    #conv.print_pchrg(outfile='pot_no-symm.xyz')
    
    if not group=='c1':
        conv.symmetrize(group=group)
        conv.summ_chrg()
        #conv.print_pchrg(outfile='pot_symm.xyz')
        #conv.reduce_both(C=C/2.,r0=r0)
        conv.reduce_both(C=C,r0=r0)
        conv.summ_chrg()
    
    # if a symmetrization is performed before
    #   the charges are discarded here all the information is used
    #   it should be enough to use the corresponding fraction of frames then
    if not sector==[0,0,0]:
        conv.take_sector(sector=sector,mode=1)
        conv.take_sector(sector=sector,mode=-1)
        conv.summ_chrg()
        print

    
    conv.print_pchrg(formout=formout)

def expand(grad_file):
    """
    Expand a gradient.
    """
    grex = grad_expand(grad_file=grad_file)
    grex.grad_expand()
    grex.print_grad()
    
def plot_pchrg():
    """
    Plot the point charges
    """
    conv = converter(orig=orig,scale=scale,read=True)
    conv.plot_pchrg(fig_name='pchrg_z.png',axis=3,sector=[math.pi/8.,2*math.pi/8.],plot_fact=10.)
    conv.plot_pchrg(fig_name='pchrg_y.png',axis=2,sector=[7*math.pi/16.,9*math.pi/16.],plot_fact=10.,figsize=(8,8))

class converter:
    def __init__(self,orig=[0.,0.,0.],scale=[1.,1.,1.],read=False):
        """
        Origin should be the center of the molecule considered in a.u.
        """
        # aside from an origin it should also be easy to add an ellipticity. but this would affect both r0 and C
        
        self.orig = orig
        self.scale = scale
        print "Scale: ", scale
        if read: self.read_pchrg()
        
    def read_pchrg(self,infile='pot_old.xyz'):
        """
        Read the point charges from potential.xyz and convert them to spherical coordinates.
        """
        print "Reading data ... (%s)"%(time.asctime())
        pchrg = open(infile,'r')
        line = pchrg.readline()
        self.num_charg_in = int(line.split()[0])
        print "%i Point charges read in. (%s)"%(self.num_charg_in,time.asctime())
        
        self.chc_pos = []
        self.chc_neg = []
        # Information is arranged as: [r, x, y, z, charge, symb, at_inds]
        # at_inds gives to the original indices corresponding to a reduced point charge
        (x0,y0,z0) = self.orig
       
        at_ind = 1 
        line = pchrg.readline()
        while(line!=''):
            words = line.split()
            symb = words[0]
            (ch,x,y,z) = [float(word) for word in words[1:5]]
            xeff = (x-x0)/self.scale[0]
            yeff = (y-y0)/self.scale[1]
            zeff = (z-z0)/self.scale[2]
            
            r = ((xeff)**2+(yeff)**2+(zeff)**2)**(.5)
            #print x,x0
            at_inds = [at_ind]
            at_ind += 1
                
            if ch > 0:
                self.chc_pos.append([r, xeff, yeff, zeff, ch, symb, at_inds])
            elif ch < 0:
                self.chc_neg.append([r, xeff, yeff, zeff, ch, symb, at_inds])
            
            line = pchrg.readline()
            
        pchrg.close()
        
        self.chc_pos.sort()
        self.chc_neg.sort()
        #print self.ch_pos
        print "Sorting finished. (%s)\n"%time.asctime()
        
    def reduce_both(self,C=.01,r0=0.):
        self.reduce_pchrg(mode=-1,C=C,r0=r0)
        self.reduce_pchrg(mode= 1,C=C,r0=r0)
        
    def reduce_pchrg(self,mode,C=.01,r0=0.):
        """
        This is the main routine which reduce the point charges. Usually this routine should be called separately for positive and negative charges to recover large scale electrostatic moments.
        C is the fineness parameter
        r0 (a.u.) a cutoff radius
        """
        # important for implementation: that pop() and append() at the end of a list a fast, whereas insdie they are slower.
        
        # in this implementation interactions between unchanged point charges are computed several times
        #ch_list.sort()
        
        red=True
        while(red):
            in_list = []
            if mode == 1:
                for charge in self.chc_pos:
                    in_list.append(charge)
                self.chc_pos = []
            elif mode == -1:
                for charge in self.chc_neg:
                    in_list.append(charge)
                self.chc_neg = []
                
            print "Starting reduction cycle (%+1i) with %8i point charges... (%s)"%(mode,len(in_list),time.asctime())
            red=False

            while(len(in_list)>0):
                #print '----------'
                ch_coll1 = in_list.pop()
                [r1,x1,y1,z1,ch1,symb1,at_inds1] = ch_coll1
                
                for ind,ch_coll2 in enumerate(reversed(in_list)):
                    #print ind,ch_coll2
                    
                    rred = max(0., (r1+ch_coll2[0])/2. - r0)
                    dthres = C * rred**2
                    [r2,x2,y2,z2,ch2,symb2,at_inds2] = ch_coll2
                    
                    # (r-r1)**2 is a lower bound for d2 of this iteration and all the following d2's
                    #   dthres will increase only slower (assuming that C>0)
                    #   in other words: charges inside this sphere will definitely not have to be considered
                    #if symb1 == symb2 == "Cl":
                     #   print r1,r1
                    if (r1-r2)**2 > dthres: break
                    
                    d2 = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
                    
                    if d2<=dthres:
#                        print '+++++++++++++++++++'
                        red = True
                        del in_list[len(in_list)-ind-1]
                         
                        chn = ch1 + ch2
                        symbn = 'RD'
                        at_indsn = at_inds1 + at_inds2
                        # -> this has to be done with weighting!
                        xn = (x1*ch1 + x2*ch2) / (ch1+ch2)
                        yn = (y1*ch1 + y2*ch2) / (ch1+ch2)
                        zn = (z1*ch1 + z2*ch2) / (ch1+ch2)
                        
                        rn = (xn**2 + yn**2 + zn**2)**(.5)
                        
                        #self.ch_pos.append([rn, costhn, phin, chn, symbn])
                        ch_coll1 = [rn, xn, yn, zn, chn, symbn, at_indsn]
                        break
                    
                if mode == 1:
                    self.chc_pos.append(ch_coll1)
                elif mode == -1:
                    self.chc_neg.append(ch_coll1)
                    
            if mode == 1:
                self.chc_pos.sort()
            elif mode == -1:
                self.chc_neg.sort()
            
        if mode == 1:
            num_act = len(self.chc_pos)
        elif mode == -1:
            num_act = len(self.chc_neg)
            
        print 'Reduction (%+1i) finished, %8i point charges remaining. (%s)\n'%(mode,num_act,time.asctime())
        
    def symmetrize(self,group):
        """
        Symmetrize point charges according to symmetry group <group>.
        """
        # an additional reduction after symmetrization will not necessarily keep the symmetry (e.g. if there are atoms on the symmetry
        #   elements but it usually does - quadruples are retained, or pairs at a symmetry plane
        print "Starting symmetrization (%s)... (%s)"%(group,time.asctime())
        idop  =  [[( 1,1), ( 1,2), ( 1,3)]]
        csop  =  [[( 1,1), ( 1,2), (-1,3)]]
        
        c2ops =  idop  + \
                 [[(-1,1), (-1,2), ( 1,3)]]
        c2vops = c2ops + \
                 [[(-1,1), ( 1,2), ( 1,3)],\
                  [( 1,1), (-1,2), ( 1,3)]]
        csops = idop + csop
                 
        d2hops = c2vops + csop + \
                 [[( 1,1), (-1,2), (-1,3)],\
                  [(-1,1), ( 1,2), (-1,3)],\
                  [(-1,1), (-1,2), (-1,3)]]
                  
        if group == 'cs':
            ops = csops
        elif group == 'c2':
            ops = c2ops
        elif group == 'c2v':
            ops = c2vops
        elif group == 'd2h':
            ops = d2hops
        else:
            print 'Symmetry group %s not available!'%group
            exit(1)

        num_op = len(ops)
                   
        old_chcs = [chc for chc in self.chc_pos] + [chc for chc in self.chc_neg]
        self.chc_pos = []
        self.chc_neg = []
        
        for chc in old_chcs:
            new_chrg = chc[4] / float(num_op)
            for op in ops:
                new_chc = [chc[0]] + [sign*chc[ind] for (sign,ind) in op] + [new_chrg, chc[5], chc[6]]
                if new_chrg > 0:
                    self.chc_pos.append(new_chc)
                else:
                    self.chc_neg.append(new_chc)
            
            
        print "Symmetrization finished. (%s)\n"%(time.asctime())
        
    def take_sector(self,sector,mode):
        """
        Take only a symmetry unique sector of the charge distribution.
        sector = [sx,sy,sz]
        0 ... take all values, 1 ... take only positive values, -1 ... take only negative values
        """
        print "Taking only the symmetry unique sector (%i,%i,%i)"%(sector[0],sector[1],sector[2])
        in_list = []
        
        if mode == 1:
            for charge in self.chc_pos:
                    in_list.append(charge)
            self.chc_pos = []
        elif mode == -1:
            for charge in self.chc_neg:
                    in_list.append(charge)
            self.chc_neg = []
        
        while(len(in_list)>0):
            #print '----------'
            ch_coll1 = in_list.pop()
            coors = ch_coll1[1:4]
            
            # atoms on the symmetry elements are always retained
            keep = True
            for i in xrange(3):
                if   sector[i] < 0 and coors[i] > 0:
                    keep = False
                elif sector[i] > 0 and coors[i] < 0:
                    keep = False
            
            if keep:
                if mode == 1:
                    self.chc_pos.append(ch_coll1)
                elif mode == -1:
                    self.chc_neg.append(ch_coll1)

        
    def summ_chrg(self):
        """
        Give the sum of positive and negative charges. This should always stay the same.
        """
        # This is fast enough. Theoretically it could be faster using "reduce"
        print "Summing charges (%s)"%(time.asctime())
        ch_pos = 0.
        ch_neg = 0.
        for chc in self.chc_pos:
            ch_pos += chc[4]
        for chc in self.chc_neg:
            ch_neg += chc[4]
            
        print "Pos: %12.8f, neg: %12.8f, net: %12.8f\n"%(ch_pos,ch_neg,ch_pos+ch_neg)
    
    def print_moments(self):
        pass
    
    def plot_dist(self):
        pass
    
    def plot_pchrg(self,fig_name='pchrg.png',axis=3,sector=[math.pi/8.,2*math.pi/8.],plot_fact=10.,figsize=(8,8)):
        """
        Plot point charges lying in a plane (defined by the normal vector), or next to it with a tolerance of tol
        
        <plot_fact> determines how big the dots are plotted. It is assumed that the charge scales with
            the third power of the radius.
        """
        plot_coors = [1,2,3]
        del(plot_coors[plot_coors.index(axis)])
        
        pylab.figure(figsize=figsize)
        
        rint_dict = {}
        zint_dict = {}
        for charge in self.chc_pos:
            [r,x,y,z,ch,symb,at_inds] = charge
            xint = charge[plot_coors[0]]
            yint = charge[plot_coors[1]]
            zint = charge[axis]
            phi = math.atan2(yint, xint)
            
            rint = None
            if sector[0] <= phi <= sector[1]:
                rint = math.sqrt(xint**2. + yint**2.)
            elif (sector[0] <= phi + math.pi <= sector[1]) or (sector[0] <= phi - math.pi <= sector[1]):
                rint = -math.sqrt(xint**2. + yint**2.)
            if not rint == None:
                if ch in rint_dict:
                    rint_dict[ch].append(rint)
                    zint_dict[ch].append(zint)
                else:
                    rint_dict[ch] = [rint]
                    zint_dict[ch] = [zint]
        
        #print rint_dict
        #print zint_dict
        # use lists to store charges of the same value
        #   -> too many pylab.plot calls would slow everything down too much
        for ch,rints in rint_dict.iteritems():
            pylab.plot(rints, zint_dict[ch],'ro',markersize=abs(ch)**(1./3.)*plot_fact,markeredgecolor='r')
                
        print "Positive charges plotted. (%s)\n"%time.asctime()
                
        rint_dict = {}
        zint_dict = {}
        for charge in self.chc_neg:
            [r,x,y,z,ch,symb,at_inds] = charge
            xint = charge[plot_coors[0]]
            yint = charge[plot_coors[1]]
            zint = charge[axis]
            phi = math.atan2(yint, xint)
            
            rint = None
            if sector[0] <= phi <= sector[1]:
                rint = math.sqrt(xint**2. + yint**2.)
            elif (sector[0] <= phi + math.pi <= sector[1]) or (sector[0] <= phi - math.pi <= sector[1]):
                rint = -math.sqrt(xint**2. + yint**2.)
            if not rint == None:
                chabs = abs(ch)
                if chabs in rint_dict:
                    rint_dict[chabs].append(rint)
                    zint_dict[chabs].append(zint)
                else:
                    rint_dict[chabs] = [rint]
                    zint_dict[chabs] = [zint]
                    
        #print rint_dict
        #print zint_dict
        for chabs,rints in rint_dict.iteritems():
            pylab.plot(rints, zint_dict[chabs],'bo',markersize=chabs**(1./3.)*plot_fact,markeredgecolor='b')
                
        print "Negative charges plotted. (%s)\n"%time.asctime()
        
        #pylab.show()
        pylab.savefig(fig_name)
        
    def print_pchrg(self,indfile='pot_indices',formout='col',outfile=None):
        """
        Print out the point charges to a file.
        A file with the original indices of the atoms relating to the point charges is printed as well.
        """
        (x0,y0,z0) = self.orig
        
        out_list = self.chc_pos+self.chc_neg
        out_list.sort()
        num_out = len(out_list)
        
        if formout == 'col':
            if outfile == None: outfile = 'potential.xyz'
        elif formout == 'tmol':
            if outfile == None: outfile = 'pointcharges'
        else:
            print 'Output format %s not included'%formout
            exit(1)
        
        print 'Writing %s with %i point charges\n'%(outfile,num_out)

        opchrg = open(outfile,'w')
        if formout == 'col':
            opchrg.write(str(num_out)+'\n')
        elif formout == 'tmol':
            opchrg.write('$point_charges\n')
            
        indf = open(indfile,'w')
        indf.write(str(self.num_charg_in)+'\n')
        for ch_coll in out_list:
            (x,y,z) = (ch_coll[1]*self.scale[0] + x0, ch_coll[2]*self.scale[1] + y0, ch_coll[3]*self.scale[2] + z0)
            ch = ch_coll[4]
            if formout == 'col':
                opchrg.write(' %3s %12.8f %12.8f %12.8f %12.8f\n'%(ch_coll[5], ch, x, y, z))
            elif formout == 'tmol':
                opchrg.write(' %12.8f %12.8f %12.8f %12.8f\n'%(x, y, z, ch))
            # using only 5 decimal digits for the charge caused noticable numerical differences
            
            ind_str = ''
            for at_ind in ch_coll[6]:
                ind_str += ' %i'%(at_ind)
            indf.write(ind_str+'\n') 
            
        if formout == 'tmol':
            opchrg.write('$end\n')
        opchrg.close()

class grad_expand:
    """
    Class to expand a gradient.
    """ 
    def __init__(self,grad_file):
        self.grad_file = grad_file
        self.grad_out = '%s.expanded'%self.grad_file

    def grad_expand(self):
        """
        Expand the gradient according to pot_indices.
        """
        self.read_pot_old()

        poti = open('pot_indices','r')
        chinchk = int(poti.readline())
        assert(chinchk == self.num_charg_in)
        self.grad_list = [[] for i in xrange(self.num_charg_in)]

        old_grad = open(self.grad_file,'r')

        line_og = old_grad.readline()
        while(line_og!=''):
           line_pi = poti.readline()
           if lvprt>=1:
              print "input: ", line_og,
           #print line_pi,
           assert(line_pi!='')

           og_comps = [float(string.translate(word,DEtrans)) for word in line_og.split()]

           ch_inds = [int(word)-1 for word in line_pi.split()]
           curr_chs = [self.charges[ch_ind] for ch_ind in ch_inds]
           sum_ch = sum(curr_chs)

           for i,ch_ind in enumerate(ch_inds):
               fac = curr_chs[i] / sum_ch
               self.grad_list[ch_ind] = [fac*og_comp for og_comp in og_comps]
               if lvprt>=1: print '->  % 14.6E % 14.6E % 14.6E'%(self.grad_list[ch_ind][0],self.grad_list[ch_ind][1],self.grad_list[ch_ind][2])

           line_og = old_grad.readline()

        poti.close()
        old_grad.close()

    def read_pot_old(self,infile='pot_old.xyz'):
        """
        Read the point charges from pot_old.xyz.
        """
        print "Reading data ... (%s)"%(time.asctime())
        pchrg = open(infile,'r')
        line = pchrg.readline()
        self.num_charg_in = int(line.split()[0])
        print "%i Point charges read in. (%s)"%(self.num_charg_in,time.asctime())

        self.charges = []

        line = pchrg.readline()
        while(line!=''):
            words = line.split()
            ch = float(words[1])

            self.charges.append(ch)

            line = pchrg.readline()

        pchrg.close()
    
    def print_grad(self):
        """
        Print the expanded gradient.
        """
        new_grad = open(self.grad_out,'w')
        for comp in self.grad_list:
            new_grad.write('% 14.6E % 14.6E % 14.6E\n'%(comp[0],comp[1],comp[2]))
        new_grad.close()

if __name__=='__main__':
    import sys
    inpfile='pot_reduce.in'
    if len(sys.argv) == 1:
        print "No argument given, using reduction mode"
        mode = "reduce"
    else:
        mode = sys.argv[1]

    C = 0.0
    r0 = 0.
    group = 'c1'
    lvprt = 0
    sector = [0,0,0]
    formout = 'col'
    orig = [0.,0.,0.]
    scale = [1.,1.,1.]
    try:
        execfile(inpfile)
    except IOError:
        print "Error: Please create input file %s, example:"%inpfile
        print "# --- %s ---"%inpfile
        print "C=0.1\nr0=0."
        print "orig=[3.0,-4.7,5.0]\nformout='col'"
        exit(2)

    if mode == "reduce":
        print "Reducing a set of pointcharges"
        reduce_pot()
    elif mode == "expand":
        import string
        DEtrans=string.maketrans('D','E') # for reading Fortran doubleprecision format

        #grad_file = sys.argv[2]
        grad_file = "cartgrd.pointcharges"
        print "Expanding gradient %s"%grad_file
        if not group == 'c1':
           print "Expansion only supported for c1"
           exit(1)
        expand(grad_file=grad_file)
        print "Expansion finished"
    elif mode == "plot":
        print "Plotting pointcharges"
        import pylab
        plot_pchrg()
    else:
        print "Mode %s not supported"%s
        exit(3)

    bummer = open('bummer','w')
    bummer.write("normal termination\n")
    bummer.close()
   
