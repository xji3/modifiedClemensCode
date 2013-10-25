# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# make yield compatible with Python 2.2
from __future__ import generators

#cl: use numpy
#cl import gzip
import gzip
#from Numeric import array, sum, sqrt
from numpy import array, sum, sqrt, asmatrix
import tempfile
import os
import sys

from Bio.PDB import *
from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap

__doc__="""
Calculation of residue depth (using Michel Sanner's MSMS program for the
surface calculation).

Residue depth is the average distance of the atoms of a residue from 
the solvent accessible surface.

Residue Depth:

    rd=ResidueDepth(model, pdb_file)

    print rd[(chain_id, res_id)]

Direct MSMS interface:

    Typical use:

        surface=get_surface("1FAT.pdb")

    Surface is a Numeric array with all the surface 
    vertices.  

    Distance to surface:

        dist=min_dist(coord, surface)

    where coord is the coord of an atom within the volume
    bound by the surface (ie. atom depth).

    To calculate the residue depth (average atom depth
    of the atoms in a residue):

    rd=residue_depth(residue, surface)
"""

def _read_vertex_array(filename):
    """
    Read the vertex list into a Numeric array.
    """
    fp=open(filename, "rb")
    vertex_list=[]
    for l in fp.readlines():
        sl=l.split()
#cl: length is 10
#        if not len(sl)==9:
        if not len(sl)==10:
            # skip header
            continue
        vl=map(float, sl[0:3])
        vertex_list.append(vl)
    fp.close()
    return array(vertex_list)




def min_dist(coord, surface):
    """
    Return minimum distance between coord
    and surface.
    """
    d=surface-coord
    d2=sum(d*d, 1)
    return sqrt(min(d2))

def residue_depth(residue, surface):
    """
    Return average distance to surface for all
    atoms in a residue, ie. the residue depth.
    """
    atom_list=residue.get_unpacked_list()
    length=len(atom_list)
    d=0
    for atom in atom_list:
        coord=atom.get_coord()
        d=d+min_dist(coord, surface)
    return d/length

def ca_depth(residue, surface):
    if not residue.has_id("CA"):
        return None
    ca=residue["CA"]
    coord=ca.get_coord()
    return min_dist(coord, surface)

class ResidueDepth(AbstractPropertyMap):
    """
    Calculate residue and CA depth for all residues.
    """
    def __init__(self, model, pdb_file):
        
        depth_dict={}
        depth_list=[]
        depth_keys=[]
        self.terminate=False
        # get_residue
        residue_list=Selection.unfold_entities(model, 'R')
        # make surface from PDB file
        surface=self.get_surface(pdb_file)

        if not self.terminate:
        # calculate rdepth for each residue
            for residue in residue_list:
                if not is_aa(residue):
                    continue
                rd=residue_depth(residue, surface)
                ca_rd=ca_depth(residue, surface)
                # Get the key
                res_id=residue.get_id()
                chain_id=residue.get_parent().get_id()
                depth_dict[(chain_id, res_id)]=(rd, ca_rd)
                depth_list.append((residue, (rd, ca_rd)))
                depth_keys.append((chain_id, res_id))
                # Update xtra information
                residue.xtra['EXP_RD']=rd
                residue.xtra['EXP_RD_CA']=ca_rd
            AbstractPropertyMap.__init__(self, depth_dict, depth_keys, depth_list)
        else:
            return None

    def get_surface(self,pdb_file, PDB_TO_XYZR="./aux/pdb_to_xyzr", MSMS="./aux/msms"):
        """
        Return a Numeric array that represents 
        the vertex list of the molecular surface.

        PDB_TO_XYZR --- pdb_to_xyzr executable (arg. to os.system)
        MSMS --- msms executable (arg. to os.system)
        """
        if self.terminate==False:
            # extract xyz and set radii
            xyz_tmp=tempfile.mktemp()
            PDB_TO_XYZR=PDB_TO_XYZR+" %s > %s"
            make_xyz=PDB_TO_XYZR % (pdb_file, xyz_tmp)
            if os.system(make_xyz): #Xiang Check
                self.terminate=True
            else:
                os.system(make_xyz)
            # make surface
            surface_tmp=tempfile.mktemp()
            MSMS=MSMS+" -probe_radius 1.5 -if %s -of %s > "+tempfile.mktemp()
            make_surface=MSMS % (xyz_tmp, surface_tmp)
            if os.system(make_surface):#Xiang Check
                self.terminate=True
            else:
                os.system(make_surface) # might calculate twice, live with it for now
            surface_file=surface_tmp+".vert"
            # read surface vertices from vertex file
            if not self.terminate:
                surface=_read_vertex_array(surface_file)
                return surface
            else:
                return None
            # clean up tmp files
            # ...this is dangerous
            #os.system("rm "+xyz_tmp)
            #os.system("rm "+surface_tmp+".vert")
            #os.system("rm "+surface_tmp+".face")
        else:
            return None
        

#if __name__=="__main__":
#
#    import sys
#
#    p=PDBParser()
#    s=p.get_structure("X", sys.argv[1])
#    model=s[0]
#
#    rd=ResidueDepth(model, sys.argv[1])
#
#    for item in rd:
#        print item
