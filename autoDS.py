#!/usr/bin/python

"Automatically creates a PDB:DNA data set."
from optparse import OptionParser, OptionGroup
from data import *
#Xiaqng's change
##import optparse



#-------------------------------------------------------------------------------
# Command line options
#-------------------------------------------------------------------------------
optpars = OptionParser(usage="usage: %prog [options] filename",
                      version="%prog 0.1")

office_Mac_address='/Users/xji3/clemensCode'
#home_pc_address='G:/Dropbox/My Files/BRC/Small Project/Clemens PY code'
#laptop_address='E:/Dropbox/Dropbox/My Files/BRC/Small Project/Clemens PY code'
address=office_Mac_address

CnFolder='/CnDataOutput'
CaFolder='/CaOutput'
OutputFolder=CaFolder
if OutputFolder==CaFolder:
    CNorCA='Ca-Ca'
elif OutputFolder==CnFolder:
    CNorCA='C-N'
else:
    print
    print '============Please select the output file================'
    print
    
group = OptionGroup(optpars, "Mandatory options")
group.add_option("-d", "--dsdir", action="store", dest="ds_dir", default=(address+OutputFolder),
                   help='Mandatory: data-set output directory')
group.add_option("-f", "--file", action="store", dest="xmlfilename", default=(address+'/input/pdb_search.xml'),
                   help='Mandatory: XML file for PDB query')
group.add_option("-p", "--pisadir", action="store", dest="local_pisadir", default=(address+OutputFolder+'/pisa'),
                   help='Mandatory: local dir where Pisa files should be stored')
group.add_option("-c", "--criterion", action="store", dest="DistCriterion", default=(CNorCA),# Choose between C-N or Ca-Ca
                   help='Mandatory: distance criterion')
optpars.add_option_group(group)

#Xiaqng's change
(opts,args) = optpars.parse_args()
##print opts
##print args

mandatories = ['ds_dir','xmlfilename','local_pisadir']
for m in mandatories:
    if not opts.__dict__[m]:
        print "mandatory option is missing\n"
        optpars.print_help()
        exit(-1)
#################################################



optpars.add_option("-i", "--initpdb", action="store", dest="initpdb_filename", default=(address+'/input/init_pdbids.txt'),
                   help='Skip PDB search, read PDB-IDs from this (single-column) file instead')
optpars.add_option("-m", "--initmp", action="store", dest="initmp_filename", default=(address+'/input/init_mpids.txt'),
                   help='Skip mpstruc search, read PDB-IDs from this (single-column) file instead')
optpars.add_option("-o", "--initmono", action="store", dest="initmono_filename", default=(address+'/input/monomers.txt'),
                   help='Skip Pisa step, read monomer PDB-IDs from this (single-column) file instead')

##args=['-d','aaa','-f','bbb','-p','ccc']
##optpars.set_defaults(xmlfilename='aaa')
(options, args) = optpars.parse_args()

#-------------------------------------------------------------------------------
if __name__ == '__main__':

    mandatory = ['xmlfilename', 'local_pisadir', 'ds_dir']
    for m in mandatory:


        if not options.__dict__[m]:
            optpars.print_help()
            optpars.error('At least one mandatory option is missing.\n')

    #---------------------------------------------------------------------------
    # Initialize Data object
    #---------------------------------------------------------------------------
    if options.ds_dir[-1] != '/': options.ds_dir = options.ds_dir + '/'
    data = Data(options.initpdb_filename, options.xmlfilename, options.ds_dir,options.DistCriterion)
    
    #---------------------------------------------------------------------------
    # Apply additional filters (no membrane proteins, only monomers, ...)
    #---------------------------------------------------------------------------
    data.filterMembraneProteins(options.initmp_filename)
    data.filterNonMonomers(options.initmono_filename, options.local_pisadir)
    
    #---------------------------------------------------------------------------
    # Get the AA sequences for the PDB files and remove those with non-unique
    # chains
    #---------------------------------------------------------------------------
    data.fasta4pdb()
    data.removeNonUniqueChainPDBs()
    
    #---------------------------------------------------------------------------
    # Map PDB-IDs to CCDS-IDs (and keep track of all other mappings)
    #---------------------------------------------------------------------------
    data.mapPDB2CCDS()
    data.filterByCCDSalignments(50, 97.0)
    
    #---------------------------------------------------------------------------
    # Print results to files
    #---------------------------------------------------------------------------
    data.printResults()
    
#-------------------------------------------------------------------------------
