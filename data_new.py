from __future__ import division
from protein_new import *
from pdbseq import *
import ftplib, gzip, os, StringIO, urllib2
from xml.dom.minidom import parseString
# Instructions for simple XML parsing found at
# http://www.travisglines.com/web-coding/python-xml-parser-tutorial
from collections import defaultdict


#-------------------------------------------------------------------------------
class Data:
#-------------------------------------------------------------------------------
    def __init__(self, init_pdblist, init_xml, ds_dir, DistCriterion):
        self.out_dir = ds_dir
        self.pdb_url = 'http://www.rcsb.org/pdb/rest/search'
        self.pdb_ftp = 'ftp.wwpdb.org'
        self.pdb_ftp_basedir = '/pub/pdb/data/structures/divided/pdb'
        self.mpr_url = 'http://blanco.biomol.uci.edu/mpstruc/listAll/mpstrucTblXml'
        self.pisa_basedir = '/pub/databases/msd/pisa/data'
        self.pisa_ftp = 'ftp.ebi.ac.uk'
        self.pdbIDs = []
        self.proteins = {} # { pdbID : Protein }
        self.pdb_seqs = defaultdict(list)
        self.ccds_seqs = defaultdict(list) # { ccdsID : [ AA seq, DNA seq ] }
        self.pdbid2uniprot = defaultdict(list)
        self.uniprot2pdb = defaultdict(list)
        self.uniprot2geneid = defaultdict(list)
        self.ccds_table = defaultdict(list)
        self.pdb2ccds = defaultdict(list)
        self.crit=DistCriterion

        if init_pdblist is None:
            self.queryPDB(init_xml)
        else:
            print "List of PDB-IDs read from file: ", init_pdblist
            self.pdbIDs = open(init_pdblist).read().splitlines()
        if self.pdbIDs:
            print "Found number of initial PDB entries:", len(self.pdbIDs)
        else:
            print "Failed to retrieve initial results."
        print

#-------------------------------------------------------------------------------
    def printResults(self):
        print 'Writing output files and calculating distance matrices & solvent accessibilities.'
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        total_pdb_num=len(self.pdbIDs)
        finishednum=1
        with gzip.open(self.out_dir + '/pdb/ErrorPDBList.solerr.gz', 'w') as f:
            f.write( 'ResidueDepth Checking Issue for [PDB ID nchck nres]' + '\n')
        with gzip.open(self.out_dir + '/pdb/ErrorPDBList.nacerr.gz', 'w') as f2:
            f2.write( 'nAccess Checking Issue for [PDB ID nchck nres]' + '\n')
            
        with open(self.out_dir + 'PDB_ID_list_final.txt', 'w') as f_pdb:
            f_pdb.write(str(len(self.proteins.keys())) + '\n')
            for pdbid, p in self.proteins.items():
                
                print '\t', pdbid, '\t', 'Calculation done:',finishednum/total_pdb_num*100,'%'
                finishednum = finishednum+1

                f_pdb.write(pdbid + '\n')
                p.printInfo()
                with gzip.open(self.out_dir + 'pdb/' + pdbid + '/' + pdbid + '.map.gz', 'w') as id_pdb:
                    for a in self.pdbid2uniprot[pdbid]:
                        for b in self.uniprot2geneid[a]:
                            ccdsids = []
                            for ccds in self.ccds_table[b]:
                                if ccds[0:4] == 'CCDS':
                                    ccdsids.append(ccds)
                            id_pdb.write(pdbid + '\t' + a + '\t' + b + '\t' + ','.join(ccdsids) + '\n')

#-------------------------------------------------------------------------------
    def queryPDB(self, xml):
        f = open(xml, 'r')
        queryXML = f.read() # the complete xml query string
        f.close()
        # Extract query filters and write to screen
        dom = parseString(queryXML)
        descriptions = dom.getElementsByTagName("description")
        print "Querying PDB ..."
        print "\nApplying the following filters:"
        for id in descriptions:
            print "\t*", id.toxml().replace('<description>','').replace('</description>','')
        print
        req = urllib2.Request(self.pdb_url, data=queryXML)
        f = urllib2.urlopen(req)
        self.pdbIDs = f.read().splitlines()
        f.close()

#-------------------------------------------------------------------------------
    def filterMembraneProteins(self, init_mplist):
        pdbIDs_membrane = []
        if init_mplist is None:
            print "No list of membrane proteins provided. Querying mpstruc ..."
            pdbIDs_membrane = self.queryMPSTRUC()
        else:
            print "List of membrane proteins read from file: ", init_mplist
            pdbIDs_membrane = open(init_mplist).read().splitlines()

        if pdbIDs_membrane:
            print 'Found number of membrane protein IDs:', len(pdbIDs_membrane)
        else:
            print 'Failed to retrieve membrane protein IDs.'
        print

        pdbIDs_membrane = list(set(self.pdbIDs) & set(pdbIDs_membrane))

        print 'The list of PDB-IDs contains', len(pdbIDs_membrane), 'membrane proteins:'
        for id in pdbIDs_membrane:
            print '\t' + id
        print

        for i in pdbIDs_membrane:
            self.pdbIDs.remove(i)
        print '=====|', len(pdbIDs_membrane), 'PDB-IDs removed;',
        print len(self.pdbIDs), 'remaining.', '|=====\n'

#-------------------------------------------------------------------------------
    def filterNonMonomers(self, init_monolist, local_pisadir):
        pdbIDs_monomers = []
        if init_monolist is None:
            print "No list of monomeric proteins provided. Querying Pisa ..."
            print 'Using Pisa to find most probable assembly (n-mer):'
            self.queryPisa(local_pisadir, pdbIDs_monomers)
        else:
            print "List of monomeric proteins read from file: ", init_monolist
            pdbIDs_monomers = open(init_monolist).read().splitlines()

        if pdbIDs_monomers:
            print 'Found number of (probably) monomeric-protein IDs:', len(pdbIDs_monomers)
        else:
            print 'Failed to retrieve monomeric-protein IDs.'
        print
        pdbIDs_nomonos = list(set(self.pdbIDs) - set(pdbIDs_monomers))
        pdbIDs_monomers = list(set(self.pdbIDs) & set(pdbIDs_monomers))

        print 'The list of PDB-IDs contains', len(pdbIDs_monomers), 'monomeric proteins:'

        for i in pdbIDs_nomonos:
            self.pdbIDs.remove(i)
        print '=====|', len(pdbIDs_nomonos), 'PDB-IDs removed;',
        print len(self.pdbIDs), 'remaining.', '|=====\n'

#-------------------------------------------------------------------------------
    def queryMPSTRUC(self):
        file = urllib2.urlopen(self.mpr_url)
        d = file.read()
        file.close()

        # Parse the xml
        dom = parseString(d)

        # Retrieve xml tags (<tag>data</tag>)
        membrane_ids = dom.getElementsByTagName("pdbCode")
        for i, id in enumerate(membrane_ids):
            membrane_ids[i] = id.toxml().replace('<pdbCode>','').replace('</pdbCode>','')
        return membrane_ids

#-------------------------------------------------------------------------------
    def queryPisa(self, local_pisadir, monomer_ids):
        ftp_open = False
        remove_id = []
        for pdbf in self.pdbIDs:
            s = pdbf.lower()
            pisa_subdir = s[1:3]
            local_dir = local_pisadir + '/' + pisa_subdir
            local_filename = local_dir + '/' + s + '_assembly.xml.gz'
            if not os.path.exists(local_filename):
                if ftp_open is False:
                    print 'Trying to download Pisa files from', self.pisa_ftp
                    try:
                        ftp = ftplib.FTP(self.pisa_ftp)
                        print ftp.login()
                        ftp.cwd(self.pisa_basedir)
                        ftp_open = True
                    except ftplib.all_errors, e:
                        errorcode_string = str(e).split(None, 1)
                        print 'Remote server Warning:', errorcode_string[1]
                        exit()
                if not os.path.exists(local_dir):
                    os.makedirs(local_dir)
                rmote_filename = pisa_subdir + '/' + s + '/' + s + '_assembly.xml.gz'
                try:
                    ftp.retrbinary('RETR ' + rmote_filename, open(local_filename, 'wb').write)
                except ftplib.all_errors, e:
                    errorcode_string = str(e).split(None, 1)
                    print 'Remote server warning:', errorcode_string[1], rmote_filename
                    print '\tRemoving PDB-ID', pdbf, 'from list.'
                    remove_id.append(pdbf)
                    os.remove(local_filename)
                    continue
            f = open(local_filename)
            compresseddata = f.read()
            f.close()
            buf = StringIO.StringIO(compresseddata)
            gzipper = gzip.GzipFile(fileobj=buf)
            tmp = gzipper.read()
            dom = parseString(tmp)
            asms = dom.getElementsByTagName("total_asm")
            mmsizes = dom.getElementsByTagName("mmsize")
            mmsize = '?'
            if mmsizes:
                mmsize = str(dom.getElementsByTagName("mmsize")[0].toxml().replace('<mmsize>','').replace('</mmsize>',''))
            elif asms:
                asm = str(asms[0].toxml().replace('<total_asm>','').replace('</total_asm>',''))
                if asm == '0':
                    mmsize = '1'
            if mmsize == '1':
                monomer_ids.append(pdbf)
            else:
                os.remove(local_filename)
        if ftp_open is True:
            print 'Closing FTP connection.'
            ftp.close()
        # Finally remove the IDs we said we would remove due to Pisa/FTP issues
        for i in remove_id:
            self.pdbIDs.remove(i)

#-------------------------------------------------------------------------------
    def fasta4pdb(self):
        url = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz'
        local_file = 'pdb_seqres.txt.gz'
        if not os.path.isfile(local_file):
            self.downloadFile(url, local_file)
        f = gzip.open(local_file, 'rb')
        seqname = True
        id_chn = []
        for line in f:
            if line.strip(): # makes sure empty lines are not included
                if seqname:
                    id_chn = line.strip().split()
                    id_chn = id_chn[0].replace('>', '').upper().split('_')
                    seqname = False
                else:
                    if id_chn[0] in self.pdbIDs:
                        self.pdb_seqs[id_chn[0]].append(pdbseq(line.strip()))
                    seqname = True
        f.close()

#-------------------------------------------------------------------------------
    def removeNonUniqueChainPDBs(self):
        remove_id = []
        print 'PDB IDs with more than one chain:\n\t',
        n_mtoc = 0
        for i in self.pdb_seqs:
            if len(self.pdb_seqs[i]) > 1:
                n_mtoc += 1
                print i + ' (' + str(len(self.pdb_seqs[i])) + ')',
            seq = list(set([a.seq for a in self.pdb_seqs[i]]))
            if len(seq) > 1:
                remove_id.append(i)
            self.pdb_seqs[i] = [self.pdb_seqs[i][0]]
        print '\nTotal:', n_mtoc, '(of', str(len(self.pdbIDs)) + ' remaining PDB IDs)'
        
        # Remove PDB-IDs with non-unique chains
        print 'Of these, the following PDB entries have non-unique chains:'
        for i in remove_id:
            print '\t' + i
            self.pdbIDs.remove(i)
        print '=====|', len(remove_id), 'PDB-IDs removed;', len(self.pdbIDs), 'remaining.', '|=====\n'

#-------------------------------------------------------------------------------
    def mapPDB2CCDS(self):
        self.mapPDB2UniProt()
        self.mapUniProt2GeneID()
        all_geneids = []
        self.cleanupMappingPDB2UniProt2GeneID(all_geneids)
        all_ccds = []
        self.readCCDStable(all_ccds, all_geneids)
        rm_ccds = []
        self.getSeqs4CCDS(all_ccds, rm_ccds)
        self.updateMappingDicts(all_geneids, rm_ccds)
        self.directlyMapPDB2CCDS()
        self.alignPDB2CCDSlocal()

#-------------------------------------------------------------------------------
    def mapPDB2UniProt(self):
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst'
        local_file = 'pdb_chain_uniprot.lst'
        if not os.path.isfile(local_file):
            self.downloadFile(url, local_file)
        f = open(local_file, 'rb')
        for line in f:
            if line.strip(): # makes sure empty lines are not included
                tmp = line.strip().split()
                p_id = tmp[0].upper()
                if p_id in self.pdbIDs:
                    self.pdbid2uniprot[p_id].append(tmp[2])
        f.close()
        
        # Use an additional source
        url = 'http://www.uniprot.org/uniprot/?query=organism:9606+AND+database:pdb&format=tab&compress=yes&columns=id,database(PDB)'
        local_file = 'uniprot-organism%3A9606+AND+database%3Apdb.tab.gz'
        if not os.path.isfile(local_file):
            self.downloadFile(url, local_file)
        f = gzip.open(local_file, 'rb')
        for line in f:
            if line.strip(): # makes sure empty lines are not included
                tmp = line.strip().split()
                up_id = tmp[0]
                plist = tmp[1].split(';')
                plist.pop(-1)
                for p_id in plist:
                    if p_id.upper() in self.pdbIDs:
                        self.pdbid2uniprot[p_id.upper()].append(up_id)
        f.close()

        # Remove duplicates
        for e in self.pdbid2uniprot:
            self.pdbid2uniprot[e] = list(set(self.pdbid2uniprot[e]))

        # Reverse mapping -> UniProtKB : PDB
        for key in self.pdbid2uniprot: # we know we only have pdb-IDs that are in self.pdbIDs
            for val in self.pdbid2uniprot[key]: # list of uniprot values for a pdb entry
                self.uniprot2pdb[val].append(key)

#-------------------------------------------------------------------------------
    def mapUniProt2GeneID(self):
        url = 'http://www.uniprot.org/uniprot/?query=organism:9606+AND+database:GeneID&format=tab&compress=yes&columns=id,database(GeneID)'
        local_file = 'uniprot-organism%3A9606+AND+database%3AGeneID.tab.gz'
        if not os.path.isfile(local_file):
            self.downloadFile(url, local_file)
        f = gzip.open(local_file, 'rb')
        for line in f:
            if line.strip(): # makes sure empty lines are not included
                tmp = line.strip().split()
                up_id = tmp[0]
                # Only use the UniProt ID if we need it
                if up_id in self.uniprot2pdb:
                    glist = tmp[1].split(';')
                    glist.pop(-1)
                    for g_id in glist:
                        self.uniprot2geneid[up_id].append(g_id)
        f.close()

#-------------------------------------------------------------------------------
    def cleanupMappingPDB2UniProt2GeneID(self, all_geneids):
        multiple_uniprot = []
        multiple_geneids = []
        no_geneid = []
        no_uniprot = []
        print 'Non-unique UniProt mapping found for:'
        for i in self.pdbIDs:
            if i in self.pdbid2uniprot:
                up_list = self.pdbid2uniprot[i]
        
                # Check if we have multiple UniProt IDs
                if len(up_list) > 1:
                    print '\t' + i + ':', self.pdbid2uniprot[i]
                    multiple_uniprot.append(i)
        
                # Check if we have a GeneID for at least one of the UniProt IDs
                # and only keep the UniProt Ids for which we do
                haveid = set(up_list) & set(self.uniprot2geneid)
                havent = set(up_list) - haveid
                for no_id in havent:
                    self.pdbid2uniprot[i].remove(no_id)
                if haveid:
                    geneids = []
                    for j in up_list:
                        geneids = geneids + self.uniprot2geneid[j]
                        for k in geneids:
                            all_geneids.append(k)
                    if len(set(geneids)) > 1:
                        multiple_geneids.append(i)
                else:
                    no_geneid.append(i)
            else:
                no_uniprot.append(i)
        
        # All GeneIDs up until this point
        for i in list(set(all_geneids)):
            all_geneids.append(i)
        all_geneids=list(set(all_geneids))

        print 'Multiple GeneIDs found for (PDB-ID: UniProt-ID: [ \'GeneIDs\' ]):'
        for i in multiple_geneids:
            print '\t' + i + ':'
            for j in self.pdbid2uniprot[i]:
                print '\t\t' + j + ':', self.uniprot2geneid[j]
        print 'No Uniprot-ID found for PDB-ID:'
        for i in no_uniprot:
            print '\t' + i
            self.pdbIDs.remove(i)
        print '=====|', len(no_uniprot), 'PDB-IDs removed;', len(self.pdbIDs), 'remaining.', '|=====\n'
        print 'No GeneID found for PDB-ID:'
        no_geneid = list(set(no_geneid))
        for i in no_geneid:
            print '\t' + i
            del self.pdbid2uniprot[i]
            self.pdbIDs.remove(i)
        print '=====|', len(no_geneid), 'PDB-IDs removed;', len(self.pdbIDs), 'remaining.', '|=====\n'

#-------------------------------------------------------------------------------
    def readCCDStable(self, all_ccds, all_geneids):
        tmp_ccds = []
        # 0          1            2    3       4
        # chromosome nc_accession gene gene_id ccds_id ccds_status	cds_strand	cds_from	cds_to	cds_locations	match_type
        url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt'
        local_file = 'CCDS.current.txt'
        if not os.path.isfile(local_file):
            self.downloadFile(url, local_file)
        f = open(local_file, 'rb')
        for line in f:
            if line.strip(): # makes sure empty lines are not included
                s = line.strip().split()
                gene_id = s[3]
                if gene_id in all_geneids:

                    chromosome = s[0]
                    ccds_id = s[4]

                    # Chromosome will be first in the list, followed by the CCDS-IDs
                    if not self.ccds_table[gene_id]:
                        self.ccds_table[gene_id].append(chromosome)

                    self.ccds_table[gene_id].append(ccds_id)
                    tmp_ccds.append(ccds_id)
        f.close()
        for i in list(set(tmp_ccds)):
            all_ccds.append(i)
        all_ccde=list(set(all_ccds))

#-------------------------------------------------------------------------------
    def getSeqs4CCDS(self, all_ccds, rm_ccds):
        # Read sequences into one-line strings
        url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS_protein.current.faa.gz'
        local_file = 'CCDS_protein.current.faa.gz'
        if not os.path.isfile(local_file):
            self.downloadFile(url, local_file)
        f = gzip.open(local_file, 'rb')
        self.readCCDSseqs(all_ccds, f, 0)
        f.close()

        url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS_nucleotide.current.fna.gz'
        local_file = 'CCDS_nucleotide.current.fna.gz'
        if not os.path.isfile(local_file):
            self.downloadFile(url, local_file)
        f = gzip.open(local_file, 'rb')
        self.readCCDSseqs(all_ccds, f, 1)
        f.close()

        print len(all_ccds), 'of unique CCDS-IDs mapped to', len(self.pdbIDs), 'PDB-IDs.'
        print 'Of these,', len(self.ccds_seqs), 'sequences found in CCDS protein fasta file.'

        # Remove the CCDS that have not been found
        for i in list(set(all_ccds) - set(self.ccds_seqs.keys())):
            rm_ccds.append(i)

        print 'CCDS sequences not found for CCDS-IDs:'
        print '\t',
        for i in rm_ccds:
            print i,
            all_ccds.remove(i)
        print

#-------------------------------------------------------------------------------
    def readCCDSseqs(self, all_ccds, filehandle, idx):
        assert(idx >= 0 and idx <= 1) # ccds_seqs[key][0]=pro seq, ccds_seqs[key][1]=nucleotide seq@Xiang
        key = ''
        readme = False
        for line in filehandle:
            if line.strip(): # makes sure empty lines are not included
                line = line.rstrip('\n')

                if line[0] == '>':
                    line = line.split('|')
                    key = line[0].replace('>', '')
                    if key in all_ccds:
                        readme = True
                    else:
                        readme = False
                elif readme:
                    if key in self.ccds_seqs and len(self.ccds_seqs[key]) == (idx+1):
                        self.ccds_seqs[key][idx] = self.ccds_seqs[key][idx] + line
                    else:
                        self.ccds_seqs[key].append(line)

#-------------------------------------------------------------------------------
    def updateMappingDicts(self, all_geneids, rm_ccds):
        rm_pdb = []
        for i in self.pdbIDs: # i = PDB ID
            uniprot_no_ccds_id = []

            # Loop over UniProtIDs
            for j in self.pdbid2uniprot[i]: # j = UniProt ID
                geneid_no_ccds_id = []

                # Loop over GeneIDs
                for k in self.uniprot2geneid[j]: # k = GeneID
                    if k in self.ccds_table:

                        # Loop over ccds_aa
                        ccds_no_seq = []
                        for c in self.ccds_table[k]:
                            if c in rm_ccds:
                                ccds_no_seq.append(c)
                        for c in ccds_no_seq:
                            self.ccds_table[k].remove(c)

                    # Do we have any CCDS-IDs left?
                    # defaultdict returns [] if key does not exists
                    # Otherwise this would be a key error
                    if not self.ccds_table[k]:
                        geneid_no_ccds_id.append(k)
                if geneid_no_ccds_id:
                    # Remove geneids without CCDS
                    for k in geneid_no_ccds_id:
                        self.uniprot2geneid[j].remove(k)
                        all_geneids.remove(k)

                # If this leaves us with an empty list for a UniProtID
                # store UniProtID and remove it from the dictionary
                # This should also handle the case where multiple PDBs map to the same
                # UniProtID
                if not self.uniprot2geneid[j]:
                    uniprot_no_ccds_id.append(j)
                    del self.uniprot2geneid[j]
                    self.pdbid2uniprot[i].remove(j)
            if not self.pdbid2uniprot[i]:
                rm_pdb.append(i)
                del self.pdbid2uniprot[i]

        print 'No CCDS-ID or CCDS sequence found for PDB-ID:'
        for i in rm_pdb:
            self.pdbIDs.remove(i)
            print '\t' + i
        print '=====|', len(rm_pdb), 'PDB-IDs removed;', len(self.pdbIDs), 'remaining.', '|=====\n'

#-------------------------------------------------------------------------------
    def directlyMapPDB2CCDS(self):
        for i in self.pdbIDs:
            for up in self.pdbid2uniprot[i]:
                for gd in self.uniprot2geneid[up]:
                    for ccds in self.ccds_table[gd]:
                        if ccds[0:4] == 'CCDS':
                            if not ccds in self.pdb2ccds[i]:
                                self.pdb2ccds[i].append(ccds)

#-------------------------------------------------------------------------------
    def alignPDB2CCDSlocal(self):
        print '\nConstructing all local pairwise PDB:CCDS amino acid Alignments (using the program \'water\' from the EMBOSS tools).'
        print 'Recording longest ungapped segment of alignment and its percent identity, as well as',
        print 'all mismatches & gaps in the entire local alignment.\n'
        self.fetchPDBs()
        for pdb in self.pdb2ccds:
            # Create a new instance of protein
            self.proteins[pdb] = Protein(pdb, self.pdb_seqs[pdb][0], self.out_dir + 'pdb/',self.crit)

            # Assign CCDS
            for ccds in self.pdb2ccds[pdb]:
                self.proteins[pdb].ccds_match.append(PDB_CCDS(ccds, self.ccds_seqs[ccds])) # create a new instance of PDB_CCDS

            # Align PDB2CCDS
            self.proteins[pdb].alignCCDSlocal()

#-------------------------------------------------------------------------------
    def fetchPDBs(self):
        ftp_open = False
        remove_id = []
        for pdbf in self.pdbIDs:
            s = pdbf.lower()
            subdir = s[1:3]
            local_dir = self.out_dir + '/pdb/' + pdbf
            local_filename = local_dir + '/' + pdbf + '.pdb.gz'
            if not os.path.exists(local_filename):
                if ftp_open is False:
                    print 'Trying to download PDB files from', self.pdb_ftp
                    try:
                        ftp = ftplib.FTP(self.pdb_ftp)
                        print ftp.login()
                        ftp.cwd(self.pdb_ftp_basedir)
                        ftp_open = True
                    except ftplib.all_errors, e:
                        errorcode_string = str(e).split(None, 1)
                        print 'Remote server Warning:', errorcode_string[1]
                        exit()
                if not os.path.exists(local_dir):
                    os.makedirs(local_dir)
                rmote_filename = subdir + '/pdb' + s + '.ent.gz'
                try:
                    print '\t' + pdbf
                    ftp.retrbinary('RETR ' + rmote_filename, open(local_filename, 'wb').write)
                except ftplib.all_errors, e:
                    errorcode_string = str(e).split(None, 1)
                    print 'Remote server warning:', errorcode_string[1], rmote_filename
                    print '\tRemoving PDB-ID', pdbf, 'from list.'
                    remove_id.append(pdbf)
                    os.remove(local_filename)
                    continue
        if ftp_open is True:
            print 'Closing FTP connection.'
            ftp.close()

#-------------------------------------------------------------------------------
    def filterByCCDSalignments(self, min_alignment_length = 50, min_pct_identity = 97):
        print 'Filtering alignments by longest ungapped segment length and percent identity:\n',
        print '\tMinimum length of ungapped segment in local alignment =', min_alignment_length
        print '\tMinimum percent identity over length of this segment =', min_pct_identity
        remove_id = []
        for k, v in self.proteins.iteritems(): #k is pdbID, v is class protein @Xiang
            if not v.checkAlignmentThresholds(min_alignment_length, min_pct_identity):
                remove_id.append(k)

        if remove_id:
            print '\nThe following PDB-IDs did not pass the PDB:CCDS alignment filter for any CCDS and will be removed.'
            for i in remove_id:
                print '\t' + i
                self.pdbIDs.remove(i)
                del self.proteins[i]
                del self.pdbid2uniprot[i]
                del self.pdb2ccds[i]
            print '=====|', len(remove_id), 'PDB-IDs removed;', len(self.pdbIDs), 'remaining.', '|=====\n'

#-------------------------------------------------------------------------------
    def downloadFile(self, url, local_filename):
        try:
            print 'Downloading', local_filename, 'from:'
            print url
            req = urllib2.Request(url)
            f = urllib2.urlopen(req)
            local_file = open(local_filename, "wb")    
            local_file.write(f.read())
            local_file.close()
            f.close()
        except urllib2.HTTPError, e:
            print "HTTP Error:",e.code , url
        except urllib2.URLError, e:
            print "URL Error:",e.reason , url

#-------------------------------------------------------------------------------
