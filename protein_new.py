from __future__ import division
import Bio.PDB
from collections import defaultdict
import os.path, gzip
from geneticCode import *
import itertools, numpy, string, subprocess, re
import math as math
import residueDepth as ResDepth
from NACCESS import *
from pdbseq import *

#-------------------------------------------------------------------------------
class Protein:
#-------------------------------------------------------------------------------
    def __init__(self, id, aa, ds_dir,criterion):
        print 'PDB:', id,
        self.out_dir = ds_dir
        self.pdb_ID = id
        self.pdb_AA = aa # the fasta sequence
        self.pdb_structure_AA = pdbres('') # the structure seq aligned to the fasta seq
        self.ccds_match = []
        self.hasMatch = False
        self.crit=criterion
        self.Non_Stand_pos=defaultdict(list) # recording where non-standard AA or un-classified cases like 'H_PTR' is
        self.Strange_res=[]
        self.Excluded_res=[]
        self.parsePDBfile()
        self.special = False

#-------------------------------------------------------------------------------
    def printFastaSequences(self):
        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '_ccds.fasta.gz', 'w') as f_fasta:
            for a in self.ccds_match:
                f_fasta.write('>' + a.ccds_ID + '_protein\n')
                f_fasta.write(a.ccds_AA + '\n')

                f_fasta.write('>' + a.ccds_ID + '_dna\n')
                f_fasta.write(a.ccds_DNA + '\n')

#-------------------------------------------------------------------------------
    def printPDBfasta2StructureAlignment(self):
        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '_fasta2struct.fasta.gz', 'w') as f_fasta:
            f_fasta.write('>' + self.pdb_ID + '_fasta\n')
            f_fasta.write(self.pdb_AA.seq + '\n')

            f_fasta.write('>' + self.pdb_ID + '_pdb_structure_residues\n')
            f_fasta.write(self.pdb_structure_AA.seq + '\n')

#-------------------------------------------------------------------------------
    def printCCDS2PDBAlignments(self):
        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '_best_ungapped.fasta.gz', 'w') as f_fasta:
            ccds = self.ccds_match[0]
            
            tmpstr = ccds.ccds_AA[:ccds.ccds_local_alignment_start] + string.replace(ccds.ccds_alignedAA[:ccds.ungapped_segment_start], '-', '')
            idx = 3*len(tmpstr)
            offset = str(idx + 1)
            f_fasta.write('>' + ccds.ccds_ID + '_to_' + self.pdb_ID + '_fasta_longest_ungapped_segment' + '_Starting_DNA_CCDS_pos_' + offset + '\n')
            f_fasta.write(ccds.ccds_DNA[idx:idx+3*ccds.ungapped_segment_length] + '\n')

            offset = str(ccds.ungapped_segment_start + ccds.ccds_local_alignment_start + 1)
            f_fasta.write('>' + ccds.ccds_ID + '_to_' + self.pdb_ID + '_fasta_longest_ungapped_segment' + '_Starting_CCDS_pos_' + offset + '\n')
            f_fasta.write(ccds.ccds_alignedAA[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length] + '\n')

            i1 = string.find(self.pdb_AA.seq, ccds.pdb_alignedAA.seq[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length])
            i2 = i1 + ccds.ungapped_segment_length
            f_fasta.write('>' + self.pdb_ID + '_to_' + ccds.ccds_ID + '_longest_ungapped_segment' + '_Starting_aligned_PDB_struct_pos_' + str(i1 + 1) + '\n')
            f_fasta.write(self.pdb_structure_AA.seq[i1:i2] + '\n')

#-------------------------------------------------------------------------------
    def printLocalCCDSalignments(self):
        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '_fasta2CCDS.fasta.gz', 'w') as f_fasta:
            for ccds in self.ccds_match:
                offset = str(ccds.ccds_local_alignment_start + 1)
                f_fasta.write('>' + ccds.ccds_ID + '_locally_aligned_to_' + self.pdb_ID + '_fasta' + '_starting_CCDS_pos_' + offset + '\n')
                f_fasta.write(ccds.ccds_alignedAA + '\n')
    
                offset = str(ccds.pdb_local_alignment_start + 1)
                f_fasta.write('>' + self.pdb_ID + '_fasta_locally_aligned_to_' + ccds.ccds_ID + '_starting_fasta_pos_' + offset + '\n')
                f_fasta.write(ccds.pdb_alignedAA.seq + '\n')

#-------------------------------------------------------------------------------
    def printDMat(self):
        ccds = self.ccds_match[0]
        r_from = string.find(self.pdb_AA.seq, ccds.pdb_alignedAA.seq[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length])
        r_to = r_from + ccds.ungapped_segment_length
        self.printSubDMatrix(r_from, r_to)

#-------------------------------------------------------------------------------
    def printMismatches(self):
        ccds = self.ccds_match[0]
        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.mm.gz', 'w') as f:
            start = ccds.ungapped_segment_start
            end = start + ccds.ungapped_segment_length
            ctr = 0
            for k, v in sorted(ccds.mismatch_or_gap.items()):
                if k >= start and k < end:
                    ctr = ctr + 1
            f.write(str(ctr) + ' ' + self.pdb_ID + ' ' + ccds.ccds_ID + '\n')
            for k, v in sorted(ccds.mismatch_or_gap.items()):
                if k >= start and k < end:

                    f.write(str(k - ccds.ungapped_segment_start + 1) + '\t')
                    for i in range(len(v)-1):
                        if isinstance(v[i],str):
                            f.write(v[i] + '\t')
                        elif isinstance(v[i],pdbseq):
                            f.write(v[i].seq + '\t')
                        else:
                            f.write(v[i] + '\t')
                    if isinstance(v[i+1],str):
                        f.write(v[i+1] + '\n')
                    elif isinstance(v[i+1],pdbseq):
                        f.write(v[i+1].seq + '\n')
                    else:
                        f.write(v[i+1] + '\n')

#-------------------------------------------------------------------------------
    def printSubDMatrix(self, r_from, r_to):
        assert(r_from <= r_to)

        dmat = self.calcDistMatrix()

        # Test version
#        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '_' + 'v2.pw.gz', 'w') as f:
#            for i in range(r_from, r_to):
#                line = [ dmat[i, j] for j in range(r_from, r_to) ]
#                f.writelines([ "%8.4f\t" % i for i in line ])
#                f.write('\n')

        # Version for Jeff
        ccds = self.ccds_match[0] #alright, no loops here @Xiang
        tmpstr = ccds.ccds_AA[:ccds.ccds_local_alignment_start] + string.replace(ccds.ccds_alignedAA[:ccds.ungapped_segment_start], '-', '')
        idx = 3*len(tmpstr) # Maybe should record this location in the future, it's used several times @Xiang
        tmpstr = ccds.ccds_DNA[idx:idx+3*ccds.ungapped_segment_length]

        tmpstr2 = ccds.pdb_aligned_structure_AA[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length]
        nres = len(string.replace(tmpstr2.seq, '-', ''))

        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.pw.gz', 'w') as f:
            f.write(self.pdb_ID + ' ' + str(nres) + '\n')
            idx1 = 1
            for i in range(r_from, r_to-1): # not that necessary to use -1 for the loop here @Xiang range(r_to, r_to)=[], anyhow, it's good habit @Xiang
                idx2 = idx1+1
                for j in range(i+1, r_to):
                    if not math.isnan(dmat[i, j]):
                        f.write(str(idx1) + ' ' + str(idx2) + ' ' + tmpstr[(idx1-1)*3:(idx1-1)*3+3] + ' ')
                        f.write(tmpstr[(idx2-1)*3:(idx2-1)*3+3] + ' ' + str(dmat[i, j]) + '\n')
                    idx2 = idx2+1
                idx1 = idx1+1

#-------------------------------------------------------------------------------
    def printSeqFile(self):
        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.seq.gz', 'w') as f:
            ccds = self.ccds_match[0]
            f.write(self.pdb_ID + ' ' + str(ccds.ungapped_segment_length) + '\n')
            tmpstr = ccds.ccds_AA[:ccds.ccds_local_alignment_start] + string.replace(ccds.ccds_alignedAA[:ccds.ungapped_segment_start], '-', '')
            idx = 3*len(tmpstr)
            f.write(ccds.ccds_DNA[idx:idx+3*ccds.ungapped_segment_length] + '\n')
            f.write(ccds.ccds_alignedAA[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length] + '\n')        
            i1 = string.find(self.pdb_AA.seq, ccds.pdb_alignedAA.seq[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length])
            i2 = i1 + ccds.ungapped_segment_length
            f.write(self.pdb_structure_AA.seq[i1:i2] + '\n')

#-------------------------------------------------------------------------------
    def printSolvAcc(self):
        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.sa.gz', 'w') as f:
            ccds = self.ccds_match[0]
        
            pdb_file = self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.pdb'
            handle = gzip.open(pdb_file + '.gz', "r")
            structure = Bio.PDB.PDBParser().get_structure(self.pdb_ID, handle)
            model = structure[0]
            chain_list = [ a.get_id() for a in model ]
            chain = model[chain_list[0]]
        
            rd = ResDepth.ResidueDepth(chain, pdb_file + '.gz')

            if rd.terminate:
                print '===========Warning: MSMS or pdb_to_xyzr doesnot work'
                return False

            # CCDS string up until the beginning of the ungapped segment without gaps
            tmpstr = ccds.ccds_AA[:ccds.ccds_local_alignment_start] + string.replace(ccds.ccds_alignedAA[:ccds.ungapped_segment_start], '-', '')
            # First codon of ungapped sequence in DNA sequence
            idx = 3*len(tmpstr)
            # The DNA sequence for the ungapped segment
            strng = ccds.ccds_DNA[idx:idx+3*ccds.ungapped_segment_length]

            
            # Start of ungapped segment in the PDB fasta sequence
            r_from = self.pdb_AA.find( ccds.pdb_alignedAA[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length])
            # End of ungapped segment in the PDB fasta sequence
            r_to = r_from + ccds.ungapped_segment_length
            # The structure sequence is aligned to the fasta sequence, may have gaps in the ungapped segment
            nres = len(self.pdb_structure_AA[r_from:r_to]) - self.pdb_structure_AA[r_from:r_to].count('-')
            f.write(self.pdb_ID + ' ' + str(nres) + '\n')

            idx = 0

            #Xiang Version start
            nchck=0
            sortloc_list = self.pdb_structure_AA.loc_list[:]
            for item in rd:
                resseq=item[0].get_id()[1]
                
                index_list = self.returnindexes(sortloc_list,resseq)
                if index_list:
                    idx = index_list[0]
                    sortloc_list[idx] = None
                    #assert (len(index_list)==1)
                    #it should be one to one but like pdbID = 1CS8, resseq 1p to 96p then 1 to 749. It's converted to 1-96 followed by 1-749 in Bio.PDB 
                else:
                    continue

                
                if item[0].get_id()[0][0:2]=='H_' and item[0].get_id()[0]!='H_W': # exclude the non standard MSE case for now
                    if idx<r_to and idx >= r_from and self.pdb_structure_AA[idx] != '-':
                        nres = nres-1
                    else:
                        print 'non-standard AA happens at position ',idx, 'but out of comparison region. resseq position saved for nAccess analysis'
                        print
                    nonStandID=item[0].get_id()[0]
                    self.Non_Stand_pos[nonStandID].append(resseq)
                    self.Excluded_res.append(resseq)
                    self.special = True
                if idx<r_to and idx>=r_from:
                    if not item[0].has_id('CA'):
                        if  self.pdb_structure_AA[idx] != '-':
                            nres = nres-1
                        self.Strange_res.append(resseq)
                        self.special = True
                        #self.Excluded_res.append(resseq) only store Nonstandard for now

                    if item[0].has_id('CA') and item[0].get_id()[0] == ' ': # No Het
                        f.write(str(idx+1) + '\t' + strng[(idx-r_from)*3:(idx-r_from)*3+3] + '\t' + str(item[1][0]) + '\t' + str(item[1][1]) + '\t'
                                + str(self.pdb_structure_AA.seq[idx]) + '\t' + str(item[0].get_resname()) + '\t' + str(resseq) + '\n')
                        nchck = nchck + 1
                        
                
            #End of Xiang Version
            self.Excluded_res = list(set(self.Excluded_res))
            self.Strange_res = list(set(self.Strange_res))

            if nchck != nres:
                with gzip.open(self.out_dir + 'ErrorPDBList.solerr.gz', 'a') as f2:
                    print '===============ResidueDepth Checking Issue for PDB ID:', self.pdb_ID
                    print 'nchck = ', nchck, 'nres = ', nres
                    f2. write(str(self.pdb_ID) + '\t' + str(nchck) + '\t' + str(nres) + '\n')
            else:
                print
            
            #assert(nchck == nres)

#-------------------------------------------------------------------------------
    def printnAccess(self):        
        
        with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.sanac.gz', 'w') as f:
            ccds = self.ccds_match[0]
        
            pdb_file = self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.pdb'
            
            handle = gzip.open(pdb_file + '.gz', "r")
            structure = Bio.PDB.PDBParser().get_structure(self.pdb_ID, handle)
            model = structure[0]
            chain_list = [ a.get_id() for a in model ]
            chain = model[chain_list[0]]
            
            pdb_file_dir=self.out_dir + self.pdb_ID + '/'
            naccess_dir='/Users/xji3/nAccess/Naccess/naccess'
        
            rd = NACCESS(chain, pdb_file_dir, self.pdb_ID, naccess_dir )


            # CCDS string up until the beginning of the ungapped segment without gaps
            tmpstr = ccds.ccds_AA[:ccds.ccds_local_alignment_start] + string.replace(ccds.ccds_alignedAA[:ccds.ungapped_segment_start], '-', '')
            # First codon of ungapped sequence in DNA sequence
            idx = 3*len(tmpstr)
            # The DNA sequence for the ungapped segment
            strng = ccds.ccds_DNA[idx:idx+3*ccds.ungapped_segment_length]

            
            # Start of ungapped segment in the PDB fasta sequence
            r_from = self.pdb_AA.find( ccds.pdb_alignedAA[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length])
            # End of ungapped segment in the PDB fasta sequence
            r_to = r_from + ccds.ungapped_segment_length
            # The structure sequence is aligned to the fasta sequence, may have gaps in the ungapped segment
            nres = len(self.pdb_structure_AA[r_from:r_to]) - self.pdb_structure_AA[r_from:r_to].count('-')
            f.write(self.pdb_ID + ' ' + str(nres) + '\n')

            idx = 0

            #Xiang Version start
            nchck=0
            nres = nres - len(list(set([i for i in self.Excluded_res if i>=r_from and i<r_to])))
            sortloc_list = self.pdb_structure_AA.loc_list[:]
            for item in rd:
                resseq=item[0].get_id()[1]
                
                index_list = self.returnindexes(sortloc_list,resseq)
                if index_list:
                    idx = index_list[0]
                    sortloc_list[idx] = None
                    #assert (len(index_list)==1)
                    #it should be one to one but like pdbID = 1CS8, resseq 1p to 96p then 1 to 749. It's converted to 1-96 followed by 1-749 in Bio.PDB 
                else:
                    continue

                
                if idx<r_to and idx>=r_from:

                    if item[0].has_id('CA') and item[0].get_id()[0] == ' ': # No Het
                        f.write(str(idx+1) + '\t' + strng[(idx-r_from)*3:(idx-r_from)*3+3] + '\t' + str(item[1]['all_atoms_rel']) + '\t'
                                + str(item[1]['all_atoms_abs']) + '\t' + str(self.pdb_structure_AA.seq[idx]) + '\t' + str(item[0].get_resname()) + '\t' + str(resseq) + '\n')
                        nchck = nchck + 1
            if nchck != nres:
                with gzip.open(self.out_dir + 'ErrorPDBList.nacerr.gz', 'a') as f2:
                    print '===============nAccess Checking Issue for PDB ID:', self.pdb_ID
                    print 'nchck = ', nchck, 'nres = ', nres
                    f2. write(str(self.pdb_ID) + '\t' + str(nchck) + '\t' + str(nres) + '\n')
            else:
                print

            #End of Xiang Version
            #assert(nchck == nres)

#-------------------------------------------------------------------------------
    def printExcludedRes(self):
        if self.special:
            ccds = self.ccds_match[0]
            # CCDS string up until the beginning of the ungapped segment without gaps
            tmpstr = ccds.ccds_AA[:ccds.ccds_local_alignment_start] + string.replace(ccds.ccds_alignedAA[:ccds.ungapped_segment_start], '-', '')
            # First codon of ungapped sequence in DNA sequence
            idx = 3*len(tmpstr)
            # The DNA sequence for the ungapped segment
            strng = ccds.ccds_DNA[idx:idx+3*ccds.ungapped_segment_length]       
            # Start of ungapped segment in the PDB fasta sequence
            r_from = string.find(self.pdb_AA, ccds.pdb_alignedAA[ccds.ungapped_segment_start:ccds.ungapped_segment_start+ccds.ungapped_segment_length])
            # End of ungapped segment in the PDB fasta sequence
            r_to = r_from + ccds.ungapped_segment_length
            idx=0
            
            with gzip.open(self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.er.gz', 'w') as f:
                if self.Non_Stand_pos:
                    print '=============Warning: non-standard Amino Acid occuring'
                    print 'pdb ID: ', self.pdb_ID
                for nonStand in self.Non_Stand_pos:
                    print 'non-standard AA: ', nonStand, 'occurred at resseq positions: ', self.Non_Stand_pos[nonStand]
                    for resseq in self.Non_Stand_pos[nonStand]:
                        index_list = self.returnindexes(self.pdb_structure_AA,resseq)
                        idx = index_list[0]
                        assert (len(index_list)==1)
                        f.write(str(idx+1) + '\t' + strng[(idx-r_from)*3:(idx-r_from)*3+3] + '\t' + nonStand + '\t' + str(self.pdb_structure_AA.seq[idx]) + '\n')
                if self.Strange_res:
                    print '=============Warning: non-CA Detected'
                    print 'pdb ID: ', self.pdb_ID, ' at resseq position: ', self.Strange_res
                    for resseq in self.Strange_res:
                        index_list = self.returnindexes(self.pdb_structure_AA,resseq)
                        idx = index_list[0]
                        assert (len(index_list)==1)
                        f.write(str(idx+1) + '\t' + strng[(idx-r_from)*3:(idx-r_from)*3+3] + '\t' + 'non-CA' + '\t' + str(self.pdb_structure_AA.seq[idx]) + '\n')
                        
                    

#-------------------------------------------------------------------------------
    def printInfo(self):
        self.printFastaSequences()
        self.printPDBfasta2StructureAlignment()
        self.printCCDS2PDBAlignments()
        self.printLocalCCDSalignments()
        self.printDMat()
        self.printMismatches()
        self.printSeqFile()
        self.printSolvAcc()
        self.printExcludedRes()
        self.printnAccess()

#-------------------------------------------------------------------------------
    def returnindexes(self,a,b): # a is pdbseq instance @Xiang
        if isinstance(a,pdbseq):
            return [i for i,x in enumerate(a.loc_list) if x==b]
        elif isinstance(a,list):
            return [i for i,x in enumerate(a) if x==b]
        else:
            print 'Please Define returnindexes function in protein class'
            return None
        
       

#-------------------------------------------------------------------------------
    def alignCCDSlocal(self):
        for pdb2ccds in self.ccds_match: # here pdb2ccds is class PDB_CCDS, not the dict in class data class @Xiang
            out = self.alignLocalWater(pdb2ccds.ccds_AA) # use water for local alignment
            iter = out.__iter__()
            # Parse water's output
            for line in iter:
                if line.strip() and line.split() and line.split()[0] == 'asis':
                    break
            aligned_subseq1 = line.split()
            pdb2ccds.pdb_local_alignment_start = int(aligned_subseq1[1]) - 1
            s1 = pdb2ccds.pdb_alignedAA = pdbseq(aligned_subseq1[2]) # PDB sequence
            iter.next() # skip the vertical bars
            aligned_subseq2 = iter.next().split()
            assert(len(aligned_subseq1[2]) == len(aligned_subseq2[2]))
            pdb2ccds.ccds_local_alignment_start = int(aligned_subseq2[1]) - 1
            s2 = pdb2ccds.ccds_alignedAA = aligned_subseq2[2] # CCDS sequence
            algnmt_length = pdb2ccds.local_alignment_length = len(aligned_subseq2[2])

            # Keep track of the gaps also in the actual structure sequence (not just the fasta of the structure)
            s3 = self.pdb_structure_AA[pdb2ccds.pdb_local_alignment_start:]
            gaps = [ i.start() for i in re.finditer('-', s1.seq) ] # record all gap positions in s1 @Xiang
            #string.replace(s3, '-', '') #comment for now @Xiang
            for i in gaps:
                s3 = self.insertInString(i, s3, '-')
            s3 = pdb2ccds.pdb_aligned_structure_AA = s3[:algnmt_length] # PDB file sequence

            # Find longest ungapped segment, mismatches, gaps, percent identity as compared to fasta seq
            seg_len = longest_seg_len = mtch_ctr = ccds_gap_ctr = 0
            for idx in range(algnmt_length):
                p_res = s1[idx]
                c_res = s2[idx]
                if p_res == '-' or c_res == '-':
                    if seg_len > longest_seg_len:
                        longest_seg_len = pdb2ccds.ungapped_segment_length = seg_len
                        pdb2ccds.percent_identity = 100 * mtch_ctr/longest_seg_len
                        # if we hit a gap we don't need to add one to get the start pos
                        pdb2ccds.ungapped_segment_start = idx - seg_len
                    # Gap => reset counters
                    mtch_ctr = seg_len = 0
                elif p_res == c_res:
                    seg_len += 1
                    mtch_ctr += 1
                else: # mismatch
                    seg_len += 1
            if seg_len > longest_seg_len:
                longest_seg_len = pdb2ccds.ungapped_segment_length = seg_len
                pdb2ccds.ungapped_segment_start = 1 + idx - seg_len
                pdb2ccds.percent_identity = 100 * mtch_ctr/longest_seg_len

            # In longest ungapped segment find mismatches, gaps, percent identity as compared to pdb file seq
            ccds_gap_ctr = 0
            for idx in range(algnmt_length):
                p_res = s3[idx]
                c_res = s2[idx]
                tmp = 3*(pdb2ccds.ccds_local_alignment_start + idx - ccds_gap_ctr)
                codon = pdb2ccds.ccds_DNA[tmp:tmp+3]
                if p_res == '-' or c_res == '-':
                    if c_res == '-':
                        ccds_gap_ctr = ccds_gap_ctr+1
                        # { aligned AA sequences offset : [ pdb_aa_res, ccds_aa_res, ccds_dna_triplet, triplet transl. ] }
                        pdb2ccds.mismatch_or_gap[idx] = [ p_res.seq, '-', '---', '-' ]
                    else:
                        pdb2ccds.mismatch_or_gap[idx] = [ '-', c_res, codon, universalCode[codon] ]
                elif p_res != c_res or c_res != universalCode[codon]: # mismatch or translation exception
                    pdb2ccds.mismatch_or_gap[idx] = [ p_res.seq, c_res, codon, universalCode[codon] ]
                

        # Sort the ccds alignments by longest ungapped segment length & percent identity
        # Sorting algorithm is stable; order of the first sort is preserved if second criterium is equal
        self.ccds_match = sorted(self.ccds_match, key=lambda prot: prot.percent_identity, reverse=True)
        self.ccds_match = sorted(self.ccds_match, key=lambda prot: prot.ungapped_segment_length, reverse=True)
        # test in XiangDraft3.py @Xiang

#-------------------------------------------------------------------------------
    def insertInString(self, idx, seq, insert):
        return seq[:idx] + pdbseq(insert) + seq[idx:]

#-------------------------------------------------------------------------------
    def alignLocalWater(self, ccds_AA):
        water_localfile_dir='/Users/xji3/emBoss/EMBOSS-6.6.0/emboss'
        if os.path.isfile(water_localfile_dir+'/water'):
            water = water_localfile_dir+'/water'
        else:
            water = 'water'
        return subprocess.check_output([
                        water,
                        '-stdout',
                        '-auto',
                        '-awidth3=100000',
                        '-asequence=asis:' + self.pdb_AA.seq,
                        '-bsequence=asis:' + ccds_AA
                    ]).split('\n')

#-------------------------------------------------------------------------------
    def parsePDBfile(self):
        pdb_file = self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.pdb'
        handle = gzip.open(pdb_file + '.gz', "r")
        structure = Bio.PDB.PDBParser().get_structure(self.pdb_ID, handle)

        model = structure[0]
        chain_list = [ a.get_id() for a in model ]
        print 'Using chain:', chain_list[0]
        chain = model[chain_list[0]]

        if self.crit=='C-N':
            ppb = Bio.PDB.PPBuilder()
            print 'Using C-N distance criterion'
            # pdBuilder() uses C-N distance criterion
        elif self.crit=='Ca-Ca':
            ppb = Bio.PDB.CaPPBuilder()
            print 'Using Ca-Ca distance criterion'
            # CaPPBuilder() uses Ca-Ca distance criterion
        else:
            ppb = Bio.PDB.CaPPBuilder()
            print 'Using Ca-Ca distance criterion'
            # CaPPBuilder() uses Ca-Ca distance criterion
            # Use Ca-Ca criterion by default
        

        
        # Include non-standard residues
        # more info. in biopdb_faq.pdf, page 11 @Xiang
        tmpstructure = [ pp for pp in ppb.build_peptides(chain, aa_only=False) ] # used for tracking position
        tmploc = []
        tmp = [ str(pp.get_sequence()) for pp in ppb.build_peptides(chain, aa_only=False) ]

        assert(not self.pdb_structure_AA)
        self.pdb_structure_AA = pdbseq('-' * len(self.pdb_AA)) #assign all gaps
        edges = [None] * len(tmpstructure)
        multiple = []
        
        i = 0
        for bb in tmpstructure:
            # Only consider fragments that occur only once
            # This doesn't cause a lot of loss since only very short (~2 res) fragments occur more than once
            a = str(bb.get_sequence())
            loclist = [c.get_id()[1] for c in bb]
            tmploc.append(loclist)
            n = self.pdb_AA.count(a)
            if(n == 1):
                idx = string.find(self.pdb_AA, a) # return the first element in a's notion in pdb_AA, could be 0 @Xiang
                self.pdb_structure_AA = self.pdb_structure_AA[:idx] + pdbseq(a,loclist) + self.pdb_structure_AA[idx+len(a):]
                edges[i] = [idx, idx+len(a)]
            elif n > 1 :
                multiple.append(i)
            else:
                print 'There is no match of this fragment in the pdb AA sequence', a # May need more operation here @Xiang
            i = i+1

        # Give fragments we couldn't place unambiguously before a second try
        for i in multiple:
            b = e = -1
            if i == 0 and edges[i+1]:
                b = 0
                e = edges[i+1][0]
            elif i == (len(tmp)-1) and edges[i-1]:
                b = edges[i-1][1]
                e = len(tmp)-1
            elif edges[i-1] and edges[i+1]:
                b = edges[i-1][1]
                e = edges[i+1][0]
            else:
                print 'I don\'t know where to place this fragment. It appears >= 2x in the structure.'
                print 'It will be ignored in the distance matrix:', tmp[i]
                continue
            if self.pdb_AA[b:e].count(tmp[i]) == 1:
                idx = string.find(self.pdb_AA[b:e], tmp[i])
                self.pdb_structure_AA = self.pdb_structure_AA[:b+idx] + pdbseq(tmp[i],tmploc[i]) + self.pdb_structure_AA[b+idx+len(tmp[i]):]
            else:
                print 'I don\'t know where to place this fragment. It appears >= 2x in the structure.',
                print 'It will be ignored in the distance matrix:', tmp[i]

#-------------------------------------------------------------------------------
    def calcDistMatrix(self):

        pdb_file = self.out_dir + self.pdb_ID + '/' + self.pdb_ID + '.pdb'
        handle = gzip.open(pdb_file + '.gz', "r")
        structure = Bio.PDB.PDBParser().get_structure(self.pdb_ID, handle)
        model = structure[0]
        chain_list = [ a.get_id() for a in model ]
        chain = model[chain_list[0]]

        if self.crit=='C-N':
            ppb = Bio.PDB.PPBuilder()
            # Include non-standard residues
            print 'Using C-N distance criterion for DistMatrix Calculation'
        elif self.crit=='Ca-Ca':

            ppb = Bio.PDB.CaPPBuilder()
            print 'Using Ca-Ca distance criterion for DistMatrix Calculation'
            # CaPPBuilder() uses Ca-Ca distance criterion
            #@Xiang
        else:
            ppb = Bio.PDB.CaPPBuilder()
            print 'Using Ca-Ca distance criterion for DistMatrix Calculation'
            # Use Ca-Ca criterion by default
        
        tmp = [ str(pp.get_sequence()) for pp in ppb.build_peptides(chain, aa_only=False) ]

        # Join the lists to get nRes
        structure_residues = [None] * len(list(itertools.chain(*tmp)))
        idx = 0
        for r in chain: #.get_residues() seems like r=residues in chain, which doesn't need that function @Xiang
            if r.has_id('CA') and r.get_id()[0] == ' ': # No Het
                structure_residues[idx] = r
                idx = idx + 1

        dmat = numpy.empty((len(self.pdb_AA), len(self.pdb_AA)), numpy.float)
        dmat[:] = numpy.NAN
        nres = len(structure_residues)
        r_idx1 = 0
        tot = len(self.pdb_AA)
        for i in range(tot-1): #why -1 range(3)=[0,1,2] Stop Codon? @Xiang 
            if self.pdb_structure_AA[i] != '-' and structure_residues[r_idx1]:
                dmat[i, i] = 0.0
                r1 = structure_residues[r_idx1]
                r_idx1 = r_idx1 + 1
            else:
                continue
            r_idx2 = r_idx1
            for j in range(i+1, tot):
                if self.pdb_structure_AA[j] != '-' and structure_residues[r_idx2]:
                    r2 = structure_residues[r_idx2]
                    r_idx2 = r_idx2 + 1
                else:
                    continue
                if (r1.has_id('CA') and r2.has_id('CA')):
                    dmat[i, j] = self.calcCAdist(r1, r2)
                    dmat[j, i] = dmat[i, j]
        return dmat

#-------------------------------------------------------------------------------
    def printDMatrix(self):
        tot = len(self.pdb_AA)
        for i in range(tot):
            for j in range(tot):
                print "%6.3f" % self.dmat[i, j], ' ',
            print
        print

#-------------------------------------------------------------------------------
    def calcCAdist(self, a, b):
        return a['CA'] - b['CA']

#-------------------------------------------------------------------------------
    def checkAlignmentThresholds(self, min_alignment_length = 50, min_pct_identity = 97):
        self.hasMatch = False
        for i in self.ccds_match:
            if (i.ungapped_segment_length >= min_alignment_length) and (i.percent_identity >= min_pct_identity):
                i.aboveThresholds = True
                self.hasMatch = True
        return self.hasMatch

#-------------------------------------------------------------------------------
class PDB_CCDS:
#-------------------------------------------------------------------------------
    def __init__(self, id, seqs):
        self.ccds_ID = id
        self.ccds_AA = seqs[0]
        self.ccds_DNA = seqs[1]

        # { ccds_AA sequence index : [ ccds_AA_res, ccds_DNA_triplet, translated_AA ] }
        # only exceptions are recorded
        # To some extent this is redundant, because translation exceptions in the ungapped segment are
        # now also recorded in mismatch_or_gap
        self.translationExceptions = defaultdict(list)
        self.checkSequences()
        self.pdb_alignedAA = '' # CCDS specific, missing residues in structure gapped out after alignment
        self.pdb_aligned_structure_AA = '' # as above but w. gaps for residues that are missing in the structure
        self.ccds_alignedAA = ''
        self.pdb_local_alignment_start = -1  # index, 0-based
        self.ccds_local_alignment_start = -1  # AA sequence index, 0-based
        self.local_alignment_length = 0
        self.ungapped_segment_length = 0
        self.ungapped_segment_start = -1  # index, 0-based
        self.percent_identity = 0.0

        # { aligned AA sequences offset : [ pdb_aa_res, ccds_aa_res, ccds_dna_triplet, ccds_dna_triplet_translation ] }
        # the triplet translation should match the ccds_aa_res; else it should also be reflected in translationExceptions
        self.mismatch_or_gap = defaultdict(list)

        # Flag for selective printing of entries; all entries are sorted by seg length
        # entries with equal seg lengths are ordered by percent identity
        aboveThresholds = False

#-------------------------------------------------------------------------------
    def checkSequences(self):
        assert(len(self.ccds_AA)*3+3 == len(self.ccds_DNA)) # They all end with a stop codon
        j = 0
        for i in range(len(self.ccds_AA)):
            aa_res = self.ccds_AA[i]
            codon = self.ccds_DNA[j:j+3]
            if universalCode[codon] != aa_res:
                self.translationExceptions[i] = [ aa_res, codon, universalCode[codon] ]
            j = j+3

#-------------------------------------------------------------------------------
