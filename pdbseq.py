import Bio.PDB
import gzip

class pdbres(str):
    def __new__(cls, seqstr, position=None):
        
        if seqstr:
            obj = str.__new__(cls, seqstr[0])
        else:
            #print 'input a residue name for class pdbres '
            obj = str.__new__(cls, seqstr)

        if position:
            obj.loc = position
            obj.hasloc = True
        else:
            obj.loc = None
            obj.hasloc = False
        return obj
    
class pdbseq(object):
    
    def __init__(self, seq, position_list=None):
        self.AA = []
        self.allhasloc = False
        self.partialloc = False
        self.seq = seq
        self.loc_list=[]
        
        if not position_list:
            for i in range(0,len(seq)):
                tmpres = pdbres(seq[i])
                self.AA.append(tmpres)
            self.loc_list = [None]*len(seq)
            
        elif len(seq) == len(position_list):
            if None in list(set(position_list)):
                if not [None] == list(set(position_list)):
                    self.partialloc = True
            else:
                self.allhasloc = True
                self.partialloc = True
            self.loc_list = position_list
            for i in range(0,len(seq)):
                tmpres = pdbres(seq[i],position_list[i])
                self.AA.append(tmpres)

        else:
            
            print 'when initializing pdbseq class, no partial location list allowed, treated as none'
            assert (len(seq) == len(position_list))
            for i in range(0,len(seq)):
                tmpres = pdbres(seq[i])
                self.AA.append(tmpres)

        
##    def getloc(self,pdb_ID='4A14',pdb_file_dir='/Users/xji3/clemensCode/CaOutput/pdb/'):
##        handle = gzip.open(pdb_file_dir + pdb_ID + '/' + pdb_ID + '.pdb.gz','r')
##        structure = Bio.PDB.PDBParser().get_structure(pdb_ID, handle)
##        model = structure[0]
##        return None

                
    def __getitem__(self,x): 
        newseq = self.seq[x]
        newloclist = self.loc_list[x]
        if isinstance(newloclist, int):
            newloclist = [self.loc_list[x]]
        ans = pdbseq(newseq,newloclist)
        return ans

    def __len__(self):
        if len(self.AA) == len(self.loc_list) == len(self.seq):
            return len(self.seq)
        else:
            print 'please check len() function in class pdbseq'
            return len(self.seq)

    def __add__(self,x):
        if isinstance(x, pdbres):
            
            #ans = new.instance(pdbseq(self.seq,self.loc_list))
            newseq=self.seq + x
            newloclist = self.loc_list + [x.loc]
            ans = pdbseq(newseq,newloclist)
            return ans
            
        elif isinstance(x, pdbseq):
            ans = pdbseq(self.seq,self.loc_list)
            for i in x.AA:
                ans = ans + i            
            return ans

        elif isinstance(x, str):
            newseq=self.seq + x
            newloclist = self.loc_list + [None]*len(x)
            ans = pdbseq(newseq,newloclist)
            return ans

        
    def __radd__(self,x):
        if isinstance(x, pdbres):
            newseq=x + self.seq
            newloclist = [x.loc] + self.loc_list
            ans = pdbseq(self.seq,self.loc_list)
            return ans
        elif isinstance(x,pdbseq):
            ans = pdbseq(self.seq,self.loc_list)
            for i in reversed(x.AA):
                ans = i + ans
            return ans
        elif isinstance(x, str):
            newseq= x + self.seq
            newloclist = [None]*len(x) + self.loc_list
            ans = pdbseq(newseq,newloclist)
            return ans
       

    def __str__(self):
        return ''.join(self.AA)

    def count(self,x):
        return self.seq.count(x)

    def find(self,x):
        if isinstance(x, str):
            return str.find(self.seq,x)
        elif isinstance(x,pdbseq):
            return str.find(self.seq,x.seq)
        else:
            print 'Please define the find function in class pdbseq for type', type(x)
            return None

    def __eq__(self,x):
        if isinstance(x, str):
            return self.seq == x
        elif isinstance(x, pdbseq):
            return (self.seq == x.seq and self.loc_list == x.loc_list
                    and self.AA == x.AA and self.allhasloc == x.allhasloc
                    and self.partialloc == x.partialloc)
        else:
            print 'Please define the __eq__ function in class pdbseq'
            return False

    def __ne__(self,x):
        if isinstance(x, str):
            return self.seq != x
        elif isinstance(x, pdbseq):
            return (self.seq != x.seq or self.loc_list != x.loc_list
                    or self.AA != x.AA or self.allhasloc != x.allhasloc
                    or self.partialloc != x.partialloc)
        else:
            print 'Please define the __eq__ function in class pdbseq'
            return False
##    def replace(x, y, z):
##        new_AA = [z if item==y else item for item in x.AA]
##        return 
