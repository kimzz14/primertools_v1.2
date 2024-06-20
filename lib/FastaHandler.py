import os

class FastaHandler:
    def __init__(self, fastaName):
        self.fastaName = fastaName

        self.index_DICT = {}
        
        self.rcNucl_DICT = {}
        self.rcNucl_DICT['A'] = 'T'
        self.rcNucl_DICT['T'] = 'A'
        self.rcNucl_DICT['G'] = 'C'
        self.rcNucl_DICT['C'] = 'G'
        self.rcNucl_DICT['N'] = 'N'

        if os.path.exists(self.fastaName + '.idx') == False:
            print("Index file dosen't exist!")
            print("Making index file ...")
            self.build_index(self.fastaName)
        
        self.open_fasta(self.fastaName)
        self.read_index(self.fastaName)

    def build_index(self, fastaFile):
        fout_seq = open(fastaFile + '.seq', 'w')

        fin = open(fastaFile)
        ref_DICT = {}
        seqName_LIST = []
        for line in fin:
            if line.startswith('>') == True:
                seqName = line[1:].rstrip('\n').split('\t')[0].split(' ')[0]
                ref_DICT[seqName] = 0
                seqName_LIST += [seqName]
            else:
                sequence = line.rstrip('\n').upper()
                ref_DICT[seqName] += len(sequence)
                fout_seq.write(sequence)
        fin.close()
        fout_seq.close()
        
        fout_index = open(fastaFile + '.idx', 'w')
        pos = 0
        for seqName in seqName_LIST:
            fout_index.write(seqName + '\t' + str(ref_DICT[seqName]) + '\t' + str(pos) + '\n')
            pos += ref_DICT[seqName]
        fout_index.close()

    def read_index(self, fastaFile):
        self.index_DICT = {}
        fin = open(fastaFile + '.idx')
        for line in fin:
            seqName, size, pos = line.rstrip('\n').split('\t')
            self.index_DICT[seqName] = (int(size), int(pos))
        fin.close()
    
    def open_fasta(self, fastaFile):
        self.fin = open(fastaFile + '.seq')

    def reverse_complementary(self,sequence):
        result = []
        for nucl in sequence[::-1]:
            result += [self.rcNucl_DICT[nucl]]
        return ''.join(result) 

    def get_seq(self, seqName, strand, sPos, ePos):
        if not seqName in self.index_DICT:
            return None
        seqSize, pos = self.index_DICT[seqName]
        if sPos < 1 or seqSize < ePos:
            return None

        self.fin.seek(pos + sPos - 1)
        sequence = self.fin.read(ePos - sPos + 1)

        if strand == '+':
            pass
        elif strand == '-':
            sequence = self.reverse_complementary(sequence)
        else:
            return None
            
        return sequence

if __name__ == '__main__':
    fastaHandler = FastaHandler('../../db/Fxananassa_675_v1.0.fa')
    seq = fastaHandler.get_seq('Fvb1-1', '-', 1, 100)
    print(seq)