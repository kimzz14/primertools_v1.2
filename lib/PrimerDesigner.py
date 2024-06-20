#from FastaHandler import FastaHandler
from lib.Primer import Primer, PrimerPair
from lib.SequenceTools import SequenceTools
from lib.Parameter import Parameter

class PrimerDesigner:
    def __init__(self):
        self.sequenceTools = SequenceTools()

    def gernate_primer_from_fasta(self, parameter, fileName):
        ref_DICT = {}
        seqName_LIST = []

        fin = open(fileName)
        for line in fin:
            if line.startswith('>') == True:
                seqName = line.rstrip('\n')[1:].split('\t')[0].split(' ')[0]
                seqName_LIST += [seqName]
                ref_DICT[seqName] = []
            else:
                sequence = line.rstrip('\n').upper()
                ref_DICT[seqName] += [sequence]
        fin.close()
        primerPair_LIST = []
        for seqName in seqName_LIST:
            sequence = ''.join(ref_DICT[seqName])
            primerPair_LIST += self.gernate_primer_from_sequence(parameter, seqName, sequence)
        
        return primerPair_LIST

    def gernate_primer_from_sequence(self, parameter, seqName, sequence):
        fPrimer_LIST = []
        rPrimer_LIST = []
        primerPair_LIST = []

        sequence = sequence.upper()
        seq_length = len(sequence)

        for idx in range(1, seq_length):
            for primer_len in range(parameter.min_primer_len, parameter.max_primer_len + 1):
                sPos = idx
                ePos = idx + primer_len - 1
                if ePos > seq_length: break

                primer_seq = sequence[sPos - 1:ePos]
                
                #Filter N base
                if primer_len != (primer_seq.count('A') + primer_seq.count('T') + primer_seq.count('G') + primer_seq.count('C')): continue

                #Filter GC
                primer_gc = self.sequenceTools.calc_gc(primer_seq)
                if primer_gc < parameter.min_primer_GC or parameter.max_primer_GC < primer_gc: continue

                #Filter Tm
                primer_tm = self.sequenceTools.calc_tm(primer_seq)
                if primer_tm < parameter.min_primer_Tm or parameter.max_primer_Tm < primer_tm: continue
                
                
                #Filter is_gc_at_3_prime
                if parameter.is_gc_at_3_prime == True and (primer_seq.endswith('A') or primer_seq.endswith('T')): 
                    pass
                else:
                    fPrimer = Primer(primer_seq)
                    fPrimer.gc = primer_gc
                    fPrimer.tm = primer_tm
                    fPrimer.add_hit(seqName, '+', sPos)
                    fPrimer_LIST += [(fPrimer, sPos)]

                primer_seq_rev = self.sequenceTools.reverse_complementary(primer_seq)
                if parameter.is_gc_at_3_prime == True and (primer_seq_rev.endswith('A') or primer_seq_rev.endswith('T')): 
                    pass
                else:
                    rPrimer = Primer(primer_seq_rev)
                    rPrimer.gc = primer_gc
                    rPrimer.tm = primer_tm
                    rPrimer.add_hit(seqName, '-', sPos)
                    rPrimer_LIST += [(rPrimer, ePos)]
        
        for fPrimer, sPos in fPrimer_LIST:
            for rPrimer, ePos in rPrimer_LIST:
                product_len = ePos - sPos + 1
                if product_len < parameter.min_product_len or parameter.max_product_len < product_len: continue
                if parameter.max_primer_TmDiff < abs(fPrimer.tm - rPrimer.tm) : continue
                
                primerPair_LIST += [PrimerPair(fPrimer, rPrimer)]
        return primerPair_LIST

    def write_file(self, fileName, primerPair_LIST):
        fout = open(fileName, 'w')
        legend_LIST = ['F.Primer', 'R.Primer', 'F.Len', 'R.Len', 'F.Tm', 'R.Tm', 'F.GC', 'R.GC', 'ProductN', 'seqName', 'start', 'end', 'size']
        fout.write('\t'.join(map(str, legend_LIST)) + '\n')
        for primerPair in primerPair_LIST:
            primerPair.find_hit(100000)
            fout.write(primerPair.text() + '\n')
        fout.close()

if __name__ == '__main__':
    seqName = 'Fvb1-1'
    sequence = 'AGGTGATTTTAAATCTTTCAGTGAGGAGTATAAGATAGATATGGTGGAACCATTCAATTCAATAAAAAAAAAGATGTCAGTGCAAGTGGCTCTTCCTGGTAGTGATATGTCTCGAGCATTCTGCAAAGGCGCATCAGAAATAGTTTTAGGAATGTGTGACAAGGTTGTCAATACTGATGGTGAAGCTTTGCCTTTGTCTGAAGAGCAAAGAAACAAAATATCTGATGTCATAAATGGTTTTTCTTGTGAAGCTTTAAGAACTCTATGCGTAGCTTTCAAGGATATAGAGAGCCCCTCTGGTGCGGAAAGTATTCCTGAAGATGGATATACATTAATAGCTGTTGTGGGAATTAAGGATCCTGTGCGCCCTGCCGTGAGGGAAGCAGTGAAGACTTGTTTAGATGCTGGGATCACTGTGCGGATGGTCACTGGTGATAATATCAATACAGCTAAAGCCATAGCTAAAGAGTGTGGCATTTTGACAGAAGATGGTTTAGCTATAGAAGGGCCAGATTTCCGAAAAATGAGCGAACAAGAGAT'
    
    parameter = Parameter()
    primerDesigner = PrimerDesigner()

    primerDesigner.generate_primer(parameter, seqName, sequence)
    primerDesigner.write_file("test.primer")