from Bio.SeqUtils import MeltingTemp as mt


class SequenceTools:
    def __init__(self):
        self.rcNucl_DICT = {}
        self.rcNucl_DICT['A'] = 'T'
        self.rcNucl_DICT['T'] = 'A'
        self.rcNucl_DICT['G'] = 'C'
        self.rcNucl_DICT['C'] = 'G'
        self.rcNucl_DICT['N'] = 'N'

    def calc_gc(self, sequence):
        return float(sequence.count('G') + sequence.count('C'))/len(sequence)
    
    def calc_tm(self, sequence):
        return mt.Tm_NN(sequence)

    def reverse_complementary(self,sequence):
        result = []
        for nucl in sequence[::-1]:
            result += [self.rcNucl_DICT[nucl]]
        return ''.join(result) 