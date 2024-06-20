from lib.BwaHandler import BwaHandler
from lib.FastaHandler import FastaHandler
from lib.SequenceTools import SequenceTools
from lib.Primer import Primer, PrimerPair
import os

class PCRSimulator:
    def __init__(self):
        self.sequenceTools = SequenceTools()

        self.primerPair_DICT = {}
        self.primerPair_LIST = []

        self.primer_DICT = {}
        self.primer_LIST = []

    def read_primerFile(self, fileName):
        print('reading primer ...')
        fin = open(fileName)
        for line in fin:
            if line.startswith('F.Primer'): continue
            fPrimerSeq, rPrimerSeq = line.rstrip('\n').upper().split('\t')[0:2]
            if not fPrimerSeq in self.primer_DICT:
                primer_gc = self.sequenceTools.calc_gc(fPrimerSeq)
                primer_tm = self.sequenceTools.calc_tm(fPrimerSeq)

                fPrimer = Primer(fPrimerSeq)
                fPrimer.gc = primer_gc
                fPrimer.tm = primer_tm
                
                self.primer_DICT[fPrimerSeq] = fPrimer
                self.primer_LIST += [fPrimer]
            else:
                fPrimer = self.primer_DICT[fPrimerSeq]

            if not rPrimerSeq in self.primer_DICT:
                primer_gc = self.sequenceTools.calc_gc(rPrimerSeq)
                primer_tm = self.sequenceTools.calc_tm(rPrimerSeq)

                rPrimer = Primer(rPrimerSeq)
                rPrimer.gc = primer_gc
                rPrimer.tm = primer_tm

                self.primer_DICT[rPrimerSeq] = rPrimer
                self.primer_LIST += [rPrimer]
            else:
                rPrimer = self.primer_DICT[rPrimerSeq]
            
            if not (fPrimerSeq, rPrimerSeq) in self.primerPair_DICT:
                primerPair = PrimerPair(fPrimer, rPrimer)
                self.primerPair_DICT[(fPrimerSeq, rPrimerSeq)] = primerPair
            else:
                primerPair = self.primerPair_DICT[(fPrimerSeq, rPrimerSeq)]

            self.primerPair_LIST += [primerPair]
        fin.close()
    
    def primer2fastq(self, fileName):
        print('Converting primer file to fastq format ...')
        fout = open(fileName, 'w')
        for primerSeq in self.primer_DICT.keys():
            fout.write('@' + str(primerSeq) + '\n')
            fout.write(primerSeq + '\n')
            fout.write('+' + '\n')
            fout.write('I'*len(primerSeq) + '\n')
        fout.close()

    def read_bwaResult(self, maxProductSize, fileName,):
        print('reading bwa result ...')
        
        fin = open(fileName)
        for line in fin:
            if line.startswith('@PG') == True:
                break
        for line in fin:
            data_LIST = line.rstrip('\n').split('\t')
            qname, flag, rname, pos = data_LIST[0:4]
            
            primerSeq = qname
            seqName = rname

            tPrimer = self.getPrimer(primerSeq)
            
            pos = int(pos)
            flag = int(flag)

            strand = '+'
            if flag&16 != 0:
                pos *= -1
                strand = '-'
            
            tPrimer.add_hit(seqName, strand, pos)

            for other in data_LIST[11:]:
                if other.startswith('XA') == False: continue
                for mapped in [seqName + ',' + str(pos)] + other[5:].rstrip(';').split(';'):
                    seqName, pos = mapped.split(',')[0:2]
                    pos = int(pos)

                    if pos > 0:
                        strand = '+'
                    else:
                        strand = '-'
                        pos *= -1
                    tPrimer.add_hit(seqName, strand, pos)
        fin.close()

        print('finding primerPair ...')
        for key, primerPair in self.primerPair_DICT.items():
            primerPair.find_hit(maxProductSize)


    def getPrimer(self, primerSeq):
        return self.primer_DICT[primerSeq]
    
    
    
    def writeFile(self, fileName):
        fout = open(fileName, 'w')

        for primerPair in self.primerPair_LIST:
            fout.write(primerPair.text() + '\n')
        fout.close()

    def debug(self):
        print(self.primerPair_LIST[0].fPrimer.seq)

    def run(self, threadN, maxMisN, maxProductSize, workingDir, refFile, inFile, outFile):
        if not os.path.exists(workingDir):
            os.makedirs(workingDir)

        self.read_primerFile(inFile)
        self.primer2fastq(workingDir + '/' + 'primer.fastq')
        
        #fastaHandler = FastaHandler(refFile)

        bwaHandler = BwaHandler()
        if bwaHandler.check_index(refFile) == False:
            bwaHandler.run_index(refFile)

        bwaHandler.run_aln(threadN, maxMisN, workingDir, refFile)
        
        self.read_bwaResult(maxProductSize, workingDir + '/' + 'primer.sam')

        self.writeFile(outFile)

if __name__ == '__main__':
    workingDir = 'tmp'
    refFile = '../db/Fxananassa_675_v1.0.fa'
    threadN = 128
    maxMisN = 1
    outFile = "out" 
    pcrSimulator = PCRSimulator()
    pcrSimulator.run(threadN, maxMisN, workingDir, refFile, outFile)