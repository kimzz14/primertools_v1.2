from lib.PrimerDesigner import PrimerDesigner
from lib.PCRSimulator import PCRSimulator
from lib.Parameter import Parameter
import argparse


parser = argparse.ArgumentParser(description="Commands:")
subparsers = parser.add_subparsers(dest='commandType', help='Sub-command help')

# Parser for 'type1' command
parser_type1 = subparsers.add_parser('design', help='design PCR primer from user-entered sequence')
parser_type1.add_argument('-i', '--input', metavar="input.fa", type=str, help='input file name [fasta format]')
parser_type1.add_argument('-o', '--out', metavar="output.primer", type=str, help='output file name [primer format]')

# Parser for 'type2' command
parser_type2 = subparsers.add_parser('pcr', help='run PCR simulation based on reference sequences')
parser_type2.add_argument('-r', '--ref', metavar="test.fa", type=str, help='reference file name [fasta format]')
parser_type2.add_argument('-i', '--input', metavar="input.primer", type=str, help='input file name [primer format]')
parser_type2.add_argument('-o', '--out', metavar="output.primer", type=str, help='output file name [primer format]')

# Parser for 'type3' command
parser_type3 = subparsers.add_parser('info', help='calcuate primer length, GC content, Tm value')
parser_type3.add_argument('-i', '--input', metavar="input.primer", type=str, help='input file name [primer format]')
parser_type3.add_argument('-o', '--out', metavar="output.primer", type=str, help='output file name [primer format]')


args = parser.parse_args()

if args.commandType == 'design':
    if args.input == None or args.out == None:
        parser_type1.print_help()
        print()
    else:
        inFile = args.input
        outFile = args.out

        parameter = Parameter()
        primerDesigner = PrimerDesigner()
        primerPair_LIST = primerDesigner.gernate_primer_from_fasta(parameter, inFile)
        primerDesigner.write_file(outFile, primerPair_LIST)

elif args.commandType == 'pcr':
    if args.ref == None or args.input == None or args.out == None:
        parser_type2.print_help()
        print()
    else:
        refFile = args.ref #'../db/Fxananassa_675_v1.0.fa'
        inFile = args.input
        outFile = args.out #"out" 

        workingDir = 'tmp'
        threadN = 128
        maxMisN = 1
        maxProductSize = 2000
        pcrSimulator = PCRSimulator()
        pcrSimulator.run(threadN, maxMisN, maxProductSize, workingDir, refFile, inFile, outFile)

elif args.commandType == 'info':
    #print(f"Handling type2 with argument: {args.arg2}")
    parser_type2.print_help()
    print()
else:
    parser.print_help()