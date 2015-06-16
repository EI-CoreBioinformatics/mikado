import sys,argparse,os
sys.path.append(
                os.path.dirname(
                                os.path.dirname(__file__)
                                ))
import csv
from shanghai_lib.loci_objects.transcript import transcript
import shanghai_lib.exceptions
from shanghai_lib.loci_objects.abstractlocus import abstractlocus
from shanghai_lib.parsers.GTF import GTF
from shanghai_lib.parsers.GFF import GFF3
import bisect # Needed for efficient research

'''This is still an embryo. Ideally, this program would perform the following functions:

1- Define precision/recall for the annotation
2- Use a "flexible" mode to determine accuracy
3- Detect gene *fusions* as well as gene *splits*  

'''

class gene:
    pass




def main():
    
    def to_gff(string):
        '''Function to recognize the input file type and create the parser.'''
    
    parser=argparse.ArgumentParser('Tool to define the spec/sens of predictions vs. references.')
    input_files=parser.add_argument_group('Prediction and annotation files.')
    input_files.add_argument('-r', '--reference', type=to_gff, help='Reference annotation file.', required=True)
    input_files.add_argument('-p', '--prediction', type=to_gff, help='Prediction annotation file.', required=True)
    parser.add_argument('--distance', type=int, default=2000, 
                        help='Maximum distance for a transcript to be considered a polymerase run-on. Default: %(default)s')
    
    args=parser.parse_args()

    



if __name__=='__main__': main()