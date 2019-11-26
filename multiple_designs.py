'''
Author: Mudra Hegde
Email: mhegde@broadinstitute.org
This script runs the base editor design script for multiple BE type inputs
'''

import pandas as pd
import argparse
import os


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-file', type=str, help='File with Ensembl transcript IDs or fasta file with nucleotide sequences')
    parser.add_argument('--ip-type', type=str, help='Indicate if input file contains transcripts or sequences; tid for transcripts, nuc for nucleptides')
    parser.add_argument('--be-file', type=str, help='.txt file with list of BE types for design or individual run parameters')
    parser.add_argument('--be-type', type=str, help='be-file contains BE types or specific input parameters; be for be-types; '
                                                    'The file with parameters should have the columns, PAM, Window, Sgrna length, Edit, Intron buffer')
    parser.add_argument('--output-prefix', type=str, help='Prefix for output folders')
    return parser


if __name__ == '__main__':
    args = get_parser().parse_args()
    be_input = pd.read_table(args.be_file)
    for i,b in be_input.iterrows():
        print(b[0])
        output_prefix = args.output_prefix+'_'+b[0]
        if args.be_type == 'be':
            cmd = 'python base_editing_guide_designs.py --input-file '+args.input_file+' --input-type '+ args.ip_type +' --be-type '+b[0]+' --output-name '+output_prefix
        else:
            pam = b['PAM']
            window = b['Window']
            sglen = b['Sgrna length']
            edit = b['Edit']
            intron_buffer = b['Intron buffer']
            filter_gc = b['Filter GC']
            cmd = 'python base_editing_guide_designs.py --input-file '+args.input_file+' --input-type '+ args.ip_type +\
                  ' --pam '+pam+' --edit-window '+window+' --sg-len '+str(sglen)+' --edit '+edit+' --intron-buffer '\
                  +str(intron_buffer)+'--filter-gc '+filter_gc+ '--output-name '+output_prefix
        os.system(cmd)



