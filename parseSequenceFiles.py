'''
Created on 2018/10/08

@author: yonezawa
'''

import os

def parseSequences(aligned_dir):
    sequences = {}
    
    if aligned_dir[-1] != '/':
        aligned_dir += '/'
        
    aligned_files = os.listdir(aligned_dir)
    for aligned_file in aligned_files:
        elm = aligned_file.split('_')
        gene = '_'.join(elm[:4])
        sequences.setdefault(gene, {})
        
        for line in open(aligned_dir + aligned_file):
            if line[0] == '>':
                transcript = line.strip()[1:]
                sequences[gene].setdefault(transcript, '')
            else:
                sequences[gene][transcript] += line.strip()

    return sequences


if __name__ == '__main__':
    aligned_dir = 'aligned_contigs'