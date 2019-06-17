'''
Created on 2018/10/09

@author: yonezawa
'''

def geneTransMap(seqfile, result):
    writefile = open(result, 'w')
    for line in open(seqfile):
        if line[0] == '>':
            elm = line.strip()[1:].split('_')
            gene = '_'.join(elm[:4])
            transcript = '_'.join(elm)
            
            writefile.write('{0:s}\t{1:s}\n'.format(gene, transcript))
            
    writefile.close()


if __name__ == '__main__':
    seqfile = 'focused_ranges.fasta'
    result = seqfile + '.gene_trans_map'
    
    geneTransMap(seqfile, result)