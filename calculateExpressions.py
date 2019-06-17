'''
Created on 2018/09/12

@author: yonezawa
'''
import subprocess
import os

def executeKallisto(kallisto_dir, seqfile, readfile_pairs, bootstrap=100):
    if trinity_util_dir[-1] != '/':
        trinity_util_dir += '/'

    calcindex_list = ['kallisto', 'index', '-i', 'region_index', seqfile]
    subprocess.call(calcindex_list)    
    
    for i in range(len(readfile_pairs)):
        readfile1 = readfile_pairs[i][0]
        readfile2 = readfile_pairs[i][1]
        calcex_list = [kallisto_dir + 'kallisto', 'quant', '-i', 'region_index', '-o', 'condition{0:02d}_kallisto'.format(i + 1), '-b', bootstrap, readfile1, readfile2]
        subprocess.call(calcex_list)


def executeRSEM(trinity_util_dir, bowtie2_dir, seqfile, readfile_pairs, cpu=4):
    if trinity_util_dir[-1] != '/':
        trinity_util_dir += '/'
    if bowtie2_dir[-1] != '/':
        bowtie2_dir += '/'
    
    bowtie2_list = [bowtie2_dir + 'bowtie2-build', '-f', seqfile, seqfile]
    subprocess.call(bowtie2_list)
    
    for cond in readfile_pairs.keys():
        for rep in readfile_pairs[cond].keys():
            for i in range(len(readfile_pairs[cond][rep])):
                if os.path.exists(readfile_pairs[cond][rep][i]):
                    gunzip_list = ['gunzip', readfile_pairs[cond][rep][i]]
                    subprocess.call(gunzip_list)
            rsem_list = [trinity_util_dir + 'align_and_estimate_abundance.pl', '--seqType', 'fq', '--thread_count', str(cpu), '--left', readfile_pairs[cond][rep][0][:-3], '--right', readfile_pairs[cond][rep][1][:-3], '--transcripts', seqfile, '--est_method', 'RSEM', '--aln_method', 'bowtie2', '--prep_reference', '--output', 'rsem_{0:s}_{1:s}'.format(cond, rep)]
            subprocess.call(rsem_list)    
            copy_list = ['cp', 'rsem_{0:s}_{1:s}/RSEM.isoforms.results'.format(cond, rep), 'rsem_{0:s}_{1:s}.dat'.format(cond, rep)]
            subprocess.call(copy_list)
        
    abundance_list = [trinity_util_dir + 'abundance_estimates_to_matrix.pl', '--est_method', 'RSEM', '--gene_trans_map', 'none', '--out_prefix', 'abundance_matrix']
    for cond in readfile_pairs.keys():
        for rep in readfile_pairs[cond].keys():
            abundance_list.append('rsem_{0:s}_{1:s}.dat'.format(cond, rep))
        
    subprocess.call(abundance_list)


if __name__ == '__main__':
    seqfile = 'focused_ranges.fasta'
    trinity_util_dir = '/usr/local/bin/trinityrnaseq-Trinity-v2.6.6/util'
    bowtie2_dir = '/home/aogura/bowtie2-2.3.4.1-linux-x86_64/'
    readfile_pairs = [['trinity_out_dir/reads.left.fq.gz.P.qtrim.gz', 'trinity_out_dir/reads.right.fq.gz.P.qtrim.gz'], ['trinity_out_dir/reads2.left.fq.gz.P.qtrim.gz', 'trinity_out_dir/reads2.right.fq.gz.P.qtrim.gz']]
    cpu = 16
    
    executeRSEM(trinity_util_dir, bowtie2_dir, seqfile, readfile_pairs, cpu)
