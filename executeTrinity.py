'''
Created on 2018/10/03

@author: yonezawa
'''

import subprocess

def executeTrinity(trinity_dir, sample_file, max_memory='10G', cpu=4, trimmomatic_flag=True):
    if trinity_dir[-1] != '/':
        trinity_dir += '/'
    
    trinity_list = [trinity_dir + 'Trinity', '--seqType', 'fq', '--max_memory', max_memory, '--samples_file', sample_file, '--CPU', str(cpu), '--output', 'trinity_out_dir']
    if trimmomatic_flag:
        trinity_list.append('--trimmomatic')
        
    subprocess.call(trinity_list)
    
if __name__ == '__main__':
    trimmomatic_flag = True
    max_memory = '8G'
    cpu = 4
    trinity_dir = '/Users/yonezawa/tools/trinityrnaseq-2.0.5'
    sample_file = 'sample_file.txt'
    
    executeTrinity(trinity_dir, sample_file, max_memory, cpu, trimmomatic_flag)