'''
Created on 2016/05/30

@author: yonezawa
'''

def parseExonFiles (expressionfiles, file_format='kallisto', threshold=1, seq_dir=''):
    from math import log
    import sys, os
    
    if len(expressionfiles) <= 1:
        sys.stderr.write('Input more than one expression files.')
        sys.exit(-1)
    
    ex_values = {}
    log_ex = {}
    sequences = {}

    if len(seq_dir) == 0:
        seq_dir = 'sequences/'
    if seq_dir[-1] != '/':
        seq_dir += '/'
        
    alignedfiles = os.listdir(seq_dir)
    for alignedfile in alignedfiles:
        transcript = alignedfile[:-14]
        
        sequences.setdefault(transcript, {})
        contig = ''
        
        for line in open(seq_dir + alignedfile):
            if line[0] == '>':
                contig = line.strip()[1:]
                sequences[transcript].setdefault(contig, '')
            else:
                sequences[transcript][contig] += line.strip()
                                
    for datafile in expressionfiles:
        for line in open(datafile, 'r'):
            elm = line.strip().split('\t')
            transcript_id = ''
            parts = []
            if line[:7] == 'TRINITY':
                parts = elm[0].split(':')[0].split('_')
                transcript_id = '_'.join(parts[:-2])
            elif line[:2] == 'TR':
                division = elm[0].split('|')
                parts = division[1].split('_')
                transcript_id = division[0] + '_' + '_'.join(parts[:-2])
            elif line[:4] == 'comp':
                parts = elm[0].split('_')
                transcript_id = '_'.join(parts[:-2])            
                
            if len(parts) > 0:
                contig_id = parts[-2]
                region_id = transcript_id + '_' + str(parts[-1])
                ex_values.setdefault(region_id, {})
                ex_values[region_id].setdefault(contig_id, [])
                if file_format == 'kallisto':
                    ex_values[region_id][contig_id].append(float(elm[-1]))                    
                elif file_format == 'trinity':
                    ex_values[region_id][contig_id].append(float(elm[-2]))
                elif file_format == 'rsem':
                    ex_values[region_id][contig_id].append(float(elm[-3]))
                    
    for transcript in ex_values.keys():
        for contig in ex_values[transcript].keys():
            sum_fpkms = 0
            for i in range(len(expressionfiles)):
                sum_fpkms += log(ex_values[transcript][contig][i] + 1, 2)
            if sum_fpkms >= threshold:
                log_ex.setdefault(transcript, {})
                log_ex[transcript].setdefault(contig, [])
                for i in range(len(expressionfiles)):
                    log_ex[transcript][contig].append(log(ex_values[transcript][contig][i] + 1, 2))
    
        if transcript in log_ex.keys() and len(list(log_ex[transcript].keys())) <= 1:
            del(log_ex[transcript])
            
    return log_ex, sequences

def parseFilesWithDifferentSpecies (expressionfiles, file_format='kallisto', threshold=1, seq_dir=''):
    from math import log
    import sys, os
    
    if len(expressionfiles) <= 1:
        sys.stderr.write('Input more than one expression files.')
        sys.exit(-1)
    
    ex_values = {}
    log_ex = {}
    sequences = {}
    og_inverse = {}

    if len(seq_dir) == 0:
        seq_dir = 'orthofinder_summary/'
    if seq_dir[-1] != '/':
        seq_dir += '/'

    alignedfiles = os.listdir(seq_dir)
    for alignedfile in alignedfiles:
        transcript = alignedfile[:-14]
        
        sequences.setdefault(transcript, {})
        contig = ''
        
        for line in open(seq_dir + alignedfile):
            if line[0] == '>':
                contig = line.strip()[1:]
                sequences[transcript].setdefault(contig, '')
                og_inverse[contig] = transcript
            else:
                sequences[transcript][contig] += line.strip()
                
    for datafile in expressionfiles:
        for line in open(datafile, 'r'):
            elm = line.strip().split('\t')
            if elm[0] == 'target_id' or elm[0] not in og_inverse.keys():
                continue
            
            transcript_id = og_inverse[elm[0]]
            ex_values.setdefault(transcript_id, {})
            ex_values[transcript_id].setdefault(elm[0], [])
            if file_format == 'kallisto':
                ex_values[transcript_id][elm[0]].append(float(elm[-1]))                    
            elif file_format == 'trinity':
                ex_values[transcript_id][elm[0]].append(float(elm[-2]))
            elif file_format == 'RSEM':
                ex_values[transcript_id][elm[0]].append(float(elm[-3]))
                
    for transcript in ex_values.keys():
        for contig in ex_values[transcript].keys():
            sum_fpkms = 0
            for i in range(len(expressionfiles)):
                sum_fpkms += log(ex_values[transcript][contig][i] + 1, 2)
            if sum_fpkms >= threshold:
                log_ex.setdefault(transcript, {})
                log_ex[transcript].setdefault(contig, [])
                for i in range(len(expressionfiles)):
                    log_ex[transcript][contig].append(log(ex_values[transcript][contig][i] + 1, 2))
    
        if transcript in log_ex.keys() and len(list(log_ex[transcript].keys())) <= 1:
            del(log_ex[transcript])
            
    return log_ex, sequences

if __name__ == '__main__':
    directory = '/Users/yonezawa/research/data/chlorella/'
    expresionfiles = []
    for species in ['cva', 'cvu']:
        expresionfiles.append(directory + 'kallisto_' + species + '/abundance.tsv')
        
    log_ex, sequences = parseFilesWithDifferentSpecies(expresionfiles, file_format = 'kallisto', threshold = 1, seq_dir = directory + 'orthofinder_summary')
    
    for transcript in log_ex.keys():
        print('{0:s}\t{1:d}'.format(transcript, len(list(sequences[transcript].keys()))))