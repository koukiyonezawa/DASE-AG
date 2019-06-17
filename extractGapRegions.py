'''
Created on 2018/07/31

@author: yonezawa
'''

def extractGapAndPeripheralRegions(read_length, gene, seqs, filtered_contigs, gap_threshold=8, edge_trim=True):
    first_contig = filtered_contigs[gene][0]
    seqs_without_scratches = {}

    for i in range(len(seqs[first_contig])):
        gap_flag = True
        for contig in filtered_contigs[gene]:
            if seqs[contig][i] != '-':
                gap_flag = False
                
        if not gap_flag:
            for contig in filtered_contigs[gene]:
                seqs_without_scratches.setdefault(contig, '')
                seqs_without_scratches[contig] += seqs[contig][i]

    gap_ranges = []
    gap_start = 0
    gap_end = 0
    for contig in seqs_without_scratches.keys():
        for i in range(len(seqs_without_scratches[contig])):                
            if seqs_without_scratches[contig][i] == '-':
                if i == 0 or seqs_without_scratches[contig][i - 1] != '-':
                    gap_start = i
                elif i == len(seqs_without_scratches[contig]) - 1:
                    gap_end = i
                    gap_ranges.append([gap_start, gap_end])
            elif i > 0 and seqs_without_scratches[contig][i - 1] == '-':
                gap_end = i - 1
                gap_ranges.append([gap_start, gap_end])
    
    gap_read_ranges = []
    remove_list = []
    for gap_range in gap_ranges:
        gap_start = gap_range[0]
        gap_end = gap_range[1]
        gap_read_range = [max(0, gap_start - read_length), min(len(seqs_without_scratches[first_contig]) - 1, gap_end + read_length)]        
        
        intersection = False
        for i in range(len(gap_read_ranges)):
            if gap_read_range[0] <= gap_read_ranges[i][0] and gap_read_range[1] >= gap_read_ranges[i][0] and gap_read_range[1] < gap_read_ranges[i][1]:
                gap_read_ranges[i][0] = gap_read_range[0]
                intersection = True
            elif gap_read_range[0] > gap_read_ranges[i][0] and gap_read_range[0] <= gap_read_ranges[i][1] and gap_read_range[1] >= gap_read_ranges[i][1]:
                gap_read_ranges[i][1] = gap_read_range[1]
                intersection = True
            elif gap_read_range[0] <= gap_read_ranges[i][0] and gap_read_range[1] >= gap_read_ranges[i][1]:
                gap_read_ranges[i] = gap_read_range
                intersection = True
            elif gap_read_range[0] >= gap_read_ranges[i][0] and gap_read_range[1] <= gap_read_ranges[i][1]:
                intersection = True
            elif gap_read_range[0] == gap_read_ranges[i][0] and gap_read_range[1] >= gap_read_ranges[i][1]:
                gap_read_ranges[i][1] = gap_read_range[1]
                intersection = True
            elif gap_read_range[0] == gap_read_ranges[i][0] and gap_read_range[1] < gap_read_ranges[i][1]:
                intersection = True
            elif gap_read_range[1] == gap_read_ranges[i][1] and gap_read_range[0] <= gap_read_ranges[i][0]:
                gap_read_ranges[i][0] = gap_read_range[0]
                intersection = True
            elif gap_read_range[1] == gap_read_ranges[i][1] and gap_read_range[0] > gap_read_ranges[i][0]:
                intersection = True
                
            if edge_trim and (gap_read_ranges[i][0] == 0 or gap_read_ranges[i][1] == len(seqs_without_scratches[first_contig]) - 1):
                remove_list.append(gap_read_ranges[i])
        
        if not intersection and gap_range[1] - gap_range[0] >= gap_threshold and (not edge_trim or (gap_read_range[0] > 0 and gap_read_range[1] < len(seqs_without_scratches[first_contig]) - 1)):
            gap_read_ranges.append(gap_read_range)
    
    result_ranges = []
    for cand in gap_read_ranges:
        if cand not in remove_list:
            check_flag = True
            
            for contig in seqs_without_scratches.keys():
                start = cand[0]
                end = cand[1]
                if seqs_without_scratches[contig][start:end + 1] == '-' * (end - start + 1):
                    check_flag = False
            
            if check_flag:
                result_ranges.append(cand)
        
    return result_ranges, seqs_without_scratches
                    

if __name__ == '__main__':
    import os
    
    directory = os.getcwd()
    if directory[-1] != '/':
        directory += '/'
        
    trimmed_dir = 'filtered_sequences/'
    
    if not os.path.exists(directory + trimmed_dir):
        os.mkdir(directory + trimmed_dir)
        
    alignment_dir = directory + 'aligned_contigs/'
    filtered_contigs = {}
    
    read_length = 125
    filtered_genes = 0
    excluded_genes = 0
    
    for line in open(directory + 'filtered_contigs.dat'):
        elm = line.strip().split('\t')
        filtered_contigs[elm[0]] = elm[1].split(',')
        
    range_file = open(directory + 'ranges.dat', 'w')
    allseq_file = open(directory + 'focused_ranges.fasta', 'w')
    for alignment in os.listdir(alignment_dir):
        gene = alignment[:-14]
        
        if gene in filtered_contigs.keys():
            seqs = {}
            flag = False
            for line in open(alignment_dir + alignment):
                if line[0] == '>':
                    contig = line.strip()[1:]
                    
                    if contig in filtered_contigs[gene]:
                        flag = True
                        seqs.setdefault(contig, '')
                    else:
                        flag = False
                        
                elif flag:
                    seqs[contig] += line.strip()
                    
            gap_read_ranges, seqs_without_scratches = extractGapAndPeripheralRegions(read_length, gene, seqs, filtered_contigs)
            
            if len(gap_read_ranges) > 0:
                writefile = open(directory + trimmed_dir + gene + '.fasta', 'w')
                count = 1
                for gap_read_range in sorted(gap_read_ranges, key=lambda x: x[0]):
                    range_file.write('{0:s}_{1:d}:{2:d}-{3:d}\n'.format(gene, count, gap_read_range[0], gap_read_range[1]))
                    for contig in filtered_contigs[gene]:
                        writefile.write('>{0:s}_{1:s}_{2:d}\n'.format(gene, contig, count))
                        allseq_file.write('>{0:s}_{1:s}_{2:d}\n'.format(gene, contig, count))
                        for i in range(len(seqs_without_scratches[contig])):
                            if seqs_without_scratches[contig][i] != '-':
                                writefile.write(seqs_without_scratches[contig][i])
                                if i >= gap_read_range[0] and i <= gap_read_range[1]:
                                    allseq_file.write(seqs_without_scratches[contig][i])
                        writefile.write('\n')
                        allseq_file.write('\n')
                    count += 1
                writefile.close()
                filtered_genes += 1
            else:
                excluded_genes += 1

    range_file.close()
    allseq_file.close()
    print('{0:d} genes were filtered.'.format(filtered_genes))
    print('{0:d} genes were excluded.'.format(excluded_genes))