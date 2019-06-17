'''
Created on 2018/08/03

@author: yonezawa
'''

def calcDistances(angles):
    distances = {}
    
    for transcript in angles.keys():
        contigs = sorted(list(angles[transcript].keys()), key = lambda x: int(x[1:]))
        distance = 0
        
        for i in range(len(contigs) - 1):
            for j in range(i + 1, len(contigs)):
                distance = angles[transcript][contigs[i]][contigs[j]]
                distances.setdefault(transcript, {})
                distances[transcript].setdefault(contigs[i], {})
                distances[transcript][contigs[i]][contigs[j]] = distance
                distances[transcript].setdefault(contigs[j], {})
                distances[transcript][contigs[j]][contigs[i]] = distance
    
    return distances

def rankRegions(distances, min_transcripts):
    from heap import heap, heappop  # @UnresolvedImport
    distance_list = []
    contig_list = []
    
    for transcript in distances.keys():
        contigs = list(distances[transcript].keys())
        
        if len(contigs) >= min_transcripts:        
            for i in range(len(contigs)):
                sum_dist = 0
                count = 0
                
                for j in range(len(contigs)):
                    if i != j and contigs[j] in distances[transcript][contigs[i]].keys():
                        sum_dist += distances[transcript][contigs[i]][contigs[j]]
                        count += 1
                
                parts = transcript.split('_')
                variant_id = '_'.join(parts[:-1])
                region_id = parts[-1]
                contig_list.append(variant_id + '_' + contigs[i] + '_' + region_id)
                distance_list.append(sum_dist / count)
    
    heap_distances, heap_contigs = heap(distance_list, contig_list)
        
    return heappop(heap_distances, heap_contigs)

def filterRegions(distances, max_rank, min_transcripts=2):
    ordered_list = rankRegions(distances, min_transcripts)
    
    return ordered_list[:max_rank + 1]

if __name__ == '__main__':
    from parseFiles import parseFiles # @UnresolvedImport
    from parseExonFiles import parseExonFiles  # @UnresolvedImport
    from parseSequenceFiles import parseSequences  # @UnresolvedImport
    from parseConfigFile import parseConfig  # @UnresolvedImport
    from parseSamplesFile import parseSamples  # @UnresolvedImport
    from parseReadFile import calcLength  # @UnresolvedImport
    from calculation import calcCoverages, trim_transcripts  # @UnresolvedImport
    from calculateExpressions import executeKallisto, executeRSEM # @UnresolvedImport
    from calculationForRegions import calcAngles  # @UnresolvedImport
    from createMap import geneTransMap # @UnresolvedImport
    from divideSequences import divideSequences  # @UnresolvedImport
    from executeMafft import executeMafft  # @UnresolvedImport
    from executeTrinity import executeTrinity  # @UnresolvedImport
    from selectRepresentativeSequences import summarizeRepresentative  # @UnresolvedImport
    from extractGapRegions import extractGapAndPeripheralRegions # @UnresolvedImport
    from enumerateCombination import combination # @UnresolvedImport
    import argparse, os, shutil, re
    
    directory = os.getcwd()
    usage = 'python {0:s}'.format(__file__)
    
    if directory[-1] != '/':
        directory += '/'
    
    parser = argparse.ArgumentParser(usage = usage)
    parser.add_argument('-rep', '--representative_sequences', type=str, dest='representative_file', help='fasta-format output file for the representative sequence of each transcripts. (default: representative_sequences.fasta)')
    parser.add_argument('-cf', '--config_file', type=str, dest='config', help='config file.', required=True)
    parser.add_argument('-o', '--output', type=str, dest='output_file', help='output file. (default: dase_{min_contigs}_{expression_threshold}_{coverage_threshold}_{dist_mode}.dat)')
    parser.add_argument('-nc', '--number_of_contigs', type=int, dest='min_transcripts', help='min. transcripts in one transcript considered. (default: 2)')
    parser.add_argument('-dist', '--distance_mode', type=str, dest='dist_mode', help='distance calculation. (default: cosine)')
    args = parser.parse_args()
    
    config_file = args.config
    
    min_transcripts = 2
    if args.min_transcripts:
        min_transcripts = args.min_contigs
        
    representative_file = directory + 'representative_sequences.fasta'
    if args.representative_file:
        representative_file = args.representative_file
    
    file_format = 'rsem'
        
    dist_mode = 'cosine'
    if args.dist_mode:
        dist_mode = args.dist_mode

    configs = parseConfig(config_file)
    print('Loaded the config file {0:s}.'.format(config_file))
    
    mafft_exe = configs['mafft_path']
    trinity_path = configs['trinity_path']
    if trinity_path[-1] != '/':
        trinity_path += '/'
    samples_file = configs['samples_file']
    bowtie2_dir = configs['bowtie2_path']
    cpu = int(configs['cpu'])
    max_memory = configs['max_memory']
    trimmomatic_flag = True
    if 'trimmomatic_flag' in configs.keys():
        trimmomatic_flag = configs['trimmomatic_flag']
    threshold = float(configs['expression_threshold'])
    coverage_threshold = float(configs['coverage_threshold'])
    gap_penalty = float(configs['gap_penalty'])

    resultfile = directory + 'dase_' + str(min_transcripts) + '_' + str(threshold) + '_' + str(coverage_threshold) + '_' + dist_mode + '.dat'
    if args.output_file:
        resultfile = directory + args.output_file
    filtered_contig_file = directory + 'filtered_contigs.dat'
    
    samples = parseSamples(samples_file)
    print('Loaded the samples file {0:s}.'.format(samples_file))
    
    executeTrinity(trinity_path, samples_file, max_memory, cpu, trimmomatic_flag)
    print('Executed Trinity.')

    seqfile = 'trinity_out_dir/Trinity.fasta'
    
    divideSequences(seqfile, seqdir=directory)
    print("Loaded the sequence file {0:s}.".format(seqfile))
    executeMafft(mafft_exe, directory=directory, gap_penalty=gap_penalty)
    print("Executed MAFFT for each transcript containing more than {0:d} variants.".format(min_transcripts - 1))        
    
    representatives = summarizeRepresentative(directory + 'aligned_contigs/')
    rep_file = open(representative_file, 'w')
    for transcript in sorted(representatives.keys()):
        isoform = list(representatives[transcript].keys())[0]
        rep_file.write('>' + transcript + '_' + isoform + '\n')
        rep_file.write(representatives[transcript][isoform] + '\n')
    rep_file.close()
    print("Selected representative sequences.")
    
    samples_keys = list(samples.keys())
    samples0_keys = list(samples_keys())
    rep = samples_keys[0] + '_rep00'
    if 'rep00' not in samples0_keys:
        rep = samples_keys[0] + '_' + samples0_keys[0]
    
    read_length = calcLength('trinity_out_dir/' + samples[samples_keys[0]][samples0_keys[0]][0] + '.PwU.qtrim.fq')

    trimmed_samples = samples
    if trimmomatic_flag:
        for sample in samples.keys():
            for rep in samples[sample].keys():
                for j in range(len(trimmed_samples[sample][rep])):
                    trimmed_samples[sample][rep][j] = 'trinity_out_dir/' + trimmed_samples[sample][rep][j] + '.P.qtrim.gz'
    
    sequences = parseSequences(directory + 'aligned_contigs/')
    coverages = calcCoverages(sequences)
    coverages = trim_transcripts(coverages, coverage_threshold)
    
    logfile = open(filtered_contig_file, 'w')
    for gene in coverages.keys():
        transcripts = list(coverages[gene].keys())
        if len(transcripts) >= 2:
            logfile.write(gene + '\t' + ','.join(sorted(transcripts)) + '\n')
    logfile.close()
    print("Recorded contigs to be taken into account to {0:s}.".format(filtered_contig_file))

    genes = list(coverages.keys())
    print('{0:d} genes extracted.'.format(len(genes)))
    number_of_transcripts = 0
    for gene in genes:
        number_of_transcripts += len(coverages[gene])
    print('{0:d} transcripts included.'.format(number_of_transcripts))        

    alignment_dir = directory + 'aligned_contigs/'
    filtered_contigs = {}
    trimmed_dir = directory + 'filtered_sequences/'
    if not os.path.exists(trimmed_dir):
        os.mkdir(trimmed_dir)
    range_file = open(directory + 'ranges.dat', 'w')
    allseq_file = open(directory + 'focused_ranges.fasta', 'w')

    for line in open(directory + 'filtered_contigs.dat'):
        elm = line.strip().split('\t')
        filtered_contigs[elm[0]] = elm[1].split(',')

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
                writefile = open(trimmed_dir + gene + '.fasta', 'w')
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

    geneTransMap(directory + 'focused_ranges.fasta', directory + 'focused_ranges.fasta.gene_trans_map')
    
    for sample in trimmed_samples.keys():
        for rep in trimmed_samples[sample].keys():
            for i in range(len(trimmed_samples[sample][rep])):
                if os.path.exists('rsem_{0:s}_{1:s}'.format(sample, rep)):
                    shutil.rmtree('rsem_{0:s}_{1:s}'.format(sample, rep))
                    print('Deleted rsem_{0:s}_{1:s}'.format(sample, rep))
            
    trinity_util_dir = trinity_path + 'util/'
    executeRSEM(trinity_util_dir, bowtie2_dir, directory + 'focused_ranges.fasta', trimmed_samples, cpu)
    
    expressionfiles = {}
    for cond in trimmed_samples.keys():
        expressionfiles.setdefault(cond, {})
        for rep in trimmed_samples[cond].keys():
            expressionfiles[cond][rep] = 'rsem_{0:s}_{1:s}.dat'
                
    total_angles = {}
    cond_list = list(expressionfiles.keys())
    all_rep_list = []
    comb_numbers = 1
    for cond in cond_list:
        rep_list = []
        for rep in expressionfiles[cond].keys():
            rep_list.append('rsem_{0:s}_{1:s}.dat'.format(cond, rep))
        all_rep_list.append(rep_list)
        comb_numbers *= len(rep_list)
        
    all_combinations = combination(all_rep_list, [], [], 0)    
    average_angles = {}
    for comb in all_combinations:
        log_ex, sequences = parseExonFiles(comb, file_format='rsem', threshold=threshold)
        angles = calcAngles(log_ex, dist_mode)
        for transcript in angles.keys():
            average_angles.setdefault(transcript, {})
            for c1 in average_angles[transcript].keys():
                average_angles[transcript].setdefault(c1, {})
                for c2 in average_angles[transcript][c1].keys():
                    average_angles[transcript][c1].setdefault(c2, 0)
                    average_angles[transcript][c1][c2] += angles[transcript][c1][c2] / comb_numbers
    
    distances = calcDistances(angles)
    print("Calculated the distances.")
    
    ordered_list = filterRegions(distances, len(distances) - 1, min_transcripts)
    print("Ranked gap regions according to the above distances.")
    
    matrixfile = 'abundance_matrix.isoform.counts.matrix'
    raw_expressions = {}
    
    for line in open(matrixfile):
        elm = re.split('\s+', line.strip())
        
        if len(elm[0]) > 10 and elm[0].split('_')[0] == 'TRINITY':
            gap_id = elm[0]
            raw_expressions.setdefault(gap_id, [])
            
            for i in range(1, len(elm)):
                raw_expressions[gap_id].append(float(elm[i]))
    
    writefile = open(resultfile, 'w')
    for pair in ordered_list:
        ag_region = pair[0]
        writefile.write('{0:s}\t{1:.4f}'.format(ag_region, pair[1]))
        
        for i in range(len(raw_expressions[ag_region])):
            writefile.write('\t{0:.2f}'.format(raw_expressions[ag_region][i]))
            
        writefile.write('\n')
    print('All tasks finished.')