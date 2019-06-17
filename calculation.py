'''
Created on 2016/05/30

@author: yonezawa
'''

import networkx as nx

def calcNorm (p1, p2):
    result = 0
    
    for i in range(len(p1)):
        result += p1[i] * p2[i]
    
    return result

def calcInnerProduct (p1, p2):
    return calcNorm(p1, p2) / (calcNorm(p1, p1) * calcNorm(p2, p2)) ** 0.5

def calcCosine (v1, v2):
    from math import acos, pi
    inner_product = calcInnerProduct(v1, v2)
    if inner_product > 1:
        inner_product = 1
        
    return acos(inner_product) / pi * 180

def calcEuclid (v1, v2):
    squared_sum = 0
    
    for i in range(len(v1)):
        squared_sum += (v1[i] - v2[i]) ** 2
    
    return squared_sum ** 0.5

def calcManhattan (v1, v2):
    absolute_sum = 0
    
    for i in range(len(v1)):
        absolute_sum += abs(v1[i] - v2[i])
    
    return absolute_sum

def calcMahalanobis (v1, v2, V):
    import numpy as np, sys
    from scipy.special import ellipeinc
    
    rotation_matrix = np.eye(v1.size)
    
    if V.shape[0] == 3:        
        for i in range(2):
            for _ in range(i + 1, v1.size):
                tmp_matrix = np.eye(v1.size)
                
                rotation_matrix = np.dot(rotation_matrix, tmp_matrix)

        v1 = v1.T.dot(rotation_matrix)
        v2 = v2.T.dot(rotation_matrix)
        V = V.dot(rotation_matrix)

    eigenvalues, eigenvectors = np.linalg.eig(V)
    v1new = np.transpose(np.dot(np.transpose(eigenvectors), np.transpose(v1)))
    v2new = np.transpose(np.dot(np.transpose(eigenvectors), np.transpose(v2)))
    theta = 0
    
    if V.shape[0] <= 3:
        longer_axis_length = eigenvalues[0]
        shorter_axis_length = eigenvalues[1]
        longer_axis = 0
        if eigenvalues[0] < eigenvalues[1]:
            longer_axis_length, shorter_axis_length = shorter_axis_length, longer_axis_length
            longer_axis = 1
            
        k = (1 - (shorter_axis_length / longer_axis_length) ** 2) ** 0.5
        longer_axis_vector = np.zeros(V.shape[0])
        longer_axis_vector[longer_axis] = 1
        
        arg1 = calcCosine(v1new, longer_axis_vector)
        arg2 = calcCosine(v2new, longer_axis_vector)
        
        if arg1 * arg2 >= 0:
            theta = np.abs(longer_axis_length * (ellipeinc(arg1 / longer_axis_length, k) - ellipeinc(arg2 / longer_axis_length, k)))
        else:
            theta = longer_axis_length * (np.abs(ellipeinc(arg1 / longer_axis_length, k)) + np.abs(ellipeinc(arg2 / longer_axis_length, k)))
        
    else:
        sys.stderr.write('The number of conditions should be less than 4 when the distance mode is \'Mahalanobis\'.')
        sys.exit(-1)
    
    return theta
            

def calcAngles (expressions, coverages, mode = 'Mahalanobis'):
    from anomalyDetection import medianAndCovariance # @UnresolvedImport
    import sys
    angles = {}
    
    if mode == 'Mahalanobis':
        _, V = medianAndCovariance(expressions)
    
    if mode not in ['Mahalanobis', 'cosine', 'longer_arc']:
        sys.stderr.write('Distance mode should be either of Mahalanobis, cosine or longer_arc.')
        sys.exit(-1)
    
    for transcript in expressions.keys():
        if transcript in coverages.keys():
            contigs = list(expressions[transcript].keys())
            
            for i in range(len(contigs) - 1):
                for j in range(i + 1, len(contigs)):
                    if contigs[i] in coverages[transcript].keys() and contigs[j] in coverages[transcript][contigs[i]].keys():
                        angles.setdefault(transcript, {})
                        angles[transcript].setdefault(contigs[i], {})
                        if mode == 'cosine':
                            angles[transcript][contigs[i]][contigs[j]] = calcCosine(expressions[transcript][contigs[i]], expressions[transcript][contigs[j]])
                        elif mode == 'longer_arc':
                            angles[transcript][contigs[i]][contigs[j]] = calcCosine(expressions[transcript][contigs[i]], expressions[transcript][contigs[j]]) * max(calcNorm(expressions[transcript][contigs[i]], expressions[transcript][contigs[i]]), calcNorm(expressions[transcript][contigs[j]], expressions[transcript][contigs[j]])) ** 0.5                        
                        elif mode == 'Mahalanobis':
                            angles[transcript][contigs[i]][contigs[j]] = calcMahalanobis(expressions[transcript][contigs[i]], expressions[transcript][contigs[j]], V)    
                        angles[transcript].setdefault(contigs[j], {})
                        angles[transcript][contigs[j]][contigs[i]] = angles[transcript][contigs[i]][contigs[j]]
                    
    return angles

def calcCoverages (sequences):
    coverages = {}
    for transcript in sequences.keys():
        contig_list = list(sequences[transcript].keys())
        if len(contig_list) > 1:
            coverages.setdefault(transcript, {})
            for i in range(len(contig_list) - 1):
                seq1 = sequences[transcript][contig_list[i]]
                for j in range(i + 1, len(contig_list)):
                    seq2 = sequences[transcript][contig_list[j]]
                    total_length = 0
                    same_site = 0
                    
                    for k in range(len(seq1)):
                        if seq1[k] != '-' or seq2[k] != '-':
                            total_length += 1
                            if seq1[k] == seq2[k]:
                                same_site += 1
                    
                    coverages[transcript].setdefault(contig_list[i], {})
                    coverages[transcript][contig_list[i]][contig_list[j]] = float(same_site) / total_length
                    coverages[transcript].setdefault(contig_list[j], {})
                    coverages[transcript][contig_list[j]][contig_list[i]] = float(same_site) / total_length
    return coverages

def sum_edge_weights(nodes, G):
    sum_weights = 0
    
    for edge in G.edges(data=True):
        start = edge[0]
        end = edge[1]
        weight = edge[2]['weight']
        if start in nodes and end in nodes:
            sum_weights += weight

    return sum_weights

def create_graph(coverage, coverage_threshold, gene):
    G = nx.Graph()
    
    transcripts = list(coverage[gene].keys())
    for i in range(len(transcripts) - 1):
        t1 = transcripts[i]
        for j in range(i + 1, len(transcripts)):
            t2 = transcripts[j]
            cov = coverage[gene][t1][t2]
            
            if cov >= coverage_threshold:
                G.add_nodes_from([t1, t2])
                G.add_edge(t1, t2, weight=cov)
                    
    return G

def find_max_clique(G):
    cliques = nx.find_cliques(G)
    max_clique = None
    for c in cliques:
        if max_clique is None or len(max_clique) < len(c):
            max_clique = c
        elif len(max_clique) == len(c):
            max_clique_total_weight = sum_edge_weights(max_clique, G)
            c_total_weight = sum_edge_weights(c, G)
                
            if max_clique_total_weight < c_total_weight:
                max_clique = c
    
    return max_clique

def trim_transcripts(coverages, coverage_threshold):
    remove_gene_list = []
    remaining_transcripts = {}
    for gene in coverages.keys():
        G = create_graph(coverages, coverage_threshold, gene)
        max_clique = find_max_clique(G)
        
        if max_clique is None:
            remove_gene_list.append(gene)
        else:
            remaining_transcripts.setdefault(gene, [])
            for transcript in max_clique:
                remaining_transcripts[gene].append(transcript)
                
    for gene in remove_gene_list:
        del(coverages[gene])
    
    for gene in remaining_transcripts.keys():
        remove_transcripts = []
        transcript_list = list(coverages[gene].keys())
        for t in transcript_list:
            if t not in remaining_transcripts[gene]:
                remove_transcripts.append(t)
        
        for t in transcript_list:
            for rt in remove_transcripts:
                if t in coverages[gene].keys() and rt in coverages[gene][t].keys():
                    del(coverages[gene][t][rt])
                    
        for rt in remove_transcripts:
            if rt in coverages[gene].keys():
                del(coverages[gene][rt])
                
    return coverages

if __name__ == '__main__':
    from parseFiles import parseFiles  # @UnresolvedImport
    
    directory = '/Users/yonezawa/research/data/human-PNN/'
    expresionfiles = []
    coverage_threshold = 0.8
    for ERR in ['kallisto_control1/', 'kallisto_knockdown1/']:
        expresionfiles.append(directory + ERR + 'abundance.tsv')
        
    log_ex, sequences = parseFiles(expresionfiles, file_format = 'kallisto', threshold = 2, seq_dir = directory + 'aligned_contigs')

    coverages = calcCoverages(sequences)
    angles = calcAngles(log_ex, coverages, mode='cosine')
    
    coverages = trim_transcripts(coverages, coverage_threshold)
        
    writefile = open(directory + 'PNN_comparison.dat', 'w')
    
    for transcript in angles.keys():
        contig_list = list(sorted(angles[transcript].keys(), key = lambda x: int(x[1:])))
        
        if len(contig_list) > 2:
            print(transcript)
            for i in range(len(contig_list) - 1):
                c1 = contig_list[i]
                for j in range(i + 1, len(contig_list)):
                    c2 = contig_list[j]
                    
                    print('{0:s}\t{1:s}\t{2:.3f}\t{3:.3f}'.format(c1, c2, coverages[transcript][c1][c2], angles[transcript][c1][c2]))
                    writefile.write('{0:s}\t{1:s}\t{2:.3f}\t{3:.3f}\n'.format(c1, c2, coverages[transcript][c1][c2], angles[transcript][c1][c2]))
    writefile.close()