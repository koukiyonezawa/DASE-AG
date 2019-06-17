'''
Created on 2016/05/30

@author: yonezawa
'''

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
            

def calcAngles (expressions, mode='cosine'):
    from anomalyDetection import medianAndCovariance # @UnresolvedImport
    import sys
    angles = {}
    
    if mode == 'Mahalanobis':
        _, V = medianAndCovariance(expressions)
    
    if mode not in ['Mahalanobis', 'cosine', 'longer_arc']:
        sys.stderr.write('Distance mode should be either of Mahalanobis, cosine or longer_arc.')
        sys.exit(-1)
    
    for transcript in expressions.keys():
        contigs = list(expressions[transcript].keys())
        
        for i in range(len(contigs) - 1):
            for j in range(i + 1, len(contigs)):
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


if __name__ == '__main__':
    from parseFiles import parseFiles  # @UnresolvedImport
    
    directory = '/Users/yonezawa/research/data/human-PNN/'
    expresionfiles = []
    coverage_threshold = 0.8
    for ERR in ['kallisto_control1/', 'kallisto_knockdown1/']:
        expresionfiles.append(directory + ERR + 'abundance.tsv')
        
    log_ex, sequences = parseFiles(expresionfiles, file_format = 'kallisto', threshold = 2, seq_dir = directory + 'aligned_contigs')
