'''
Created on 2017/02/13z

@author: yonezawa
'''

def oneClassSVM (data, nu = 0.1, kernel = 'rbf'):
    import numpy as np
    from sklearn.svm import OneClassSVM
    import matplotlib.pyplot as plt
    
    original_data = []    
    labels = []
    max_expression = 0
    
    for transcript in data.keys():
        for variant in sorted(data[transcript].keys()):
            labels.append(transcript + '_' + variant)
            
            this_data = data[transcript][variant]
            tmp_data = []
            samples = list(this_data.keys())
            for sample in samples:
                tmp_data.append(this_data[sample])
                if max_expression < this_data[sample]:
                    max_expression = this_data[sample]
                    
            original_data.append(tmp_data)

    data_input = np.array(original_data)
    
    clf = OneClassSVM(nu = nu, kernel = kernel)
    _ = clf.fit(data_input)
    prediction = clf.predict(data_input)
    
    x, y = np.meshgrid(np.linspace(0, max_expression, 100), np.linspace(0, max_expression, 100))
    Z = clf.decision_function(np.c_[x.ravel(), y.ravel()]).reshape(x.shape)
    max_expression = 15
    s = 40
    plt.contourf(x, y, Z, levels = np.linspace(0, Z.max(), 7))
    normal = data_input[prediction == 1]
    anomaly = data_input[prediction == -1]
    """
    for i in range(prediction.shape[0]):
        if prediction[i] == 1:
            normal.append(original_data[i])
        else:
            anomaly.append(original_data[i])
    """
    normal = np.array(normal)
    anomaly = np.array(anomaly)
    plt.scatter(normal[:, 0], normal[:, 1], c = 'blueviolet', s = s)
    plt.scatter(anomaly[:, 0], anomaly[:, 1], c = 'gold', s = s)
    plt.axis('tight')
    plt.xlim((0, max_expression))
    plt.ylim((0, max_expression))
    plt.show()
    
    return prediction, labels, Z

def mad (vector):
    from numpy import median
    tmp = vector - median(vector)
    
    for i in range(tmp.size):
        if tmp[i] < 0:
            tmp[i] *= -1
    
    return median(tmp)

def MSDestimator (data, iteration = 10):
    import numpy as np
    # import matplotlib.pyplot as plt
    from numpy.random import randn
    
    original_data = []    
    labels = []
    max_expression = 0
    
    for transcript in data.keys():
        for variant in sorted(data[transcript].keys()):
            labels.append(transcript + '_' + variant)
            
            this_data = data[transcript][variant]
            tmp_data = []
            samples = list(this_data.keys())
            for sample in samples:
                tmp_data.append(this_data[sample])
                if max_expression < this_data[sample]:
                    max_expression = this_data[sample]
                    
            original_data.append(tmp_data)

    data_input = np.array(original_data)
    
    # center the data using an L1 estimator
    center = np.median(data_input, axis = 0)
    centered_data = np.zeros(data_input.shape)
    
    for i in range(data_input.shape[0]):
        centered_data[i] = data_input[i] - center
        
    # initialize a weighting vector
    delta = np.ones(data_input.shape[0])
    
    # main procedure
    v = randn(centered_data.shape[1], centered_data.shape[1])
    for k in range(iteration):
        # first vector of v
        v0 = v[0]
        v0 /= np.linalg.norm(v0)
        v[0] = v0
        
        for i in range(1, v.shape[0]):
            vi = v[i]
            # Gram-Schmidt orthonormalization
            for j in range(i):
                vi -= np.dot(v[i], v[j]) * v[j]
            v[i] = vi / np.linalg.norm(vi)
            
        # calculate delta_k
        delta_k = np.ones(delta.size)
        
        for i in range(centered_data.shape[0]):
            tmp = np.zeros(centered_data.shape[0])
            for j in range(v.shape[0]):
                for l in range(centered_data.shape[0]):
                    tmp[l] = np.dot(v[j], centered_data[l])
                
                r_tmp = abs(np.dot(v[j], centered_data[i]) - np.median(tmp)) / mad(tmp)
                
                if r_tmp <= 4 and r_tmp > 2.5:
                    delta_k[i] *= 2.5 / r_tmp
                elif r_tmp > 4:
                    delta_k[i] = 0
                
            # update delta
            if delta_k[i] < delta[i]:
                delta[i] = delta_k[i]
        
        print('Procedure {0:d} finished.'.format(k + 1))

    # location vector
    u = np.zeros(centered_data.shape[1])
    for i in range(centered_data.shape[0]):
        u += delta[i] * centered_data[i]
    u /= np.sum(delta)
    u += center
    
    # scale covariance matrix
    V = np.zeros((centered_data.shape[1], centered_data.shape[1]))
    denominator = 0
    for i in range(delta.size):
        denominator += delta[i] ** 2
        V_tmp = delta[i] ** 2 * np.dot(data_input[i] - u, np.transpose(data_input[i] - u))
        V += V_tmp
    V /= denominator
        
    """
    normal = data_input[prediction == 1]
    anomaly = data_input[prediction == -1]
    
    xx, yy = np.meshgrid(np.linspace(0, max_expression, 100), np.linspace(0, max_expression, 100))
    _ = clf._decision_function(np.c_[xx.ravel(), yy.ravel()]).reshape(xx.shape())
                                                                      
    s = 40
    normal = np.array(normal)
    anomaly = np.array(anomaly)
    plt.scatter(normal[:, 0], normal[:, 1], c = 'blueviolet', s = s)
    plt.scatter(anomaly[:, 0], anomaly[:, 1], c = 'gold', s = s)
    plt.axis('tight')
    plt.xlim((0, max_expression))
    plt.ylim((0, max_expression))
    plt.show()
    """

    return u, V

def medianAndCovariance (data):
    import numpy as np
    original_data = []    
    labels = []
    max_expression = 0
    
    for transcript in data.keys():
        for variant in sorted(data[transcript].keys()):
            labels.append(transcript + '_' + variant)
            
            this_data = data[transcript][variant]
            tmp_data = []
            for i in range(len(this_data)):
                tmp_data.append(this_data[i])
                if max_expression < this_data[i]:
                    max_expression = this_data[i]
                    
            original_data.append(tmp_data)

    data_input = np.array(original_data)
    
    matrix = np.zeros((data_input.shape[1], data_input.shape[1]))
    
    median = np.median(data_input, axis = 0)
    
    for i in range(data_input.shape[1]):
        var_sum = 0
        for j in range(data_input.shape[0]):
            var_sum += (data_input[j, i] - median[i]) ** 2
                    
        matrix[i, i] = var_sum / data_input.shape[1]
        
        var_sum = 0
        for j in range(i + 1, data_input.shape[1]):
            for k in range(data_input.shape[0]):
                var_sum += (data_input[k, i] - median[i]) * (data_input[k, j] - median[j])
            
            matrix[i, j] = var_sum / data_input.shape[1]
            matrix[j, i] = matrix[i, j]
    
    return median, matrix

def Mahalanobis (p, u, V):
    import numpy as np
    
    return np.dot(np.dot(np.transpose(p - u), np.linalg.inv(V)), p - u) ** 0.5
    
if __name__ == '__main__':
    from math import log2
    import numpy as np, matplotlib.pyplot as plt, matplotlib
    
    directory = '/Users/yonezawa/research/data/paralvinella/'
    expressions = {}
    tissues = ['hot', 'cold']
    threshold = 2
    nu_value = {}
    iteration = 5
    
    datafile = directory + 'Tri_tra_Phe.TMM.fpkm.matrix'
        
    for line in open(datafile):
        elm = line.strip().split('\t')
        
        if elm[0] != 'result_hot':
            if log2(float(elm[1]) + 1) + log2(float(elm[2]) + 1) >= threshold:
                parts = elm[0].split('_')
                transcript = '_'.join(parts[:-1])
                variant = parts[-1]
                expressions.setdefault(transcript, {})
                expressions[transcript].setdefault(variant, [])            
                for i in range(1, 3):
                    tissue = tissues[2 - i]
                    expressions[transcript][variant].append(log2(float(elm[i]) + 1))
                    
    u, V = medianAndCovariance(expressions)
    print(u)
    distances = {}
    max_dist = 0
    for transcript in expressions.keys():
        for variant in expressions[transcript].keys():
            contig = transcript + '_' + variant
            this_data = []
            for tissue in expressions[transcript][variant].keys():
                this_data.append(expressions[transcript][variant][tissue])
                
            this_data = np.array(this_data)
            
            distances[contig] = Mahalanobis(this_data, u, V)
            if max_dist < distances[contig]:
                max_dist = distances[contig]
    
    for transcript in expressions.keys():
        for variant in expressions[transcript].keys():        
            contig = transcript + '_' + variant
            x, y = expressions[transcript][variant][0], expressions[transcript][variant][1]
            plt.scatter(x, y, c = distances[contig], cmap = matplotlib.cm.get_cmap('bwr'), vmax = max_dist, vmin = 0)
    
    plt.colorbar()
    plt.xlim(0, 20)
    plt.ylim(0, 20)
    plt.show()
    
    """
    for nu in nus:
        prediction, labels, _ = oneClassSVM(expressions, nu = nu, kernel = 'rbf')
        
        for i in range(len(labels)):
            if prediction[i] > 0:
                parts = labels[i].split('_')
                transcript = '_'.join(parts[:-1])
                variant = parts[-1]
                
                this_data = []
                for tissue in tissues:
                    this_data.append(expressions[transcript][variant][tissue])
        
        for i in range(len(labels)):
            if prediction[i] < 0:
                parts = labels[i].split('_')
                transcript = '_'.join(parts[:-1])
                variant = parts[-1]
                
                nu_value.setdefault(labels[i], nu)
                if nu_value[labels[i]] > nu:
                    nu_value[labels[i]] = nu
        
        print('nu: {0:.4f}'.format(nu))
        break
    """
    """    
    for label in sorted(nu_value.keys(), key = lambda x: nu_value[x]):
        parts = label.split('_')
        transcript = '_'.join(parts[:-1])
        variant = parts[-1]
        print('{0:s}\t{1:.4f}\t{2:.4f}\t{3:.4f}'.format(label, nu_value[label], expressions[transcript][variant]['hot'], expressions[transcript][variant]['cold']))
    """