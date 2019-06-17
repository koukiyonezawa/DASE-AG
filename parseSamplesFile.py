'''
Created on 2018/10/06

@author: yonezawa
'''

import codecs

def parseSamples(samples_file):
    samples = {}
    
    for line in codecs.open(samples_file, 'r', 'utf-8'):
        elm = line.strip().split('\t')
        
        replicate = elm[1]
        parts = replicate.split('_')
        rep_name = elm[0]
        if len(parts) >= 2:
            rep_name = '_'.join(parts[:-1])
        samples.setdefault(rep_name, {})
        rep_number = 1
        while parts[-1][-rep_number].isnumeric():
            rep_number += 1
        rep_number -= 1
        if rep_number < 0:
            rep_number = 'rep0'
        else:
            rep_int = int(parts[-1][-rep_number])
            rep_number = 'rep{0:d}'.format(rep_int)
        
        samples[rep_name].setdefault(rep_number, [])
        for i in range(2, len(elm)):
            samples[rep_name][rep_number].append(elm[i])
    
    return samples


if __name__ == '__main__':
    samples_file = '../mouse_Rett/samples_file.txt'
    
    samples = parseSamples(samples_file)
    for rep_name in samples.keys():
        print(rep_name)
        for rep_number in samples[rep_name].keys():
            print(rep_number)
            print(', '.join(samples[rep_name][rep_number]))
