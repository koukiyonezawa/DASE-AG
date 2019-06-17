'''
Created on 2018/10/04

@author: yonezawa
'''

import re, codecs

def parseConfig(config_file):
    configs = {}
    for line in codecs.open(config_file, 'r', 'utf-8'):
        if line.find('=') > 0:
            elm = re.split('\s*=\s*', line.strip())
            
            item, value = elm[0], elm[1].split()[0]
            configs[item] = value
            
    return configs

if __name__ == '__main__':
    config_file = 'config.txt'
    configs = parseConfig(config_file)
    
    for item, value in configs.items():
        print(item + '\t' + value)