'''
Created on 2019/01/04

@author: yonezawa
'''
import copy

def combination(input_list, appending_list, resulting_list, index):
    combination_number = 1
    for i in range(len(input_list)):
        combination_number *= len(input_list[i])
        
    if index == len(input_list) - 1:
        for i in range(len(input_list[index])):
            if i == 0:
                appending_list.append(input_list[index][i])
            else:
                appending_list[-1] = input_list[index][i]
            adding_list = copy.deepcopy(appending_list)
            resulting_list.append(adding_list)
        if len(resulting_list) == combination_number:
            return resulting_list
    else:
        for i in range(len(input_list[index])):
            if i == 0:
                appending_list.append(input_list[index][i])
            else:
                appending_list[-1] = input_list[index][i]
            adding_list = copy.deepcopy(appending_list)
            answer = combination(input_list, adding_list, resulting_list, index + 1)
            if answer is not None and len(answer) == combination_number:
                return answer

if __name__ == '__main__':
    input_list = [['a', 'b'], ['a', 'b', 'c'], ['a', 'b'], ['a', 'b', 'c', 'd']]
    print(combination(input_list, [], [], 0))