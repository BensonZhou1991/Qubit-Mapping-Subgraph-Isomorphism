# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 10:56:11 2019

@author: zxz58
"""

'''
this module is to check the functionality of transformed circuit
'''

def SwapMap(current_map, swap1, swap2):
    temp = current_map[swap1]
    current_map[swap1] = current_map[swap2]
    current_map[swap2] = temp
    
def CheckTypeOfCX(C_in_remain, C_out_remain, current_map):
    current_CX = C_out_remain[0]
    length = len(C_in_remain)
    current_control_phy = current_CX[0]
    current_target_phy = current_CX[1]
    current_control = current_map[current_control_phy]
    current_target = current_map[current_target_phy]
    pos = None
    for i in range(length):
        if C_in_remain[i][0] == current_control or C_in_remain[i][1] == current_control or \
           C_in_remain[i][0] == current_target or C_in_remain[i][1] == current_target:
               pos = i
               break
    if pos == None: return 'wrong', None
    if C_in_remain[pos][0] == current_control and C_in_remain[pos][1] == current_target:
        return 'CX', pos
    else:
        if C_out_remain[1][0] == current_target_phy and C_out_remain[1][1] == current_control_phy and\
           C_out_remain[2][0] == current_control_phy and C_out_remain[2][1] == current_target_phy:
               return 'swap', current_control_phy, current_target_phy
        else:
               return 'wrong', None

def CheckFunctionality(C_in, C_out, initial_map):
    current_map = initial_map.copy()
    C_in_remain = C_in.copy()
    C_out_remain = C_out.copy()
    
    while len(C_out_remain) != 0:
# =============================================================================
#         print('output C index is ', len(C_out)-len(C_out_remain))
#         print('input C index is ', len(C_in)-len(C_in_remain))
#         print('mapping is ', current_map)
# =============================================================================
        type_CX = CheckTypeOfCX(C_in_remain, C_out_remain, current_map)
        if type_CX[0] == 'CX':
            C_out_remain.pop(0)
            C_in_remain.pop(type_CX[1])
        if type_CX[0] == 'swap':
            C_out_remain.pop(0)
            C_out_remain.pop(0)
            C_out_remain.pop(0)
            SwapMap(current_map, type_CX[1], type_CX[2])
        if type_CX[0] == 'wrong':
            return False
    return True
        