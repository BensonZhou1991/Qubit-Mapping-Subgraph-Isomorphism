# Created on June 10, 2019 by Sanjiang Li || mrlisj@gmail.com

import networkx as nx
from ag import q20 # architecture graph
from inimap import _tau_bsg_, _tau_bstg_ # two initial mappings
from QTokyo_f import qubit_in_circuit, entail, swap, topgates, topgates2, topgates3, unoccupied,\
     swap_reduce_min_dist, SWAP3x, SWAP4x, greedy_solved_gates, complete_mapping, good_next_mappingx, qct
import json
import os

# method for initial mapping
# 'weighted_graph', 'empty' ot 'topsubgraph'
method_IM = 'weighted_graph'

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# define the architecture graph
# IBM Q Tokyo (Q20)

def qubit_in_circuit(D): # the set of logic qubits appeared in D, a subset of C
    ''' Return the set of qubits in a circuit D

    Args:
        D (list): a sublist of CNOT gates of the input circuit C
    Returns:
        Q (set): the set of qubits in D
    '''
    Q = set()
    for gate in D:
        Q.add(gate[0])
        Q.add(gate[1])
    return Q

H = q20()
G = nx.Graph.to_undirected(H)
EG = nx.edges(G)

res_name = []
res_original_gates = []
res_output_gates = []
res_initial_map = []
res_time = []
results = {'names':res_name, 'inout gate number':res_original_gates,
           'output gate number':res_output_gates, 'initial map':res_initial_map,
           'time consumption':res_time}

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
### The input circuit C    
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# Note: we compare three initial mappings:
    # 1. weighted graph initial mapping, 2. empty initial mapping, 3. topgraph initial mapping


path = "CNOT_lists2/" #文件夹目录
files= os.listdir(path) #得到文件夹下的所有文件名称

count = 0
for file_name in files:
    count += 1
    current_path = 'CNOT_lists2/' + file_name
    print('\nThis is the %dth circuit'%count)
    print('The name of the circuit is', file_name)
    with open(current_path, 'r') as f:
        # 98 gates, 10 qubits 
        # results: +69 (8.5s, weighted graph) +57 (5.2s, empty) +54 (12.8s, topgraph)
        # weighted graph initial mapping 
        # [0, 1, 2, 20, 20, 3, 6, 4, 5, 20, 20, 7, 8, 9, 20, 20, 20, 20, 20, 20]
        # final mapping
        # [5, 6, 4, 2, 20, 1, 7, 8, 20, 20, 20, 0, 3, 9, 20, 20, 20, 20, 20, 20]
    
    ##with open('CNOT_lists//max46.txt', 'r') as f:
        # 11844 gates, 10 qubits
        # results +4308 (130.8s) +4560 (125.5s) +4410 (132.6s)
        # initial mapping
        # [20, 2, 3, 20, 20, 4, 7, 8, 0, 20, 5, 6, 9, 1, 20, 20, 20, 20, 20, 20]
        # final mapping
        # [20, 2, 5, 20, 20, 4, 7, 6, 0, 20, 20, 3, 8, 20, 20, 20, 1, 9, 20, 20]
    
    ##with open('CNOT_lists//co14_215.qasm.txt', 'r') as f:
        # 7840 gates， 15 qubits 
        # results: 4563 (194.6s) +4629 (212.6s) +3975 (196.4s) vs. 8982
        # initial mapping
        # [3, 11, 4, 14, 0, 10, 5, 12, 13, 1, 6, 9, 2, 20, 20, 7, 8, 20, 20, 20]
        # final mapping
        # [3, 2, 14, 10, 5, 8, 13, 12, 11, 1, 7, 9, 4, 0, 20, 20, 20, 6, 20, 20]
            
        sqn = json.loads(f.read())
        res_name.append(file_name)
    C = sqn
    l = len(C)
    res_original_gates.append(l)
    print('The input circuit contains %s gates' %l)
    
    nl = len(qubit_in_circuit(C))
        
    #\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
      ### select an initial mapping ### 
    #\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
    
    ##################################################################
    ### We compare three different initial mappings
    
    ### 1. weighted graph initial mapping
    if method_IM == 'weighted_graph':
        #print('We use the weighted graph initial mapping tau:')
        tau = _tau_bsg_(C,G)
    
    ### 2. empty mapping
    if method_IM == 'empty':
        tau = [20]*20
        #print('We use the empty initial mapping tau:')
    
    ### 3. topsubgraph mapping
    if method_IM == 'topsubgraph':
        #print('We use the topgraph initial mapping tau:')
        tau = _tau_bstg_(C,G,nl)
    
    res_initial_map.append(tau)
    ##################################################################
    
    C_out, cost_time = qct(C,G,EG,tau)
    res_output_gates.append(len(C_out))
    res_time.append(cost_time)
    
    print('The output circuit contains %s gates' %len(C_out))

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
