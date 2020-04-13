# Created on June 10, 2019 by Sanjiang Li || mrlisj@gmail.com
## designed for checking how a <= 20-bit circuit can be efficiently implemented in IBM Q Tokyo （Q20）
## the architecture graph can be replaced with any other undirected graph

import networkx as nx
from vfs import Vf
import json
#import matplotlib.pyplot as plt
import time

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
start = time.time()
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# define the architecture graph
# IBM Q Tokyo (Q20) 
def q20():
    g = nx.DiGraph()
    g.add_nodes_from([0,19])
    for i in range(0,4):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)

    for i in range(5,9):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(10,14):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(15,19):
        g.add_edge(i,i+1)
        g.add_edge(i+1,i)
        
    for i in range(0,15):
        g.add_edge(i,i+5)
        g.add_edge(i+5,i)

    for i in [1,3,5,7,11,13]:
        g.add_edge(i,i+6)
        g.add_edge(i+6,i)

    for i in [2,4,6,8,12,14]:
        g.add_edge(i,i+4)
        g.add_edge(i+4,i)
    return g


H = q20()
G = nx.Graph.to_undirected(H)

# G has 20 vertices and 43 edges


##for e in H.edges():
##    print(e)
##print(len(H.edges()))
##print(len(G.edges()))

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
### The input circuit C    
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# Note: we compare three initial mappings:
    # 1. weighted graph initial mapping, 2. empty initial mapping, 3. topgraph initial mapping

    
##with open('CNOT_lists//sys6.txt', 'r') as f:
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

with open('CNOT_lists//co14_215.qasm.txt', 'r') as f:
    # 7840 gates， 15 qubits 
    # results: 4563 (194.6s) +4629 (212.6s) +3975 (196.4s) vs. 8982
    # initial mapping
    # [3, 11, 4, 14, 0, 10, 5, 12, 13, 1, 6, 9, 2, 20, 20, 7, 8, 20, 20, 20]
    # final mapping
    # [3, 2, 14, 10, 5, 8, 13, 12, 11, 1, 7, 9, 4, 0, 20, 20, 20, 6, 20, 20]
        
    sqn = json.loads(f.read())    
C = sqn
l = len(C)
    
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
print(time.asctime())
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

### GRAPH OF THE CIRCUIT
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

def qubit_in_circuit(D): # the set of logic qubits appeared in D, a subset of C
    Q = set()
    for gate in D:
        Q.add(gate[0])
        Q.add(gate[1])
    return Q

QIC = qubit_in_circuit(C) # qubits in C
nl = len(QIC) # number of qubits in C

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
print('The input circuit C has %s CNOT gates and uses %s qubits in %s.' %(l,nl,QIC))
print('  We transform C into a physical circuit executatble in IBM Q Tokyo.')
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

# Note: different from the article, tau: V -> Q is a mapping from V to Q
        # V: the 20 physical qubits in Q20), Q: the set of logic qubits in the circuit C
        # tau is encoded as a list with length 20
        # tau[i] is the logical qubit associated with the i-th physical qubit in Q20
        # tau[i] = 20 if i is not assigned/occupied
    # gate = [0]*2 represents a CNOT gate, where gate[0] and gate[1] are qubits in C

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# check if tau satisfies/entails/solves/executes a gate
def entail(tau,gate):
    u = gate[0]
    v = gate[1]
    if (u not in tau) or (v not in tau):
        return False
    else:
        p = tau.index(gate[0]) # p is the physical qubit that gate[0] occupies 
        q = tau.index(gate[1]) # q is the physical qubit that gate[1] occupies 

        if (p,q) in nx.edges(G): # G is the architecture graph
            return True    
        else:
            return False

#transform tau to tau' with the images of p and q swapped
def swap(tau,p,q): # p, q are physical qubits
    if (p,q) not in nx.edges(G):
        return tau
    else:
        taunew = tau[:]
        taunew[p] = tau[q]
        taunew[q] = tau[p]
        return taunew       
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

# D is a subset of C, LD is the corresponding index set
def topgates(LD):
    # return the index set of topgates in the *first* front layer of D
    LT = []
    N = set() # the set of qubits of gates in LT
    for i in LD:
        gate = C[i]
        p = gate[0]
        q = gate[1]

        if len(N) >= nl-1:  # nl is the number of qubits in the circuit C
            return LT

        # if one qubit of the gate has appeared in N, it's not a topgate, and all edges after it are not topgates
        if p in N or q in N:
            N.add(p)
            N.add(q)
            continue
        else: 
            N.add(p)
            N.add(q)
            LT.append(i)
    return LT

# D is a subset of C, LD is the corresponding index set
# go to the second layer
def topgates2(LD):
    # return the index set of topgates in the *second* layer of D
    LDx = LD[:]
    
    LT = topgates(LDx)
    for i in LT:
        LDx.remove(i)
        
    LT2 = topgates(LDx)
    
    return LT2

# D is a subset of C, LD is the corresponding index set
# go to the third layer
def topgates3(LD):
    # return the index set of topgates in the *third* layer of D
    LDx = LD[:]
    
    LT = topgates(LDx)
    for i in LT:
        LDx.remove(i)
        
    LT = topgates(LDx)
    for i in LT:
        LDx.remove(i)
        
    LT3 = topgates(LDx)
    
    return LT3

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

##############################################################################################
        # Note: Fallback and mapping evloving process \\
# in case tau does not solve any gate in TG, the following function outputs a swap \
  # that decreases the minimum distance from tau to TG
  # we use this function as our fallback if SWAPx3 and SWAPx4 do not give a useful
      # good new mapping,
      # but this never happens for Q20 and all circuits we tested
  # in case the initial mapping is not complete, we also need this function to
      # extend the mapping step by step
##############################################################################################

# the list of logical qubits unoccupied by tau
def unoccupied(tau):
    U = list(p for p in range(20) if tau[p] == 20)
    return U

# tau is the current mapping, TG a subset of gates in C  
def swap_reduce_min_dist(tau,TG):
    md = 5 # the diameter of G is 4
    # Note: this parameter md is adjustable w.r.t. the architecture graph
    
    MD = []
    Q = qubit_in_circuit(TG)
    U = unoccupied(tau)
    path = []
    taux = tau[:]
    for gate in TG: # gate = [x,y] is a gate in TG, where x, y are logic qubits      
        u = gate[0]
        v = gate[1]
        
        if (u not in tau) and (v not in tau): #type 0
            for p in U:
                for q in U:
                    if p == q:
                         continue
                    d = nx.shortest_path_length(G,p,q)
                    if d >= md:
                         continue
                    elif d == 1:
                         taux[p] = u
                         taux[q] = v
                         return taux, [] # tau is extended
                    else:
                         md = d
                         MD = [[u,p], [v,q], 0] # 0 denotes the type

        elif (u in tau) and (v not in tau): #type 1
            p = tau.index(u)
            for q in U:
                d = nx.shortest_path_length(G,p,q)
                if d >= md:
                     continue
                elif d == 1:
                     taux[q] = v
                     return taux, [] # tau is extended
                else:
                     md = d
                     MD = [[u,p], [v,q], 1]
                         
        elif (u not in tau) and (v in tau): #type 2
            q = tau.index(v)
            for p in U:
                d = nx.shortest_path_length(G,p,q)
                if d >= md:
                     continue
                elif d == 1:
                     taux[p] = u
                     return taux, [] # tau is extended
                else:
                     md = d
                     MD = [[u,p], [v,q], 2]

        else:  #type 3
            p = tau.index(u) # p is a physical qubit
            q = tau.index(v) # q is a physical qubit
            d = nx.shortest_path_length(G,p,q)
            if d >= md:
                continue
            elif d == 1: # tau entails (u,v)
                 return tau, [] # tau is not changed as it can execute some gates in TG
            else:
                 md = d
                 MD = [[u,p], [v,q], 3]
                 
    if TG and not MD:
        print('MD is empty. This is impossible. Please check!')
##        print(taux)
##        print(list([i,C[i]] for i in TG))
        return taux, []
                            
    u = MD[0][0]                               
    p = MD[0][1]
    v = MD[1][0]
    q = MD[1][1]
    t = MD[2] # type
    path = []
    if t == 0:
        taux[p] = u
        taux[q] = v
    elif t == 1:
        taux[q] = v
    elif t == 2:
        taux[p] = u
    else: # t == 3
        pi = nx.shortest_path(G,p,q)
        if p != pi[0]:
            print('Something wrong, p:%s and pi[0]:%s should be the same' %(p,pi[0]))  
##        print('We select swap [%s, %s] to reduce the minimum distance!' %(pi[0],pi[1]))
        taux[p] = tau[pi[1]]
        taux[pi[1]] = tau[p]
        path = [[p,pi[1]]]
    return taux, path # update tau with swaps in path

#####################################################################
  ###\__/#      GreedyV3 SWAPS         \__/#\#/~\\__/#\__/#\#/~\__/
#####################################################################

# SWAP3 returns all swaps that lead to mappings close (dist≤3) to tau
# D is a subset of the input circuit, and LD its index set
def SWAP3x(tau,LD):
    # We replace Q with Q2 in the second and third levels

    LTG = topgates(LD)
    TG = list(C[i] for i in LTG)
    Q = qubit_in_circuit(TG)

    LTG2 = topgates2(LD)
    TG2 = list(C[i] for i in LTG2)
##    Q2 = qubit_in_circuit(TG+TG2) #x00
    # if we use x00 instead of x01, the quality could be further improved while very slow
    Q2 = qubit_in_circuit(TG2) #x01

    # added @ Aug 6, 2019
    LTG3 = topgates3(LD)
    TG3 = list(C[i] for i in LTG3)
##    Q3 = qubit_in_circuit(TG+TG2+TG3)
    Q3 = qubit_in_circuit(TG3)
    
    SWAP3 = []
    for edge_1 in nx.edges(G):
        p = edge_1[0]
        q = edge_1[1]
        
        if tau[p] not in Q and tau[q] not in Q: # both p, q are irrelevant phy. qubits w.r.t. TG
            continue

        tau1 = swap(tau,p,q)
        x1 = [ [edge_1], tau1]
        if x1 not in SWAP3:
            SWAP3.append(x1)

        for edge_2 in nx.edges(G):
            if edge_1 == edge_2:
                continue
            
            p = edge_2[0]
            q = edge_2[1]

            # modified @ July 11, 2019 Q -> Q2           
            if tau1[p] not in Q2 and tau1[q] not in Q2:
                continue
            
            tau2 = swap(tau1,p,q)           
            x2 = [ [edge_1,edge_2], tau2]
            if x2 not in SWAP3:
                SWAP3.append(x2)

            for  edge_3 in nx.edges(G):
                if edge_3 in {edge_1, edge_2}:
                    continue
                p = edge_3[0]
                q = edge_3[1]
                
                if tau2[p] not in Q2 and tau2[q] not in Q2: #x00
##                if tau2[p] not in Q3 and tau2[q] not in Q3: #x01
                    continue

                tau3 = swap(tau2,p,q)
                x3 = [ [edge_1, edge_2, edge_3], tau3]
                if x3 not in SWAP3:
                    SWAP3.append(x3)
                    
    return SWAP3

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# SWAP4 returns swaps that lead to mappings 4-close (dist=4) to tau
    # By far, no tested circuits go to this level
    # Consider delete this function
    
def SWAP4x(tau,LD):

    # revised @ July 11, 2019
    LTG2 = topgates2(LD)
    TG2 = list(C[i] for i in LTG2)
    Q2 = qubit_in_circuit(TG2)
    
    SWAP4 = []
    z3 = SWAP3x(tau,LD)
    
    for x in z3:
        tau3 = x[1]
        s = x[0]      # Achtung! 1 <= len(s) <= 3 
        for edge_4 in nx.edges(G):

            p = edge_4[0]
            q = edge_4[1]

            # modified @ July 11, 2019 Q -> Q2           
            if tau3[p] not in Q2 and tau3[q] not in Q2:
                continue
                
            tau4 = swap(tau3,p,q)
            s.append(edge_4)
            x4 = [s,tau4]
            if (x4 not in SWAP4):
                SWAP4.append(x4)

    return SWAP4

#####################################################################
###             Greedy V3 SWAPS
#####################################################################

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# How many gates tau can solve greedily?
# This subroutine iteratively computes those topgates in D (could be a solvable subgraph) that tau can solve

def greedy_solved_gates(tau,LD): 
    #tau is a mapping, D a part of the given input circuit C, LD the index set of D
        # Output LTG: the index set (a subset of L) of a maximal subset of *top*gates that can be solved by tau

    LDx = LD[:]
    LSTG = [] # the index set of all solved topgates (in perhaps two or more consecutive layers that are solvable by tau)
    while True:
        # as long as there are gates that can be solved by this mapping
        if not LDx:
            return LSTG
        
        LX = topgates(LDx) # topgates returns indices of topgates
        ldx = len(LDx)
        
        # examine gates in X one by one, check if they are entailed by tau
        for i in LX:
            gate = C[i]
            if entail(tau,gate):
                LDx.remove(i)
                LSTG.append(i) # C[i] is solved by tau                   
                                    
        if len(LDx) == ldx: # tau does not solve any new topgate
            return LSTG

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# check if tau is complete w.r.t. the input circuit C, where nl is the numnber of qubits in C
def complete_mapping(tau):
    if tau.count(20) == 20 - nl:
        return True
    else:
        return False
    
# Starting from tau, what is the next best mapping for gates in D?
# Before calling this function, usually we should remove all topgates that are executable by tau
    # that is, greedy_solved_gates(tau,LD) should be empty

def good_next_mappingx(tau,LD):
    # we assume that tau cannot solve any gates in D

    LTG = topgates(LD)
    TG = list(C[i] for i in LTG)
    Q = qubit_in_circuit(TG)
    if not complete_mapping(tau):
##        print('The mapping %s is incomplete %s' %(tau,complete_mapping(tau)))
        UU = swap_reduce_min_dist(tau,TG) # UU has form [taunew, path]
        if not UU[1]:
##            print('  We extend the incomplete mapping to an unoccupied physical qubit.')
            return UU[1], UU[0]
            
    # rho in e.g. SWAP3(tau) has form [[edge_1,...],taunew], where tau is a mapping, [edge_1,...] is a sequence of swaps
    SG = [] # elements in SG has form [rho, t], where rho in SWAP3x and t the efficiency rate
    for rho in SWAP3x(tau,LD):
        u = len(rho[0])
        x = greedy_solved_gates(rho[1],LD) # x has form [solved_gates]
        if not x:
            continue
        t = len(x)/(3*u) # the efficiency rate of rho is defined as (number of solved gates)/(3*(number of swaps used))
        SG.append([rho,t])

    if not SG:
        print('Yes, we go to the fourth level')
        for rho in SWAP4x(tau,LD):
            u = len(rho[0]) # u < 3 is possible
            x = greedy_solved_gates(rho[1],LD)
            if not x:
                continue
            t = len(x)/(3*u) # the efficiency rate of rho is defined as (number of solved gates)/(3*(number of swaps used))
            SG.append([rho,t])
            
    taunew = tau[:]
    
    if not SG: # if no proper mapping in SWAP1-4, take a random neighbor of tau
        UU = swap_reduce_min_dist(tau,TG) # UU has form [taunew, path]
        print('  Oops, there are no good swaps, we select one swap that reduces the minimal distance')
        return UU[1], UU[0]
    
    # else, we select the mapping with the best efficiency rate from SG
    SG.sort(key=lambda t:t[1], reverse=True)
    sg = SG[0] # with form [rho, val], where rho has form [ [edge1,...], taunew]
    path = sg[0][0] # with form [edge1,...], where each edge1 is a swap [p,q]
    taunew = sg[0][1]
##    if len(path) == 3:
##        print('Yes, we go to the third level')
                          
    return path, taunew              

##################################################################
  #\__/#\__/#\#/~\ THE MAIN GREEDYV3 ALGORITH10M \__/#\__/#\#/~\
##################################################################
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# Weighted SUBGRAPH
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

# each circuit C induces an undirected graph
def graph_of_circuit(C):
    g = nx.Graph()
    g.add_nodes_from(QIC)
    for gate in C:
        g.add_edge(gate[0],gate[1])
    return g

def best_g(C):
    g_of_c = graph_of_circuit(C)
    Eg = list(g_of_c.edges())
    Eg.sort()

    edge_wgt_list = list([C.count([e[0],e[1]]) + C.count([e[1],e[0]]), e] for e in g_of_c.edges())
    edge_wgt_list.sort(key=lambda t: t[0], reverse=True) # q[0] weight, q[1] edge
    edge_wgt_list1 = list(q[1] for q in edge_wgt_list)

    Egw = list(edge_wgt_list1)
    legw = len(Egw)
    
    # We use the vf2 algorithm to determine if a graph is embeddable in G and compute a mapping if it does
    vf2 = Vf()
    result = {} 
    result = vf2.dfsMatch(g_of_c, G, result)
    lng = len(nx.nodes(g_of_c))
    if len(result) == lng:
##        print('The graph of the circuit is embeddable in G')
        return g_of_c

##    print('The graph of the circuit is NOT embeddable in G')

    # we search backward, remove the first edge that makes g not embeddable, and continue till we get a maximum graph
    Hard_Edge_index = 0 # the index of the first hard edge
    g = nx.Graph()
    # we don't include all nodes from QIC as sometimes my vf2 alg. is perhaps not good when g is disconnected

    # add the first edge into g
    q = Egw[0]
    g.add_edge(q[0],q[1])
    
    rp = 0
    for rp in range(legw):

        # h is the index of the last edge that can be added into g
        h = Hard_Edge_index
        if h == legw - 1: 
            return g
        
        LW = Egw[h+1:legw]

        for q in LW:           
            g.add_edge(q[0],q[1])           
            i = Egw.index(q)
            
            # find the largest i such that the first i-1 edges are embeddable
            vf2 = Vf()
            result = {} 
            result = vf2.dfsMatch(g, G, result)
            if len(result) != len(nx.nodes(g)):                    
                Hard_Edge_index = i
                g.remove_edge(q[0],q[1])
                break
            
            if i == legw-1 and len(result) == len(nx.nodes(g)):                
                return g
    return g

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# the weighted subgraph initial mapping
def _tau_bsg_(C):
    bsg = best_g(C)
    vf2 = Vf()
    result = {}
    result = vf2.dfsMatch(bsg, G, result) 
    BB = result
    tau = [20]*20
    for key in BB: # map physical qubit BB[key] to logic qubit key
        tau[BB[key]] = key
    return tau

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
##################################################################
# compute the gate dependency graph of (a part of) the circuit
def gate_dependency_graph(LX): # LX is the index set of a subset of C
    g = nx.DiGraph()
    g.add_nodes_from(LX)
    L1 = LX[:]
    while L1:
        L0 = L1[:]
        TG = topgates(L0)
        for i in TG:
            L0.remove(i)
        TGx = topgates(L0)
        for i in TG:
            gate1 = C[i]
            for j in TGx:
                gate2 = C[j]
                if set(gate1) & set(gate2):
                     g.add_edge(i,j)
                     
        L1 = L0[:]
    return g

# compute the top subgraph of C
def _topgraph_(C,L1): 
    # C is the input circuit
    # L1 is the index set of unsolved gates in C
    gdg = gate_dependency_graph(L1)

    # construct the top subgraph g from the circuit 
    g = nx.Graph() # the subgraph of G induced by nodes with indices in GATES
    GATES = [] # the indices of topgates that can be executed (i.e., put in g) in this round
    
    # we consider unsolved gate in C one by one
    Dump = set() # record all those gates that cannot be put in the solvable graph in this round
    for j in L1:
        a = L1.index(j)
        x = C[j]

        CHANGE = True        
        while CHANGE and Dump:
            ldump = len(Dump)
            # for every s in Dump, add its descedents into s
            Dump_TG = topgates(list(Dump))
            for s in Dump_TG:
                Dump = Dump | set(nx.descendants(gdg, s))

            if ldump == len(Dump):
                CHANGE = False
                
        if Dump >= set(L1[a:]):
            return g, GATES

        if j in Dump:
            continue
                        
        if (x[0],x[1]) in g.edges(): # the same gate has been considered for some j<i with C[i]==C[j]
            GATES.append(j) # C[i] is to be solved in this round
            continue
        
        # Temporalliy add x into g
        g.add_node(x[0])
        g.add_node(x[1])
        g.add_edge(x[0],x[1])
        
        #############################################################
        #### Decide if g remains solvable after we added x into g ###
        #############################################################

        vf2 = Vf()
        result = {}
        result = vf2.dfsMatch(g, G, result)

        if len(result) == len(g.nodes()): # g is embeddable
            GATES.append(j) # C[j] will be solved in this round

        # if g is not embeddable, then remove x from g
        else:               
            # remove same gates from consideration
            Dump.add(j)
                    
            g.remove_edge(x[0],x[1])
            if nx.degree(g,x[0]) == 0:
                g.remove_node(x[0])
            if nx.degree(g,x[1]) == 0:
                g.remove_node(x[1])

    return g, GATES

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# the topgraph initial mapping
def _tau_bstg_(C):
    l = len(C)
    L = list(range(l))
    g_GATES =_topgraph_(C,L)
    bstg = g_GATES[0] # the topgraph of C
    vf2 = Vf()
    result = {}
    result = vf2.dfsMatch(bstg, G, result)
    
    BB = result
    tau = [20]*20    
    for key in BB: # map physical qubit BB[key] to logic qubit key
        tau[BB[key]] = key
    return tau

#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
  ### select an initial mapping ### 
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

##################################################################
### We compare three different initial mappings

### 1. weighted graph initial mapping
##print('We use the weighted graph initial mapping tau:')
##tau = _tau_bsg_(C)

### 2. empty mapping
##tau = [20]*20
##print('We use the empty initial mapping tau:')

### 3. topsubgraph mapping
print('We use the topgraph initial mapping tau:')
tau = _tau_bstg_(C)
##################################################################

print(tau)
print(' Here tau maps a physical qubit in G to a logic qubit and the i-th physical qubit is free if tau[i] is 20') 

L = list(range(l)) # the index set of C
nsvg = 0 # the number of executed gates
COST = 0 # the number of auxiliary CNOT gates
taunew = tau[:]
L1 = L[:]
np = 0
C_out = [] # the output circuit
while nsvg < l:

    # first check what gates can be executed by taunew (obtained in the previous round),
     # then remove their indices from L1 
    GSG = greedy_solved_gates(taunew,L1)       
    for i in GSG:
        # if C[i] is executed, remove i from L1 and add CNOT gate [p,q] to C_out
        L1.remove(i)
        p = taunew.index(C[i][0])
        q = taunew.index(C[i][1])
        C_out.append([p,q])
                
    nsvg += len(GSG) # the number of solved gates

##    print('Round %s solved %s gates' %(np, len(GSG)))
##    print('Round %s solved gates in %s' %(np, list(C[i] for i in GSG)))

    if not L1:
##        print('We have solved all gates within %s rounds' %(np+1))
        break
    
    np += 1
    Y = good_next_mappingx(taunew,L1)
    path = Y[0] # a sequence of swaps that transforms the previous mapping into the mapping tuanew below
    taunew = Y[1]
    cost = len(path)*3
##    print('Round %s: The next mapping' %np)
##    print('   %s' %taunew)
##    print('   is obtained by swaps in %s are enforced with cost %s' %(path,cost))
        
    COST += cost
    for edge in path:
        # implement in C_out the swap corresponding to edge
        C_out.append([edge[0],edge[1]])
        C_out.append([edge[1],edge[0]])
        C_out.append([edge[0],edge[1]])
        
print('We need to add %s auxiliary CNOT gates' %COST)
print(l,nsvg) # check if l=nsvg
print('The final mapping is %s' %taunew)
##print('The output circuit contains %s gates' %len(C_out))   
   
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
end = time.time()
print('It takes time: %s seconds.' %(end-start))
#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
