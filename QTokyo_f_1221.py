# Created on June 10, 2019 by Sanjiang Li || mrlisj@gmail.com
# Modified on December 21, 2019 to enhance Fallback mechanism
# Modified on April 19, 2020 // introduced QFilter_type

# This module defines the main algorithm
#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# Note: 
# G is the architecture graph
# EG is the list of edges in G
# tau a mapping from V to Q (!! this is different from the article, where tau : Q -> V)
        # V: the 20 physical qubits in Q20), Q: the set of logic qubits in the circuit C
        # tau is encoded as a list with length 20
        # tau[i] is the logical qubit associated with the i-th physical qubit in Q20
        # tau[i] = 20 if i is not assigned/occupied
# gate = [0]*2 represents a CNOT gate, where gate[0] and gate[1] are qubits in C
# nl is the number of qubits in Q
# selection: #x00 (current) #x01 (second) #x10 (third)
# todo: replace 20 with nq: number of qubits in G

#####\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

import networkx as nx
import time

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

# check if tau satisfies/entails/solves/executes a gate
def entail(EG,tau,gate):
    ''' Check if a mapping entails a gate

    Args:
        EG (list): the edge list of the architecture graph G
        tau (list): a mapping from V to Q (! this is different from the article, where tau: Q -> V), where
            V: the 20 physical qubits in Q20), Q: the set of logic qubits in the input circuit C
            tau is encoded as a list with length 20
            tau[i] is the logical qubit associated with the i-th physical qubit in Q20
            tau[i] = 20 if i is not assigned/occupied
        gate (list): [p,q] represents a CNOT gate, where p=gate[0] and q=gate[1] are qubits in C
    Returns:
        Boolean: True if tau entails gate
    '''

    u = gate[0]
    v = gate[1]
    if (u not in tau) or (v not in tau):
        return False
    else:
        p = tau.index(gate[0]) # p is the physical qubit that gate[0] occupies 
        q = tau.index(gate[1]) # q is the physical qubit that gate[1] occupies 

        if (p,q) in EG: # G is the architecture graph, EG its edge list
            return True    
        else:
            return False

#transform tau to tau' with the images of p and q swapped
def swap(EG,tau,p,q): # p, q are physical qubits
    ''' Swap the images of p and q in a mapping tau: V->Q

    Args:
        EG (list): the edge list of the architecture graph G
        tau (list): a mapping from V to Q (! this is different from the article, where tau: Q -> V), where
            V: the 20 physical qubits in Q20), Q: the set of logic qubits in the input circuit C
            p,q: two physical qubits in V
    Returns:
        taunew (list): a new mapping with the images of p and q swapped
    '''

    if (p,q) not in EG:
        return tau
    else:
        taunew = tau[:]
        taunew[p] = tau[q]
        taunew[q] = tau[p]
        return taunew       
###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

# D is a subset of C, LD is the corresponding index set
def topgates(C,nl,LD):
    ''' Return the index list of CNOT gates in the toplayer (front layer) of the current circuit D, where
        D is a sublist of the input circuit C, LD is the corresponding index list of D

    Args:
        C (list): the input circuit
        nl (int): the number of qubits in C
        LD (list): the index list of D, which represents the current logical circuit
    Returns:
        LT (list): the index list of CNOT gates which are in the front layer of D
    '''

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
# the second layer of gates in LD
def topgates2(C,nl,LD):
    ''' Return the index list of CNOT gates in the 2nd layer (the 1st look-ahead layer) of the current circuit D, where
        D is a sublist of the input circuit C, LD is the corresponding index list of D

    Args:
        C (list): the input circuit
        nl (int): the number of qubits in C
        LD (list): the index list of D, which represents the current logical circuit
    Returns:
        LT2 (list): the index list of CNOT gates which are in the 2nd layer (the 1st look-ahead layer) of D
    '''

    # return the index set of topgates in the *second* layer of D
    LDx = LD[:]
    
    LT = topgates(C,nl,LDx)
    for i in LT:
        LDx.remove(i)
        
    LT2 = topgates(C,nl,LDx)
    
    return LT2

# D is a subset of C, LD is the corresponding index set
# go to the third layer
def topgates3(C,nl,LD):
    ''' Return the index list of CNOT gates in the 3rd layer (the 2nd look-ahead layer) of the current circuit D, where
        D is a sublist of the input circuit C, LD is the corresponding index list of D

    Args:
        C (list): the input circuit
        nl (int): the number of qubits in C
        LD (list): the index list of D, which represents the current logical circuit
    Returns:
        LT3 (list): the index list of CNOT gates which are in the 3nd layer (the 2nd look-ahead layer) of D
    '''
    # return the index set of topgates in the *third* layer of D
    LDx = LD[:]
    
    LT = topgates(C,nl,LDx)
    for i in LT:
        LDx.remove(i)
        
    LT = topgates(C,nl,LDx)
    for i in LT:
        LDx.remove(i)
        
    LT3 = topgates(C,nl,LDx)
    
    return LT3

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
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
    ''' Return the list of physical qubits that are not mapped by tau

    Args:
        tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
    Returns:
        U (list): the list of physical qubits that are not mapped by tau (i.e., all p with tau[p]==20)
    '''
    U = list(p for p in range(20) if tau[p] == 20)
    return U

# tau is the current mapping, TG a subset of gates in C  
def swap_reduce_min_dist(G,tau,TG):
    ''' Fallback when SWAP3x and SWAP4x are empty;
        Generate a mapping from scratch when the emtpy initial mapping is used;
        Extend a mapping when tau is not complete
        
        Method: examine all gates ion TG one by one and record the last gate with minimum distance md in the form
            - md: the initial value of md (short for minimum_distance) depends on G and could be set as md = nx.diameter(G) + 1
            - MD: [[u,p], [v,q], type], where gate=[u,v], p,q are either unmapped or u==tau[p], v==tau[q], md==nx.shortest_path_length(G,p,q)            
        whenver nx.shortest_path_length(G,p,q) is 1, we return the updated taux and an empty action;
        otherwise, we compute a shortest path pi from p to q and swap the value of tau[p] and tau[p'] with p' the 2nd qubit (i.e., pi[1]) on pi 

    Args:
        G (undirected graph): the architecture graph 
        tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
        TG (list): a sublist of the input circuits, here denotes the front layer
    Returns:
        taux (list): an extended mapping from V to Q
        action (list): list of swaps (one or zero) that transform tau into the new mapping         
    '''
    md = 5 # the diameter of G (IBM Q20) is 4
    # md = nx.diameter(G) + 1
    # Note: this parameter md is adjustable w.r.t. the architecture graph
    # The function can be revised so that we extend tau only if we have no other choice
    
    MD = []
    Q = qubit_in_circuit(TG)
    U = unoccupied(tau)
    action = []
    taux = tau[:]
    for gate in TG: # gate = [u,v] is a gate in TG, where u, v are logic qubits      
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
                         md = d # temporalily store this best distance and the corresponding extension p->u and q-> v in MD
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
                 MD = [[u,p], [v,q], 3] #
                 
    if TG and not MD:
        print('MD is empty. This is impossible. Please check!')
        return taux, []

    #print(TG,MD)                            
    u = MD[0][0]                               
    p = MD[0][1]
    v = MD[1][0]
    q = MD[1][1]
    t = MD[2] # type
    action = []
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
        action = [[p,pi[1]]]
    return taux, action # update tau with swaps in action

#####################################################################
  ###\__/#      GreedyV3 SWAPS         \__/#\#/~\\__/#\__/#\#/~\__/
#####################################################################
#@April 19, 2020
def Q_Filter(C,EG,nl,tau,LD):
    ''' Output the four possible Qubit sets used as Q-filters
        Args:
            C (list): the input circuit
            EG (list): the list of edges of the architecture graph G
            nl (int): the number of qubits in C
            tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
            LD (list): the index list of D, which represents the current logical circuit           
        Returns:
            Q0: qubits in C
            Q1: qubits in the front layer (In the paper, this is known as Q0-filter)
            Q2: qubits in the second layer (In the paper, this is Q1-filter)
            Q2x: Q1+Q2 (In the paper, this is weak Q1-filter)
        '''
    Q0 = qubit_in_circuit(C)
    LTG = topgates(C,nl,LD)
    TG = list(C[i] for i in LTG)
    Q1 = qubit_in_circuit(TG)
    LTG2 = topgates2(C,nl,LD)
    TG2 = list(C[i] for i in LTG2)
    Q2 = qubit_in_circuit(TG2) #x00
    Q2x = qubit_in_circuit(TG+TG2) #x01
    return Q0, Q1, Q2, Q2x


# SWAP3x returns all swap sequences that lead to mappings close (dist≤3) to tau
# D is a subset of the input circuit, and LD its index set
def SWAP3x(C,EG,nl,tau,LD, QFliter_type):
    ''' Return all (Q-filtered) actions (swap sequences) with length ≤3, Q0 and Q2 filters are used to exclude possibly less relevant actions
        Args:
            C (list): the input circuit
            EG (list): the list of edges of the architecture graph G
            nl (int): the number of qubits in C
            tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
            LD (list): the index list of D, which represents the current logical circuit
            QFliter_type (str): '0' (no filter), '1' (Q0-filter as in the paper), '12' (Q0+Q1-filter), '12x' (Q0+weak_Q1-filter), '2x' (weak_Q1-filter)
        Returns:
            SWAP3 (list): each item has form [action, mapping], where action is a list consists of at most 3 swaps,
                        mapping is obtained by applying the action to tau
    '''
    QFilter = Q_Filter(C,EG,nl,tau,LD)
    if QFliter_type == '0': #no filter is used
        QF1 = QFilter[0]
        QF2 = QFilter[0]
        QF3 = QFilter[0]
    elif QFliter_type == '1': #only Q1-filter (note in the paper this is called Q0-filter)
        QF1 = QFilter[1]
        QF2 = QFilter[1]
        QF3 = QFilter[1]
    elif QFliter_type == '12': #Q1-filter for the first layer and Q2-filter for the rest (Q0+Q1-filter in the paper)
        QF1 = QFilter[1]
        QF2 = QFilter[2]
        QF3 = QFilter[2]
    elif QFliter_type == '12x': #Q1-filter for the first layer and (Q1+Q2)-filter for the rest (Q0+weak_Q1-filter
        QF1 = QFilter[1]
        QF2 = QFilter[3]
        QF3 = QFilter[3]
    elif QFliter_type == '2x': #only Q2x-filter
        QF1 = QFilter[3]
        QF2 = QFilter[3]
        QF3 = QFilter[3]
    else:
        pass
        
    SWAP3 = []
    for edge_1 in EG:
        p = edge_1[0]
        q = edge_1[1]
        if tau[p] not in QF1 and tau[q] not in QF1: # both p, q are irrelevant phy. qubits w.r.t. TG
            continue

        tau1 = swap(EG,tau,p,q)
        x1 = [ [edge_1], tau1]
        if x1 not in SWAP3:
            SWAP3.append(x1)

        for edge_2 in EG:
            if edge_1 == edge_2:
                continue
            
            p = edge_2[0]
            q = edge_2[1]
            if tau1[p] not in QF2 and tau1[q] not in QF2:
                continue
            
            tau2 = swap(EG,tau1,p,q)           
            x2 = [ [edge_1,edge_2], tau2]
            if x2 not in SWAP3:
                SWAP3.append(x2)

            for  edge_3 in EG:
                if edge_3 in {edge_1, edge_2}:
                    continue
                p = edge_3[0]
                q = edge_3[1]
                
                if tau2[p] not in QF3 and tau2[q] not in QF3:
                    continue

                tau3 = swap(EG,tau2,p,q)
                x3 = [ [edge_1, edge_2, edge_3], tau3]
                if x3 not in SWAP3:
                    SWAP3.append(x3)                    
    return SWAP3

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# SWAP4 returns swaps that lead to mappings 4-close (dist=4) to tau
    # By far, no tested circuits go to this level
    # Consider delete this function
    
def SWAP4x(C,EG,nl,tau,LD, QFilter_type):

    ''' Return all (Q-filtered) actions (swap sequences) with length ≤4 by considering one more swap than SWAP3x
        Args:
            C (list): the input circuit
            EG (list): the list of edges of the architecture graph G
            nl (int): the number of qubits in C
            tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
            LD (list): the index list of D, which represents the current logical circuit
            QFliter_type (str): '0' (no filter), '1' (Q0-filter as in the paper), '12' (Q0+Q1-filter),\
                                    '12x' (Q0+weak_Q1-filter), '2x' (weak_Q1-filter)
            
        Returns:
            SWAP4 (list): each item has form [action, mapping], where action is a list consists of at most 4 swaps,
                        mapping is obtained by applying the action to tau
    '''
    
    QFilter = Q_Filter(C,EG,nl,tau,LD)
    if QFliter_type == '0':
        QF4 = QFilter[0]
    elif QFliter_type == '1':
        QF4 = QFilter[1]
    elif QFliter_type == '12':
        QF4 = QFilter[2]
    elif QFliter_type == '12x':
        QF4 = QFilter[3]
    elif QFliter_type == '2x':
        QF4 = QFilter[3]
    else:
        pass
    
    SWAP4 = []
    z3 = SWAP3x(C,EG,nl,tau,LD,QFilter_type)
    
    for x in z3:
        tau3 = x[1]
        s = x[0]      # Achtung! 1 <= len(s) <= 3 
        for edge_4 in EG:

            p = edge_4[0]
            q = edge_4[1]

            if tau3[p] not in QF4 and tau3[q] not in QF4:
                continue
                
            tau4 = swap(EG,tau3,p,q)
            s.append(edge_4)
            x4 = [s,tau4]
            if (x4 not in SWAP4):
                SWAP4.append(x4)

    return SWAP4

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# This subroutine iteratively identifies those topgates in D (could be a solvable subgraph) that tau can solve/execute

def greedy_solved_gates(C,G,EG,nl,tau,LD):

    ''' Return the list of the indices of all gates solvable by a mapping tau

    Args:
        C (list): the input circuit
        G (graph): the architecture graph
        EG (list): the list of edges of G
        nl (int): the number of qubits in C
        tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
        LD (list): the index list of D, which represents the current logical circuit
        
    Returns:
        LSTG (list): the list of indices of solvable gates in D
    '''
    
    #tau is a mapping, D a part of the given input circuit C, LD the index set of D
        # Output LTG: the index set (a subset of L) of a maximal subset of *top*gates that can be solved by tau
    LDx = LD[:]
    LSTG = [] # the index set of all solved topgates (in perhaps two or more consecutive layers that are solvable by tau)

    while True:
        # as long as there are gates that can be solved by this mapping
        if not LDx:
            return LSTG
        
        LX = topgates(C,nl,LDx) # topgates returns indices of topgates
        ldx = len(LDx)
        
        # examine gates in LX one by one, check if they are entailed by tau
        for i in LX:
            gate = C[i]
            if entail(EG,tau,gate):
                LDx.remove(i)
                LSTG.append(i) # C[i] is solved by tau                   
                                    
        if len(LDx) == ldx: # tau does not solve any new topgate
            return LSTG

###\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
# check if tau is complete w.r.t. the input circuit C, where nl is the numnber of qubits in C
def complete_mapping(C,nl,tau):
    ''' Check if the mapping tau is complete

    Args:
        C (list): the input circuit
        nl (int): the number of qubits in C
        tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits        
    Returns:
        True if tau is complete False else
    '''
    
    if tau.count(20) == 20 - nl:
        return True
    else:
        return False
    
# Starting from tau, what is the next best mapping for gates in D?
# Before calling this function, usually we should have removed all topgates that are executable by tau
    # that is, greedy_solved_gates(tau,LD) should be empty

def good_next_mappingx(C,G,EG,nl,tau,LD,QFilter_type):
    # we assume that tau cannot solve any gates in D
    ''' Select the action (a sequence of swaps) and the mapping obtained by applying the action on tau

        Method: first check if tau can be extended so that some topgate is executable;
                then, for each result rho=[action,mapping] in SWAP3x, compute the efficiency rate t of rho, and append [rho,t] in SG if t>0;
                if SG is empty, then consider one more swap (i.e., consider SWAP4x), and obtain new SG;
                if SG is still empty, then call the fallback;
                otherwise, select the result with the highest t and return the corresponding action and the new mapping.
                The efficiency rate of rho is defined as (number of solved gates)/(3*(number of swaps used))

    Args:
        C (list): the input circuit
        G (graph): the architecture graph
        EG (list): the list of edges of G
        nl (int): the number of qubits in C
        tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
        LD (list): the index list of D, which represents the current logical circuit
        QFliter_type (str): '0' (no filter), '1' (Q0-filter as in the paper), '12' (Q0+Q1-filter),\
                                    '12x' (Q0+weak_Q1-filter), '2x' (weak_Q1-filter)        
    Returns:
        best_action (list), best_mapping (list), Fallback (Bool),
        2nd_best_action (list), 2nd_best_mapping (list), 2nd_best_value >= 2/3*best_value (Bool)
    '''
    
    LTG = topgates(C,nl,LD)
    TG = list(C[i] for i in LTG)
    if not LD:
        print('This is wrong') 

    # we first check if tau can be extended so that some topgate is executable    
    if not complete_mapping(C,nl,tau): 
##        print('The mapping %s is incomplete %s' %(tau,complete_mapping(tau)))
        UU = swap_reduce_min_dist(G,tau,TG) # UU has form [taunew, action]
        if not UU[1]:
##            print('  We extend the incomplete mapping to an unoccupied physical qubit.')
            return UU[1], UU[0], False, [], [], False
            
    # rho in SWAP3x has form [[edge_1,...],taunew], where tau is a mapping, [edge_1,...] is a sequence of swaps
    SG = [] # elements in SG has form [rho, t], where rho in SWAP3x and t the efficiency rate
    for rho in SWAP3x(C,EG,nl,tau,LD,QFilter_type):
        u = len(rho[0])
        x = greedy_solved_gates(C,G,EG,nl,rho[1],LD) # x has form [solved_gates]
        if not x:
            continue
        t = len(x)/(3*u) # the efficiency rate of rho is defined as (number of solved gates)/(3*(number of swaps used))
        SG.append([rho,t])

    if not SG:
        print('Yes, we go to the fourth level')
        for rho in SWAP4x(C,EG,nl,tau,LD,QFilter_type):
            u = len(rho[0]) # u < 3 is possible
            x = greedy_solved_gates(C,G,EG,nl,rho[1],LD)
            if not x:
                continue
            t = len(x)/(3*u) # the efficiency rate of rho is defined as (number of solved gates)/(3*(number of swaps used))
            SG.append([rho,t])
            
    taunew = tau[:]
    
    if not SG: # if no proper mapping in SWAP1-4, take a random neighbor of tau
        UU = swap_reduce_min_dist(G,tau,TG) # UU has form [taunew, action]
        print('  Oops, there are no good swaps, we select one swap that reduces the minimal distance')
        return UU[1], UU[0], True, [], [], False
    
    # else, we select the mapping with the best efficiency rate from SG
    SG.sort(key=lambda t:t[1], reverse=True)
    sg = SG[0] # with form [rho, val], where rho has form [ [edge1,...], taunew]
    action = sg[0][0] # with form [edge1,...], where each edge1 is a swap [p,q]
    taunew = sg[0][1]

    #added on 2019-12-21 to address the fallback concern
        # put the next best action and mapping in the end,
        # and set a Bool parameter to indict whether Fallback would be activated
    sgx = SG[1]
    if sgx[1] < 2/3*sg[1]: #the next best action has significant worse value
        return action, taunew, False, [], [], False
    actionx = sgx[0][0]
    taunewx = sgx[0][1]
                          
    return action, taunew, False, actionx, taunewx, True              

##################################################################
      #\__/#\__/#\#/~\ THE MAIN ALGORITH10M \__/#\__/#\#/~\
##################################################################

def qct(C,G,EG,tau,QFilter_type):
    ''' Transform a quantum circuit so that it is executable in Q20 with initial mapping tau

        Method: while there are unprocessed CNOT gates in the logical circuit, first process all gates that can be solved by the current mapping,
                then compute the next mapping and apply the corresponding action to transform to the next mapping.

    Args:
        C (list): the input circuit
        G (graph): the architecture graph
        EG (list): the list of edges of G
        tau (list): a mapping from V to Q, where V (Q) is the list of physical (logical) qubits
        QFliter_type (str): '0' (no filter), '1' (Q0-filter as in the paper), '12' (Q0+Q1-filter),\
                                    '12x' (Q0+weak_Q1-filter), '2x' (weak_Q1-filter)         
    Returns:
        C_out (list): the output circuit
        cost_time: time (seconds) required
    '''

    QIC = qubit_in_circuit(C)
    nl = len(QIC) # number of qubits in C
    l = len(C)

    #\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\
    #print('The input circuit C has %s CNOT gates and uses %s qubits in %s.' %(l,nl,QIC))
    #print('  We transform C into a physical circuit executatble in IBM Q Tokyo.')
    #print(' Here tau maps a physical qubit in G to a logic qubit and the i-th physical qubit is free if tau[i] is 20') 
    #\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\#\__/#\#/\#\__/#\#/\__/--\__/#\__/#\#/~\

    L = list(range(l)) # the index set of C
    nsvg = 0 # the number of executed gates
    COST = 0 # the number of auxiliary CNOT gates
    taunew = tau[:]
    L1 = L[:]
    nr = 0
    C_out = [] # the output circuit

    start = time.time()
    while nsvg < l:

        # first check what gates can be executed by taunew (obtained in the previous round),
         # then remove their indices from L1
        GSG = greedy_solved_gates(C,G,EG,nl,taunew,L1)       
        for i in GSG:
            # if C[i] is executed, remove i from L1 and add CNOT gate [p,q] to C_out
            L1.remove(i)
            p = taunew.index(C[i][0])
            q = taunew.index(C[i][1])
            C_out.append([p,q])
                    
        nsvg += len(GSG) # the number of solved gates

##        print('Round %s solved %s gates' %(nr, len(GSG)))
##        print('Round %s solved gates in %s' %(nr, list(C[i] for i in GSG)))

        if not L1:
##        print('We have solved all gates within %s rounds' %(nr+1))
            break
        
        nr += 1
        # modified on 2019-12-21 to address the fallback concern
        # simulate the process by looking ahead one step: if the best action would activate a fallback then try the second best action

        Y = good_next_mappingx(C,G,EG,nl,taunew,L1) # returns (best_action, best_mapping, Fallback, 2nd_best_action, 2nd_best_mapping)
        taux_temp = Y[1] # the temp new mapping resulted from the best action
        tau_temp = taunew[:]
        L1_temp = L1[:]
        GSG_temp = greedy_solved_gates(C,G,EG,nl,taux_temp,L1_temp)       
        for i in GSG_temp: 
            L1_temp.remove(i)

        if L1_temp and Y[5]: #Y[5] denotes that the 2nd best (action) value is at least 2/3 of the best (action) value
            Y_tempx = good_next_mappingx(C,G,EG,nl,taux_temp,L1_temp) # look-ahead one more step by using the best action
            if Y_tempx[2] == True: # Fallback would be activated
                tauy_temp = Y[4] # the mapping resulted from the 2nd best action
                L2_temp = L1[:]
                Y_tempy = good_next_mappingx(C,G,EG,nl,tauy_temp,L2_temp) # look-ahead one more step by using the 2nd best action
                if Y_tempy[2] == False: #Fallback would not be activated for the 2nd best action, we select the 2nd action
                    action = Y[3]
                    taunew = Y[4]
        
        action = Y[0] # a sequence of swaps that transforms the previous mapping into the mapping tuanew below
        taunew = Y[1]
        cost = len(action)*3
##    print('Round %s: The next mapping' %nr)
##    print('   %s' %taunew)
##    print('   is obtained by swaps in %s are enforced with cost %s' %(action,cost))
            
        COST += cost
        for edge in action:
            # implement in C_out the swap corresponding to edge
            C_out.append([edge[0],edge[1]])
            C_out.append([edge[1],edge[0]])
            C_out.append([edge[0],edge[1]])
            
    end = time.time()
    cost_time = end-start
    #print('Search process takes time: %s seconds.' %cost_time)
    #print('We need to add %s auxiliary CNOT gates' %COST)
    #print(l,nsvg) # check if l=nsvg
    #print('The final mapping is %s' %taunew)
    ##print('The output circuit contains %s gates' %len(C_out))
    return C_out, cost_time
   
