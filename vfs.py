#-*- coding:utf-8 -*-
# AUTHOR:   yaolili
# FILE:     vf.py
# ROLE:     vf2 algorithm
# CREATED:  2015-11-28 20:55:11
# MODIFIED: 2015-12-05 11:58:12
# ADDAPTED: 2019-06-25 for Q20 Mapping by Sanjiang Li (SL)
# all comments by SL started with '#sl'

import sys
import os
import networkx as nx
from maps import Map

class Vf:

    __origin = None
    __sub = None

    #sl Achtung! Here pairs are taken from subNMN and gNMN, not from subMN and gMN
    #sl Here subNMN reads as not-in-Map-neighbor-of-the-examined-vertex-in-the-subgraph
    def candidate(self, subNMNeighbor, gNMNeighbor):
        if not (subNMNeighbor and gNMNeighbor):
            print ("Class Vf candidate() arguments value error! subNMNeighbor or gNMNeighbor is empty!")
            exit()
        if not (isinstance(subNMNeighbor, list) and isinstance(gNMNeighbor, list)):
            print ("Class Vf candidate() arguments type error! type list expected!")
            exit()
        if not all(isinstance(x, int) for x in subNMNeighbor):
            print ("Class Vf candidate() arguments type error! int in subNMNeighbor list expected!")
        if not all(isinstance(x, int) for x in gNMNeighbor):
            print ("Class Vf candidate() arguments type error! int in gNMNeighbor list expected!")        
        
        pairs = []
        for x in subNMNeighbor:
            for y in gNMNeighbor:
                pairs.append([x,y])

        return pairs

    #type = 0, pre; type = 1, succ
    #sl divide the graph neighborhood of a vertex into two disjoint parts: pre (in map) and succ (not in map)
    def preSucc(self, vertexNeighbor, map, type):
        #vertexNeighbor and map can be empty
        if not (isinstance(vertexNeighbor, list) and isinstance(map, list)):
            print ("Class Vf preSucc() arguments type error! vertexNeighbor and map expected list!")
            exit()
        if not (type == 0 or type == 1):
            print ("Class Vf preSucc() arguments value error! type expected 0 or 1!")
           
        result = []
        #succ
        if type:
            for vertex in vertexNeighbor:
                if vertex not in map:                   
                    result.append(vertex)
        #pre
        else:
            for vertex in vertexNeighbor:
                if vertex in map:
                    result.append(vertex)
        return result


    def isMeetRules(self, v1, v2, subgraph, graph, result, subMap, gMap, subMNeighbor, gMNeighbor):
            
        if not result:
            return True        
        
        v1Neighbor = list(nx.all_neighbors(subgraph, v1))
        v2Neighbor = list(nx.all_neighbors(graph, v2))
                
        v1Pre = self.preSucc(v1Neighbor, subMap, 0)
        v1Succ = self.preSucc(v1Neighbor, subMap, 1)
        v2Pre = self.preSucc(v2Neighbor, gMap, 0)
        v2Succ = self.preSucc(v2Neighbor, gMap, 1)

        if len(v1Pre) > len(v2Pre) or len(v1Succ) > len(v2Succ):
            return False
                       
        for pre in v1Pre:
            if result[pre] not in v2Pre:
                return False
                       
        len1 = len(set(v1Neighbor) & set(subMNeighbor))
        len2 = len(set(v2Neighbor) & set(gMNeighbor))

        if len1 > len2 or len(set(subMNeighbor)) > len(set(gMNeighbor)) or len(subMNeighbor) > len(gMNeighbor):
            return False
               
        return True
        

    def dfsMatch(self, subgraph, graph, result):
#        print ("in dfsMatch() result %s: " %result)
        if not isinstance(result, dict):
            print ("Class Vf dfsMatch() arguments type error! result expected dict!")
        
        curMap = Map(result)
            
        if len(result) == len(nx.nodes(subgraph)):
            return result

        subMNeighbor = curMap.neighbor(subgraph, 0)
        gMNeighbor = curMap.neighbor(graph, 1)   

        if not subMNeighbor:
            X = list(set(nx.nodes(subgraph)) - set(curMap.subMap()))
        else:
            X = subMNeighbor
            
        subNMN_deg = list([nx.degree(subgraph, v), v] for v in X)
        subNMN_deg.sort(key=lambda t: t[0], reverse=True)
        subNMNeighbor = [ subNMN_deg[0][1] ]

        if not gMNeighbor:
            gNMNeighbor = list(set(nx.nodes(graph)) - set(curMap.gMap()))
        else: 
            gNMNeighbor = gMNeighbor[:]
        
        pairs = self.candidate(subNMNeighbor, gNMNeighbor)        
        if not pairs:
            return result
                
        for pair in pairs:
            v1 = pair[0]
            v2 = pair[1]
            if(self.isMeetRules(v1, v2, subgraph, graph, result, curMap.subMap(), curMap.gMap(), subMNeighbor, gMNeighbor)):
                result[v1] = v2
                self.dfsMatch(subgraph, graph, result)

                if len(result) == len(nx.nodes(subgraph)):
                    return result
                result.pop(v1)
                
        #sl the procedure stops when it either constructs a complete mapping or finds out that the current result is incompletable
##        print(result)
        return result   
