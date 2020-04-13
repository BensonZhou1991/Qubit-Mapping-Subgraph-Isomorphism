import networkx as nx

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
