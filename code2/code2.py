import networkx as nx
import statistics
import random
import math
import matplotlib.pyplot as plt
import sys

#initializing lists for nodes and edges
node_list = []
edge_list = []
avg_edge_list = []

number_of_simulations = int(sys.argv[5])
for simulations in range(number_of_simulations):
    #initializing the graph
    grid_x = int(sys.argv[1])
    grid_y = int(sys.argv[2])
    G = nx.triangular_lattice_graph(m = grid_x, n = grid_y, periodic=True, with_positions=True, create_using=None)
    for (w, x), (y, z) in G.edges:
        if (x == z):
            G.remove_edge((w, x), (y, z))

    #computing probabilities from interaction energies between edges
    R = 8.314 # J/mol
    T = 300 #Kelvin
    E1 = -611.236*1000/(R*T) #green color edge
    E2 = -398.779*1000/(R*T) #yellow color edge
    E3 = -367.681*1000/(R*T) #red color edge
    E4 = -523.476*1000/(R*T) #blue color edge
    E_correction = 0*1000/(R*T) #turned off if not accounting for local changes in node-node interaction

    #assigning energies to the edges
    for ((w, x), (y, z)) in G.edges():
        if x == 0:
            x_max = max(w, x)
        if y == 0:
            y_max = max(y, z)

    for j in range(0, y_max + 1, 2):
        for i in range(x_max + 1):
            if (1 + i) == (x_max + 1):
                G[(0 + i, 1 + j)][(0, 0 + j)]['weight'] = E3
                G[(0 + i, 1 + j)][(0, 0 + j)]['color'] = 'r'
                if (2 + j) == (y_max + 1):
                    G[(0 + i, 1 + j)][(0, 0)]['weight'] = E4
                    G[(0 + i, 1 + j)][(0, 0)]['color'] = 'b'
                else:
                    G[(0 + i, 1 + j)][(0, 2 + j)]['weight'] = E4
                    G[(0 + i, 1 + j)][(0, 2 + j)]['color'] = 'b'
            else:
                G[(0 + i, 1 + j)][(1 + i, 0 + j)]['weight'] = E3
                G[(0 + i, 1 + j)][(1 + i, 0 + j)]['color'] = 'r'
                if (2 + j) == (y_max + 1):
                    G[(0 + i, 1 + j)][(1 + i, 0)]['weight'] = E4
                    G[(0 + i, 1 + j)][(1 + i, 0)]['color'] = 'b'
                else:
                    G[(0 + i, 1 + j)][(1 + i, 2 + j)]['weight'] = E4
                    G[(0 + i, 1 + j)][(1 + i, 2 + j)]['color'] = 'b'
            G[(0 + i, 1 + j)][(0 + i, 0 + j)]['weight'] = E1
            G[(0 + i, 1 + j)][(0 + i, 0 + j)]['color'] = 'g'
            if (2 + j) == (y_max + 1):
                G[(0 + i, 1 + j)][(0 + i, 0)]['weight'] = E2
                G[(0 + i, 1 + j)][(0 + i, 0)]['color'] = 'y'
            else:
                G[(0 + i, 1 + j)][(0 + i, 2 + j)]['weight'] = E2
                G[(0 + i, 1 + j)][(0 + i, 2 + j)]['color'] = 'y'

    weights = [G[u][v]['weight'] for u, v in G.edges()]
    colors = [G[u][v]['color'] for u,v in G.edges()]

    #calculate initial partition function
    total_sum = 0
    for weight in range(len(weights)):
        total_sum += math.exp(weights[weight])

    #initialize fragmentation
    i = 0
    total_nodes = len(list(G.nodes))
    empty_nodes = []

    #run fragmentation
    number_of_runs = int(sys.argv[3])
    for i in range(total_nodes*number_of_runs):
        graph_list_nodes = list(G.nodes)
        graph_list_edges = list(G.edges)

        num1 = random.randint(1, len(graph_list_edges) - 1)

        weights = [G[u][v]['weight'] for u, v in G.edges()]
        if (total_sum <= 0):
            for weight in range(len(weights)):
                total_sum += math.exp(weights[weight])
        proxy = math.exp(list(G.edges(data=True))[num1][2]['weight'])
        probability = (math.exp(list(G.edges(data=True))[num1][2]['weight']))/(total_sum)

        weight_color = (list(G.edges(data=True))[num1][2]['color'])
        threshold = random.uniform(0, 1)

        #breaking an edge
        if (1 - probability) <= threshold:
            total_sum += -(math.exp(list(G.edges(data=True))[num1][2]['weight']))
            G.remove_edge(graph_list_edges[num1][0], graph_list_edges[num1][1])

            node1 = graph_list_edges[num1][0]
            node2 = graph_list_edges[num1][1]

            #recalculating the partition function
            for node1_neighbor in list(G.neighbors(node1)):
                total_sum += -math.exp(G.edges[node1, node1_neighbor]['weight'])
                G.edges[node1, node1_neighbor]['weight'] = (G.edges[node1, node1_neighbor]['weight']) + E_correction
                total_sum += math.exp(G.edges[node1, node1_neighbor]['weight'])
            
            for node2_neighbor in list(G.neighbors(node2)):
                total_sum += -math.exp(G.edges[node2, node2_neighbor]['weight'])
                G.edges[node2, node2_neighbor]['weight'] = (G.edges[node2, node2_neighbor]['weight']) + E_correction
                total_sum += math.exp(G.edges[node2, node2_neighbor]['weight'])
        
        for n in G.nodes:
            if len(list(G.neighbors(n))) == 0:
                if n in empty_nodes:
                    ("")
                else:
                    empty_nodes += [n]

        #fragmentation threshold
        fragmentation_threshold_fraction = int(sys.argv[4])    
        if len(empty_nodes) >= (total_nodes/fragmentation_threshold_fraction):
            break

    edge_list += [len(list(G.edges))]
    avg_edge_list += [statistics.mean(edge_list)]
    node_list += [total_nodes - len(empty_nodes)]

plt.plot(avg_edge_list)
plt.xlabel("Number of Simulations")
plt.ylabel("Average number of remaining graph edges")
plt.title("Evolution of edges for freeing a fixed number of nodes")
plt.savefig('code2.png')

