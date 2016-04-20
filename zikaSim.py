#!/usr/bin/python3
"""
zikaSim.py is simulator for infection spreading by air travel and at major
metro areas.  The goal of this simulation is to model the Zika infection spread and
different experiments

# test edge-based quarantine strategies for the given network.

Usage:
    simulator.py -ms [c=<IATA>] [d=<start infection data]
        [r1=<start simulation>] [r2=<end simulation>]
        <airport database> <route database> <mosquito curves database>

Flags:
    -m: Show map of infection
    -s: Stats on infection at starting city

Option:
    --c         The IATA code for hub
    --d         Start day of infection
    --r1        Start of simulation
    --r2        End of simulation
    --t         Ro - Reproduction number for particular virus
"""

# Title:  zikaSim.py
# Updated Authors: Tilak Patel and Derrick Williams
# Original Authors: Nicholas A. Yager and Matthew Taylor
# Date:   2016-04-12

# python zikaSim.py -s --c ATL --d 125 --r1 120 --r2 180 ./Data/airportsMin.csv ./Data/airlineRoutesPassengerData.csv ./Data/mosCurves.csv
# python zikaSim.py -s --d 125 --r1 120 --r2 365 --t 4 ./Data/airportsMin.csv ./Data/airlineRoutesPassengerData.csv ./Data/mosCurves.csv



import copy
import getopt
import math
import networkx as nx
import matplotlib.pyplot as plt
import operator
import os
import random
import sys
from scipy import stats
import time
from mpl_toolkits.basemap import Basemap
import queue


# GLOBAL
MAP = False
TAU       = 4      # Zika virus "Ro", DEFAULT
routeInfo = dict()
approvedAirports = dict()
airportsToInfect = dict()
I = dict()
S = dict()
V = dict()
R = dict()
timeStepsTracker = list()
CITY_TO_INFECT = "ATL"
DATE_TO_INFECT = 1
RANGE_BEGIN = 1
RANGE_END   = 365
STAT = False


def main():
    """
    Primary function that initiates network creation and handles execution of
    infection simulations.

    Args:
        argv: A list of command-line arguments passed to the application.

    Returns:
        Void

    """
    global MAP, CITY_TO_INFECT, RANGE_BEGIN, RANGE_END, DATE_TO_INFECT, STAT,TAU
    global timeStepsTracker

    # Determine the parameters of the current simulation.
    opts, args = getopt.getopt(sys.argv[1:], "ms", ["c=", "d=", "r1=", "r2=",
                                                    "t="])

    # Check if the data arguments are available
    if len(args) < 3:
        print(__doc__)
        exit()

    AIRPORT_DATA = args[0]
    ROUTE_DATA = args[1]
    MOSQUITO_CURVES = args[2]

    # simulations = list()

    for opt, par in opts:
        # map of network
        if opt == "-m":
            MAP = True
        # stats on infection city
        elif opt == "-s":
            STAT = True
        # infect a particular city
        elif opt == "--c":
            CITY_TO_INFECT = par
        # date to infect
        elif opt == "--d":
            DATE_TO_INFECT = int(par)
        # Beginning part of simulation
        elif opt == "--r1":
            RANGE_BEGIN = int(par)
        # Ending part of simulation
        elif opt == "--r2":
            RANGE_END = int(par)
        # A different tau for mosquito dynamics
        elif opt == "--t":
            TAU = float(par)


    # Create the network using the command arguments
    network = create_network(AIRPORT_DATA, ROUTE_DATA, MOSQUITO_CURVES)
    # print(network.nodes(network))

    # Setup Global SIVR stat tracker
    setupGlobalSIVR()

    # Run infection simulation
    for i in range(RANGE_BEGIN,RANGE_END):
        timeStepsTracker.append(i)
        infection(network, i)

    # Visualize network
    if MAP:
        visualize(network)

    # Stats of infection
    if STAT:
        for node in I:
            if node == CITY_TO_INFECT:

                i, = plt.plot(timeStepsTracker, I[node],label="I")
                s, = plt.plot(timeStepsTracker, S[node],label='S')
                r, = plt.plot(timeStepsTracker, R[node],label='R')
                plt.legend(handles=[i,s,r])

        plt.title(CITY_TO_INFECT + " Infection Dynamics")
        plt.xlabel('Days of Year')
        plt.ylabel('People')
        plt.xlim(RANGE_BEGIN,RANGE_END)
        plt.show()







def create_network(nodes, edges, curves):
    """
    Create a NetworkX graph object using the airport and route databases.

    Args:
        nodes: The file path to the nodes .csv file.
        edeges: The file path to the edges .csv file.
        curves: The file path to the mosquito curves .csv file.

    Returns:
        G: A NetworkX Graph object populated with the nodes and edges assigned
           by the data files from the arguments.

    """
    global routeInfo
    global approvedAirports

    print("Creating network.")
    G = nx.Graph()

    print("\tLoading airports", end="")
    sys.stdout.flush()

    # Load mosquito curves
    mosquitoCurves = dict()
    with open(curves, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            entries = line.split(",")
            for i in range(2,14):
                entries[i] = float(entries[i])
            mosquitoCurves[entries[1]] = entries[2:14]

    # Populate the graph with nodes.
    with open(nodes, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")

            G.add_node(int(entries[0]),
                       name=entries[1],
                       IATA=entries[2],
                       lat=entries[3],
                       lon=entries[4],
                       pop=int(entries[5]),
                       pos=(float(entries[3]),float(entries[4])),
                       I=[0,0],  # First num - Total I, second num new I
                       Iair=0,
                       S=int(entries[5]),       #createHumans(int(entries[5])),
                       V=0,
                       R=0,
                       MOS=mosquitoCurves[entries[2]]
                       )

    print("\t\t\t\t\t[Done]")

    print("\tLoading routes",end="")
    sys.stdout.flush()

    # Populate the graph with edges.
    edge_count = 0
    error_count = 0
    duplicate_count = 0
    line_num = 1
    with open(edges, 'r', encoding="utf-8") as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")
            routeInfo[(entries[0],entries[1],entries[2],entries[3])] = \
                float(entries[4])
            approvedAirports[entries[0]] = int(entries[1])
            try:
                if G.has_edge(int(entries[1]),int(entries[3])) or \
                    G.has_edge(int(entries[3]),int(entries[1])):
                    duplicate_count += 1
                else:
                    if line_num > 1:
                        from_vertex = int(entries[1])
                        to_vertex = int(entries[3])
                        G.add_edge(from_vertex, to_vertex )
                        G.edge[from_vertex][to_vertex]['IATAFrom'] = entries[0]
                        G.edge[from_vertex][to_vertex]['IATATo'] = entries[2]
                        edge_count += 1
            except ValueError:
                # The value doesn't exist
                error_count += 1
                pass
            line_num += 1

    print("\t\t\t\t\t\t[Done]")



    # Calculate the edge weights
    # print("\tCalculating edge weights",end="")
    # G = calculate_weights(G)
    # print("\t\t\t\t[Done]")

    # Add clustering data
    print("\tCalculating clustering coefficents",end="")
    cluster_network = nx.Graph(G)
    lcluster = nx.clustering(cluster_network)
    for i,j in G.edges():
        cluster_sum = lcluster[i] + lcluster[j]
        G[i][j]['cluster'] = cluster_sum
    print("\t\t\t[Done]")

    return G






def infection(input_network, timeStep):
    global I,S,V,R
    # mosCurve = [0,0,0,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0] # phx,
    approxMonth = timeStep // 31


    # Spread disease to other Airports
    currentNodes = input_network.node
    for node in input_network.nodes_iter(input_network):
        for key in routeInfo:
            if node[1]["IATA"] == key[2]:
                nodeDetails = currentNodes[int(key[1])]

                node[1]["Iair"] += int(nodeDetails["I"][1] / \
                                   nodeDetails["pop"] * \
                                   routeInfo[key] / 365)


    #print (input_network.node)

    # Infection simulation at hubs
    for node in input_network.nodes_iter(input_network):
        #if node[1]["IATA"] == "ATL":
        #    print(node[1]["I"])
        #  Record stats
        I[node[1]["IATA"]].append(node[1]["I"][0])
        S[node[1]["IATA"]].append(node[1]["S"])
        V[node[1]["IATA"]].append(node[1]["V"])
        R[node[1]["IATA"]].append(node[1]["R"])

        # Check for recovery
        if node[1]["I"][0] > 0:
            if timeStep - node[1]["I"][2][0] >= 7:
                group = node[1]["I"].pop(2)
                node[1]["I"][0] -= group[1]
                node[1]["R"]    += group[1]

        # Infect cities
        if timeStep == DATE_TO_INFECT and node[1]["IATA"] == CITY_TO_INFECT:
            infectCity(input_network)

        newlyInfected = min(math.ceil(TAU * node[1]["MOS"][approxMonth] *
                            (node[1]["I"][1] +
                              node[1]["Iair"])),node[1]["S"])

        # if node[1]["IATA"] == "ATL":
        #     print (min(math.ceil(TAU * node[1]["MOS"][approxMonth] *
        #                     (node[1]["I"][1] + node[1]["Iair"])),
        #                     node[1]["S"]))
        if newlyInfected > 0:
            node[1]["S"] -= newlyInfected
            node[1]["I"].append((timeStep,newlyInfected))
            node[1]["I"][0] += newlyInfected
            node[1]["I"][1] = newlyInfected

        # Remove temporary airport visitors
        node[1]["Iair"] = 0

        # print ("newlyInfected",newlyInfected)
        # print ("currentlyInfected",node[1]["I"])
        # print ("currentlyRecovered",node[1]["R"])


def infectCity(input_network):
    for node in input_network.nodes_iter(input_network):
        if node[1]["IATA"] == CITY_TO_INFECT:
            node[1]["S"] -= 1
            node[1]["I"].append((DATE_TO_INFECT,1))
            node[1]["I"][0] += 1
            node[1]["I"][1] += 1

            #print("infected",node[1]["I"])

def setupGlobalSIVR():
    global I,S,V,R

    for airport in approvedAirports:
        I[airport] = list()
        S[airport] = list()
        V[airport] = list()
        R[airport] = list()




def visualize(network):
    print("-- Starting to Visualize --")

    m = Basemap(
        projection='merc',
        ellps='WGS84',
        llcrnrlon=-160,urcrnrlon=-60,llcrnrlat=10,urcrnrlat=80,
        resolution="l"
        )

    pos = dict()

    for pos_node in network.nodes():
        # Normalize the lat and lon values
        x,y = m(float(network.node[pos_node]['lon']),
                float(network.node[pos_node]['lat']))

        pos[pos_node] = [x,y]

    #m.drawmapboundary("aqua")
    #m.fillcontinents('#555555')
    m.drawlsmask(land_color='coral',ocean_color='aqua',lakes=True)
    #m.bluemarble()

    # First pass - Green lines
    nx.draw_networkx_edges(network,pos,edgelist=network.edges(),
            width=1,
            edge_color="green",
            alpha=0.5,
            arrows=False)

    nx.draw_networkx_nodes(network,
            pos,
            linewidths=1,
            node_size=10,
            with_labels=False,
            node_color = "green")

    #m.bluemarble()
    #plt.title=title

    # Adjust the plot limits
    cut = 1.05
    xmax = cut * max(xx for xx,yy in pos.values())
    xmin =  min(xx for xx,yy in pos.values())
    xmin = xmin - (cut * xmin)


    ymax = cut * max(yy for xx,yy in pos.values())
    ymin = (cut) * min(yy for xx,yy in pos.values())
    ymin = ymin - (cut * ymin)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    plt.axis('off')
    plt.show()
    plt.close()




    #seed = 100
    #random.seed(seed)

        # if o == "-b":
        #     simulations.append("betweenness")
        # elif o == "-r":
        #     simulations.append("random")
        # elif o == "-s":
        #     simulations.append("sir")
        # elif o == "-c":
        #     simulations.append("clustering")
        # elif o == "-v":
        #     VISUALIZE = True
        # elif o == "-y":
        #     RECALCULATE = False
        # elif o == "--delay":
        #     DELAY = int(a)
        # elif o == "--nsim":
        #     NUM_SIMULATIONS = int(a)




    # Generate target-selection weights, and choose target vertices to infect.
    # Dict of node id with degrees per node
    # degrees = network.degree()

    #pos = nx.get_node_attributes(network,'pos')
    # print (pos)






    #nx.draw(network, pos)
    #plt.savefig("hello.png")

    # weights = dict()
    # for airport, degree in degrees.items():
    #     weights[airport] = network.out_degree(airport) +\
    #                        network.in_degree(airport)
    # targets = list()
    # for ind in range(0,NUM_SIMULATIONS):
    #     target_round = list()
    #     while len(target_round) < 10:
    #          chosen_airport = weighted_random(weights)
    #          if chosen_airport not in target_round:
    #              target_round.append(chosen_airport)
    #     targets.append(target_round)
    #
    #
    # # Make a directory for the data, and change into that directory.
    # currenttime = time.strftime("%Y-%m-%dT%H%M%S", time.gmtime())
    # os.makedirs(currenttime)
    # os.chdir(currenttime)
    #
    # # Record relevent data about the simulation.
    # # TMP simulation_data(network, currenttime, target, seed)
    #
    #
    # # Prepare simulations
    #
    # # Remove some edges as per necessary from the pool of edges that can be
    # # used as cancelled edges.
    #
    # edgepool = network.edges(data=True)
    # # if INTERNATIONAL:
    # #     for i,j,data in edgepool:
    # #         if data["international"] == False:
    # #             edgepool.remove((i,j,data))
    # #         index += 1
    # # elif DOMESTIC:
    # #     for i,j,data in edgepool:
    # #         if data["domestic"] == False:
    # #             degrees.edgepool((i,j,data))
    # #         index += 1
    #
    #
    # for strategy in simulations:
    #
    #     print("{0} Mode.".format(strategy) )
    #
    #     index = 0
    #
    #     # Generate a list a sorted list of flights to cancel based on the
    #     # strategy.
    #
    #     cancellist = list()
    #     if strategy == "random":
    #         # Sort the edges randomly
    #
    #         cancellist = random.sample(edgepool, len(edgepool))
    #
    #     elif strategy == "clustering":
    #         # Sort the edges based on the sum of the clustering coefficent.
    #
    #         sorted_cluster = sorted(edgepool, key=lambda k: k[2]['cluster'],
    #                         reverse=True)
    #         for cluster_item in sorted_cluster:
    #             if network[cluster_item[0]][cluster_item[1]]['cluster'] < 2:
    #                 if network[cluster_item[0]][cluster_item[1]]['cluster'] > 0:
    #                     cancellist.append((cluster_item[0], cluster_item[1]))
    #
    #     elif strategy == "betweenness":
    #         # Sort the edges based on weighted edge-betweenness.
    #
    #         betweennesses = nx.edge_betweenness_centrality(network,
    #                                                        weight="weight")
    #         cancellist = sorted(betweennesses.keys(),
    #                             key=lambda k: betweennesses[k], reverse=True)
    #
    #
    #     print(cancellist[:20])
    #     # Make a new folder for the data.
    #     os.makedirs(strategy)
    #
    #     iteration = 0
    #     efforts = [0]
    #     efforts.extend(range(1,101,5))
    #     for target in targets:
    #
    #         # Open a file for this targets dataset
    #         output_file = open("{0}/{0}_{1}.csv".format(strategy,
    #                                                     pad_string(iteration,4)
    #                                                     ),"w")
    #         output_file.write('"effort","total_infected, edges_closed"\n')
    #
    #
    #         for effort in efforts:
    #             if effort != 0:
    #                 max_index = int(len(cancellist) * (effort/100))-1
    #                 cancelled = cancellist[0:max_index]
    #             else:
    #                 cancelled = None
    #
    #             title = "{0} - {1}%".format(strategy, effort/100)
    #             results = infection(network, cancelled, target, vis=VISUALIZE,
    #                                 title=title, DELAY=DELAY)
    #             total_infected = results["Infected"] + results["Recovered"]
    #             output_file.write("{0},{1}\n".format(effort/100,total_infected))
    #
    #             if total_infected == 1:
    #                 for remaining_effort in range(effort+5,101,5):
    #                     output_file.write("{0},{1}\n".format(remaining_effort/100,
    #                                                           total_infected))
    #                 break
    #
    #         iteration += 1
    #         output_file.close()



def weighted_random(weights):
    number = random.random() * sum(weights.values())
    for k,v in weights.items():
        if number <= v:
            break
        number -= v
    return k

def pad_string(integer, n):
    """
    Add "0" to the front of an integer so that the resulting string in n
    characters long.

    Args:
        integer: The number to pad.
        n: The desired length of the string

    Returns
        string: The padded string representation of the integer.
        
    """

    string = str(integer)

    while len(string) < n:
        string = "0" + string

    return string

def simulation_data(network, time, targets, seed):
    """
    Output various statistics of the nature of the network to a file, including
    the diameter, the number of verticies and edges, and the
    average in and out degrees.

    Args:
        network: A NetworkX network graph.

    Returns:
        VOID

    IO:
        network.dat: A file with all of the relevant network information.

    """

    print("\tDetermining network type.")
    # Determine if the graph is directed or undirected
    if isinstance(network,nx.DiGraph):
        network_type = "Directed"
    else:
        network_type = "Undirected"

    print("\tCalculaing edges and verticies.")
    # Number of verticies and edges
    edges = network.number_of_edges()
    verticies = network.number_of_nodes()

    
    # Not every vertex can lead to every other vertex.
    # Create a subgraph that can.
    print("\tTemporarily converting to undirected.")
    undirected = network.to_undirected()
    print("\tFinding subgraphs.")
    subgraphs = nx.connected_component_subgraphs(undirected)


    # Find the number of vertices in the diameter of the network

    print("\tFinding network diameter.")
    diameter = nx.diameter(subgraphs[0])


    print("\tStoring network parameters")

    data_file = open("network.csv", "w")
    data_file.write("Simulation name: {0}\n\n".format(time))
    data_file.write("Network properties\n===============\n")
    data_file.write("Network type: {0}\n".format(network_type))
    data_file.write("Number of verticies: {0}\n".format(verticies))
    data_file.write("Number of edges: {0}\n".format(edges))
    data_file.write("Diameter: {0}\n".format(diameter))

    data_file.close()

def create_network2(nodes, edges):
    """
    Create a NetworkX graph object using the airport and route databases.

    Args:
        nodes: The file path to the nodes .csv file.
        edeges: The file path to the edges .csv file.

    Returns:
        G: A NetworkX DiGraph object populated with the nodes and edges assigned
           by the data files from the arguments.
           
    """

    print("Creating network.")
    G = nx.DiGraph()

    print("\tLoading airports", end="")
    sys.stdout.flush()
    # Populate the graph with nodes.
    with open(nodes, 'r', encoding='utf-8') as f:
        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")

            G.add_node(int(entries[0]),
                       name=entries[1],
                       IATA=entries[2],
                       lat=entries[3],
                       lon=entries[4],
                       #pop=entries[5],
                       pos=(float(entries[3]),float(entries[4])),
                       I=[0],
                       S=int(entries[5]),       #createHumans(int(entries[5])),
                       V=0,
                       R=0
                       )
            # G.add_node(int(entries[0]),
            #            country=entries[3],
            #            name=entries[1],
            #            lat=entries[6],
            #            lon=entries[7])

    #for node in G.nodes_iter(G):
    #    print (len(node[1]["S"]))
    print("\t\t\t\t\t[Done]")
    
    print("\tLoading routes",end="")
    sys.stdout.flush()
    # Populate the graph with edges.
    edge_count = 0
    error_count = 0
    duplicate_count = 0
    line_num = 1
    with open(edges, 'r', encoding="utf-8") as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")
            try:
                # if G.has_edge(int(entries[3]),int(entries[5])):
                if G.has_edge(int(entries[1]),int(entries[3])):
                    duplicate_count += 1
                else:
                    if line_num > 1:
                        # from_vertex = int(entries[3])
                        # to_vertex = int(entries[5])
                        # G.add_edge(from_vertex, to_vertex )
                        # G.edge[from_vertex][to_vertex]['IATAFrom'] = entries[2]
                        # G.edge[from_vertex][to_vertex]['IATATo'] = entries[4]
                        from_vertex = int(entries[1])
                        to_vertex = int(entries[3])
                        G.add_edge(from_vertex, to_vertex )
                        G.edge[from_vertex][to_vertex]['IATAFrom'] = entries[0]
                        G.edge[from_vertex][to_vertex]['IATATo'] = entries[2]
                        edge_count += 1
            except ValueError:
                # The value doesn't exist
                error_count += 1
                pass
            line_num += 1
    
    print("\t\t\t\t\t\t[Done]")

    # Limit to the first subgraph
    # print("\tFinding largest subgraph",end="")
    # undirected = G.to_undirected()
    # subgraphs = nx.connected_component_subgraphs(undirected)
    # print (subgraphs)
    # subgraph_nodes = subgraphs[0].nodes()
    # to_remove = list()
    # for node in G.nodes():
    #     if node not in subgraph_nodes:
    #         to_remove.append(node)
    # G.remove_nodes_from(to_remove)
    print("\t\t\t\t[Done]")

    
    # print("\tRemoving isolated vertices",end="")
    # # Remove nodes without inbound edges
    # indeg = G.in_degree()
    # outdeg = G.out_degree()
    # to_remove = [n for n in indeg if (indeg[n] + outdeg[n] < 1)]
    # G.remove_nodes_from(to_remove)
    # print("\t\t\t\t[Done]")



    # Calculate the edge weights
    print("\tCalculating edge weights",end="")
    G = calculate_weights(G)
    print("\t\t\t\t[Done]")
    
    # Add clustering data
    print("\tCalculating clustering coefficents",end="")
    cluster_network = nx.Graph(G)
    lcluster = nx.clustering(cluster_network)
    for i,j in G.edges():
        cluster_sum = lcluster[i] + lcluster[j]
        G[i][j]['cluster'] = cluster_sum
    print("\t\t\t[Done]")

    # Flag flights as domestic or international.
    # print("\tCategorizing international and domestic flights",end="")
    # for i,j in G.edges():
    #     if G.node[i]["country"] == G.node[j]['country']:
    #         G[i][j]['international'] = False
    #     else:
    #         G[i][j]['international'] = True
    # print("\t\t[Done]")

    return G

def calculate_weights(input_network):
    """
    Add weights to the edges of a network based on the degrees of the connecting
    verticies, and return the network.

    Args:
        input_network: A NetworkX graph object
    Returns:
        G: A weighted NetworkX graph object.
    """
    
    G = input_network.copy()

    # Add weights to edges
    for node in G.nodes():
        successors = G.successors(node)
        weights = dict()

        # Calculate the total out degree of all succs
        total_degree = 0
        for successor in successors:
  
            try:
                total_degree += G.out_degree(successor)
            except TypeError:
                # Don't add anything
                pass

        # Find the weight for all possible successors
        for successor in successors:
            successor_degree = G.out_degree(successor)

            try:
                int(successor_degree)
            except TypeError:
                successor_degree = 0

            if total_degree > 0:
                probability_of_infection = successor_degree / \
                                           total_degree
            else:
                probability_of_infection = 0

            weights[successor] = probability_of_infection
        
        largest_weight = 0
        smallest_weight = 2
        for successor, weight in weights.items():
            if weight > largest_weight:
                largest_weight = weight
            elif weight < smallest_weight:
                smallest_weight = weight
        #(strat.shared_fitness - lowest_fitness) / \
        #                       (highest_fitness - lowest_fitness)

        for successor in successors:
            if largest_weight != smallest_weight:
                relative_weight = (weights[successor] - smallest_weight) /\
                                  (largest_weight - smallest_weight)
            else:
                relative_weight = 0
            G[node][successor]['weight'] = relative_weight

    return G






def infection2(input_network, vaccination, starts,DELAY=0, vis = False,
              file_name = "sir.csv", title="",  RECALCULATE = True):
    """
    Simulate an infection within network, generated using seed, and with the
    givin vaccination strategy. This function will write data from each timestep
    to "data.csv".

    Args:
        network: A NetworkX DiGraph object.
        vaccination: A list of node indices to label as recovered from the 
                     begining.

    Returns:
        state: A dictionary of the total suscceptable, infected, and recovered.

    """

    print("Simulating infection.")

    network = input_network.copy()
    
    # Recalculate the weights of the network as per necessary

    # Open the data file
    f = open(file_name, "w")
    f.write("time, s, e, i, r\n")

    # Set the default to susceptable
    sys.stdout.flush()
    for node in network.nodes():
        network.node[node]["status"] =  "s"
        network.node[node]["color"] = "#A0C8F0"
        network.node[node]["age"] = 0
    
    # Assign the infected
    for start in starts:
        infected = start
        network.node[infected]["status"] = "i"
        network.node[infected]["color"]  = "green"
        
        if isinstance(network,nx.DiGraph):
            in_degree = network.in_degree()[infected] 
            out_degree = network.out_degree()[infected]
            degree = in_degree + out_degree
        else:
            degree = network.degree()[infected]

        print("\t",network.node[infected]["name"],"[",degree,"]")


    if vaccination is not None:
        print("\tVaccinated: ", len(vaccination) )
    else: 
        print("\tVaccinated: None")

    if vis:
        pos = nx.spring_layout(network, scale=2)

    # Iterate through the evolution of the disease.
    for step in range(0,99):
        # If the delay is over, vaccinate.
        # Convert the STRING! 
        if int(step) == int(DELAY):
            if vaccination is not None:
                print(DELAY,"on step",step)
                network.remove_edges_from(vaccination)
                # Recalculate the weights of the network as per necessary
                if RECALCULATE == True:
                    network = calculate_weights(network)


        # Create variables to hold the outcomes as they happen
        S,E,I,R = 0,0,0,0

        for node in network.nodes():
            status = network.node[node]["status"]
            age = network.node[node]["age"]
            color = network.node[node]["color"]

            if status is "i" and age >= 11:
                # The infected has reached its recovery time
                network.node[node]["status"] = "r"
                network.node[node]["color"] = "purple"
                
            if status is "e" and age >= 3 and age < 11:
                # The infected has reached its recovery time
                network.node[node]["status"] = "i"
                network.node[node]["color"] = "green"

            elif status is "e":
                network.node[node]["age"] += 1

            elif status is "i":
                # Propogate the infection.
                if age > 0:
                    victims = network.successors(node)
                    number_infections = 0
                    for victim in victims:
                        infect_status = network.node[victim]["status"]
                        infect = False # Set this flag to False to start 
                                       # weighting.


                        if random.uniform(0,1) <= network[node][victim]['weight']:

                            infect = True
                            number_infections+=1

                        if infect_status == "s" and infect == True:
                            network.node[victim]["status"] = "e"
                            network.node[victim]["age"] = 0
                            network.node[victim]["color"] = "#FF6F00"
                network.node[node]["age"] += 1


        # Loop twice to prevent bias.
        for node in network.nodes():
            status = network.node[node]["status"]
            age = network.node[node]["age"]
            color = network.node[node]["color"]

            if status is "s":
                # Count those susceptable
                S += 1

            if status is "e":
                E += 1

            if status is "v":
                S += 1

            elif status is "r":

                R += 1

            elif status is "i":
                
                I += 1
        print("{0}, {1}, {2}, {3}, {4}".format(step, S, E, I, R))

        printline = "{0}, {1}, {2}, {3}, {4}".format(step, S, E, I, R)
        f.write(printline + "\n")

       # print("\t"+printline)

        if I is 0:
            break

        if vis:
            #write_dot(network, title+".dot")
            visualize(network, title, pos)
        
    print("\t----------\n\tS: {0}, I: {1}, R: {2}".format(S,I,R))

    return {"Suscceptable":S,"Infected":I, "Recovered":R}







def visualize2(network, title, pos):
    """
    Visualize the network given an array of positions.
    """
    print("-- Starting to Visualize --")
    MAP = True

    #if MAP:
    m = Basemap(
        #projection='cea',
        projection='merc',
        ellps='WGS84',
        llcrnrlat=20, urcrnrlat=50,
        llcrnrlon=-130, urcrnrlon=-90,
        resolution=None
        )

    pos = dict()

    for pos_node in network.nodes():
        # Normalize the lat and lon values
        x,y = m(float(network.node[pos_node]['lon']),
                float(network.node[pos_node]['lat']))

        pos[pos_node] = [x,y]

    #m.fillcontinents(color="brown")

    colors = []
    i_edge_colors = []
    d_edge_colors = []
    default = []
    infected = []
    for node in network.nodes():
        colors.append(network.node[node]["color"])
    for i,j in network.edges():
        color = network.node[i]["color"]
        alpha = 0.75
        if color == "#A0C8F0" or color == "#FF6F00" or color == "purple":
            color = "#A6A6A6"
            default.append((i,j))
            d_edge_colors.append(color)
        else:
            color = "#29A229"
            infected.append((i,j))
            i_edge_colors.append(color)

    #plt.figure(figsize=(10,10))
    #plt.figure(figsize=(7,7))

    # Fist pass - Gray lines
    nx.draw_networkx_edges(network,pos,edgelist=default,
            width=0.5,
            edge_color=d_edge_colors,
            alpha=0.5,
            arrows=False)
   
    # Second Pass - Colored lines
    nx.draw_networkx_edges(network,pos,edgelist=infected,
            width=0.5,
            edge_color=i_edge_colors,
            alpha=0.75,
            arrows=False)

    nx.draw_networkx_nodes(network,
            pos,
            linewidths=0.5,
            node_size=10,
            with_labels=False,
            node_color = colors)
    
    # Adjust the plot limits
    cut = 1.05
    xmax = cut * max(xx for xx,yy in pos.values())
    xmin =  min(xx for xx,yy in pos.values())
    xmin = xmin - (cut * xmin)


    ymax = cut * max(yy for xx,yy in pos.values())
    ymin = (cut) * min(yy for xx,yy in pos.values())
    ymin = ymin - (cut * ymin)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    #if MAP:
        # Draw the map
    m.bluemarble()
    plt.title=title

    plt.axis('off')

    number_files = str(len(os.listdir()))
    while len(number_files) < 3:
        number_files = "0" + number_files

    plt.savefig("infection-{0}.png".format(number_files),
                bbox_inches='tight', dpi=600 
            )
    plt.close()

if __name__ == "__main__":
    main()
