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
import matplotlib.animation as animation
import numpy as np


# GLOBAL
MAP = False
TAU = 4      # Zika virus "Ro", DEFAULT
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
START = 1
DAYS_IN_YEAR = 365
DAYS_IN_MONTH = 31
SIMULATION_LENGTH = DAYS_IN_YEAR
STAT = False
INCUBATION = 3  # 3-12 days , DEFAULT
TO_RECOVER = 7  # DEFAULT
RUN_ALL = False
VACCINATE_PERC = 1 - (1/TAU) # DEFAULT
VACCINATE = False
SCREEN_PERC = 0
MONTH_MAP = {1: "January", 2: "February", 3: "March", 4: "April", 5: "May",
             6:"June", 7: "July", 8: "August", 9: "September", 10: "October",
             11: "November", 12: "December"}

def main():
    """
    Primary function that initiates network creation and handles execution of
    infection simulations.

    Args:
        argv: A list of command-line arguments passed to the application.

    Returns:
        Void

    """
    global MAP, CITY_TO_INFECT, START, SIMULATION_LENGTH, DATE_TO_INFECT, STAT,\
        timeStepsTracker, INCUBATION, RUN_ALL, VACCINATE_PERC, VACCINATE, TAU, \
        SCREEN_PERC

    # Determine the parameters of the current simulation.
    opts, args = getopt.getopt(sys.argv[1:], "msav", ["c=", "d=", "start=",
                                                       "days=","tau=", "inc=",
                                                       "vac=", "screen="])

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
        # Run all cities through infection for every month to find the worst
        # set of city - infection month pairs for number of infections and
        # number of nodes infected
        elif opt == "-a":
            RUN_ALL = True
        # vaccinate using default percent based on TAU
        elif opt == "-v":
            VACCINATE = True
        # infect a particular city
        elif opt == "--c":
            CITY_TO_INFECT = par
        # date to infect
        elif opt == "--d":
            DATE_TO_INFECT = int(par)
        # Beginning part of simulation
        elif opt == "--start":
            START = int(par)
        # Number of days to run simulation, if not entered, default is 365
        elif opt == "--days":
            SIMULATION_LENGTH = int(par)
        # A different tau for mosquito dynamics and preset vaccination rate
        # based on TAU, if user wants to change it, they must use the --v option
        elif opt == "--tau":
            TAU = float(par)
            VACCINATE_PERC = 1 - (1/TAU)
        # Set the incubation period before symptoms show and the infected
        # person can propagate the disease
        elif opt == "--inc":
            INCUBATION = float(par)
        # Change vaccination percentage and set model to vaccinate
        elif opt == "--vac":
            VACCINATE_PERC = float(par)
            VACCINATE = True
        # Percent of passengers to screen from normal population proportion
        # dynamics
        elif opt == "--screen":
            SCREEN_PERC = float(par)


    # Create the network using the command arguments
    network = create_network(AIRPORT_DATA, ROUTE_DATA, MOSQUITO_CURVES)

    # Setup Global SIVR stat tracker
    setupGlobalSIVR()

    infectionAllStats = dict()

    # Entire network run
    if RUN_ALL:
        for airport in approvedAirports:
            # if airport == "MIA":
            CITY_TO_INFECT = airport
            print (airport)
            infectionAllStats[airport] = list()
            for i in range(1,DAYS_IN_YEAR,DAYS_IN_MONTH): # START + SIMULATION_LENGTH+1
                #print (i)
                networkCopy = network.copy()
                setupGlobalSIVR()
                DATE_TO_INFECT = i
                infectionAllStats[airport].append(0)
                # Run infection simulation
                for j in range(i,DAYS_IN_YEAR+i): # START + SIMULATION_LENGTH+1+i
                    j %= SIMULATION_LENGTH
                    if j % INCUBATION == 0 or j == DATE_TO_INFECT:
                        # timeStepsTracker.append(i)
                        infection(networkCopy, j)
                    # for node in networkCopy.nodes_iter(networkCopy):
                    #     if node[1]["IATA"] == "MIA" and airport == "MIA":
                    #         print (node)
                for node in networkCopy.nodes_iter(networkCopy):
                    # print (node[1]["R"])
                    infectionAllStats[airport][-1] += node[1]["R"]
                    #if airport == "MIA":
                    #    print (network.nodes(networkCopy))
                # print
        #print (infectionAllStats)
        for airport in infectionAllStats:
            print(airport, infectionAllStats[airport])
    else:

        for i in range(START,START+SIMULATION_LENGTH):
            if i % INCUBATION == 0 or i == DATE_TO_INFECT:
                timeStepsTracker.append(i)
                infection(network, i)

    # Print # of infected timeline
    # for airport in approvedAirports:
    #     print (airport, I[airport])

    # Print # of recovered per airport
    for airport in approvedAirports:
        print(airport, R[airport][-1])




    # Visualize network
    if MAP:
        #visualize(network)
        theIDic = updateIDic(network)
        print("Infected Dict Ratios: ", theIDic)
        print("ATL Ratio List : ", theIDic['ATL'])
        print("Length of ATL List", len(theIDic['ATL']))
        for i in range(len(theIDic[CITY_TO_INFECT])):
            month = i + START//DAYS_IN_MONTH % 12
            print("MON", month)
            updatedVisualize(network, theIDic, i, month)
    # Stats of infection
    if STAT:
        for node in I:
            # only for city to infect; can change to others by commenting out
            if node == CITY_TO_INFECT:

                i, = plt.plot(timeStepsTracker, I[node],label="I")
                s, = plt.plot(timeStepsTracker, S[node],label='S')
                r, = plt.plot(timeStepsTracker, R[node],label='R')
                v, = plt.plot(timeStepsTracker, V[node],label='V')
                plt.legend(handles=[i,s,r,v], loc = 'best')

        plt.title(CITY_TO_INFECT + " Infection Dynamics")
        plt.xlabel('Days of Year')
        plt.ylabel('People')
        plt.xlim(START,START+SIMULATION_LENGTH)
        plt.show()


# def reorderStats():
#     global timeStepsTracker, I, S, R, V
#
#     foundSpot = 0
#     for i in range(1,len(timeStepsTracker)):
#         if timeStepsTracker[i-1] > timeStepsTracker[i]:
#             foundSpot = i
#             break
#         print(foundSpot)
#     timeStepsTracker = timeStepsTracker[foundSpot:] + \
#                        timeStepsTracker[:foundSpot]
#     for airport in approvedAirports:
#         I[airport] = I[airport][foundSpot:] + I[airport][:foundSpot]
#         S[airport] = S[airport][foundSpot:] + S[airport][:foundSpot]
#         R[airport] = R[airport][foundSpot:] + R[airport][:foundSpot]
#         V[airport] = V[airport][foundSpot:] + V[airport][:foundSpot]




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
            population = int(entries[5])
            vaccinated = 0

            # if VACCINATE:
            #     vaccinated = math.ceil(population * VACCINATE_PERC)
            #     population -= vaccinated

            G.add_node(int(entries[0]),
                       name=entries[1],
                       IATA=entries[2],
                       lat=entries[3],
                       lon=entries[4],
                       pop=int(entries[5]),
                       pos=(float(entries[3]),float(entries[4])),
                       I=[0,0],  # First num - Total I, second num new I
                       Iair=0,
                       S=population,       #createHumans(int(entries[5])),
                       V=vaccinated,
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
    # print("\tCalculating clustering coefficents",end="")
    # cluster_network = nx.Graph(G)
    # lcluster = nx.clustering(cluster_network)
    # for i,j in G.edges():
    #     cluster_sum = lcluster[i] + lcluster[j]
    #     G[i][j]['cluster'] = cluster_sum
    # print("\t\t\t[Done]")

    return G






def infection(input_network, timeStep):
    global I,S,V,R
    # mosCurve = [0,0,0,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0.33,0] # phx,
    approxMonth = timeStep // DAYS_IN_MONTH % 12


    # Spread disease to other Airports
    currentNodes = input_network.node
    for node in input_network.nodes_iter(input_network):
        for key in routeInfo:
            if node[1]["IATA"] == key[2]:
                nodeDetails = currentNodes[int(key[1])]
                # if node[1]["IATA"] == "ATL":
                #     print (int(nodeDetails["I"][0] / \
                #                        nodeDetails["pop"] * \
                #                        routeInfo[key] / 365))

                node[1]["Iair"] += math.ceil((int(nodeDetails["I"][0] /
                                   nodeDetails["pop"] *
                                   routeInfo[key] / DAYS_IN_YEAR)) *
                                   (1-SCREEN_PERC))


    #print (input_network.node)

    # Infection simulation at hubs
    for node in input_network.nodes_iter(input_network):
        # if node[1]["IATA"] == "MIA":
        #    print(node[1]["I"])

        #  Record stats
        I[node[1]["IATA"]].append(node[1]["I"][0])
        S[node[1]["IATA"]].append(node[1]["S"])
        V[node[1]["IATA"]].append(node[1]["V"])
        R[node[1]["IATA"]].append(node[1]["R"])


        # Check for recovery
        if node[1]["I"][0] > 0:
            if timeStep - node[1]["I"][2][0] >= (TO_RECOVER + INCUBATION):
                group = node[1]["I"].pop(2)
                node[1]["I"][0] -= group[1]
                node[1]["R"]    += group[1]
                # If the recovered amount leaving matches the last one, then
                # clear
                if group[1] == node[1]["I"][1]:
                    node[1]["I"][1] = 0

        # Vaccinate people
        if VACCINATE:
            if node[1]["I"][0] > 0 and node[1]["S"] > 0:
                num_vaccinate = min(math.ceil(VACCINATE_PERC * node[1]["I"][0] *
                                              node[1]["S"] / (node[1]["I"][0] +
                                              node[1]["S"])), node[1]["S"])
                node[1]["V"] += num_vaccinate
                node[1]["S"] -= num_vaccinate
                # print(node[1]["IATA"], num_vaccinate)


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



def updateIDic(network):
    popIDict = {}
    print("ATL's I :", I['ATL'])
    # divide all I by the population
    for node in network.nodes_iter(network):
        popIDict[node[1]["IATA"]] = [x / node[1]["pop"] for x in I[node[1]["IATA"]]]
    print("ATL's IDIC :", popIDict['ATL'])
    print("POPULATION IDict :", popIDict)
    updatedIDict = {}
    for key in popIDict:
        updatedIDict[key] = [0] * 12
    print("BEFORE UPDATED IDIC:",updatedIDict)

    timeStepMonth = list()
    for i in timeStepsTracker:
        if ((i // DAYS_IN_MONTH) % 12 == 0):
            #print(12)
            timeStepMonth.append(12)
        else:
            #print((i // DAYS_IN_MONTH) % 12)
            timeStepMonth.append((i // DAYS_IN_MONTH) % 12)
    print("TIME STEP MONTH:" , timeStepMonth)

    for key in popIDict:
        for i in range(len(timeStepMonth)):
            if (updatedIDict[key][(timeStepMonth[i]-1)] < popIDict[key][i]):
                updatedIDict[key][(timeStepMonth[i]-1)] = popIDict[key][i]

    print("AFTER BUT STILL NOT FINAL updatedIDict", updatedIDict)
    inital = timeStepMonth[0]
    for key, value in updatedIDict.items():
        value = value[inital-1:] + value[:inital-1]
        updatedIDict[key] = value

    print("FINAL updatedIDict", updatedIDict)
    # arr = (np.linspace(0, len(popIDict[CITY_TO_INFECT]),
                       #math.ceil((SIMULATION_LENGTH)/DAYS_IN_MONTH),
                       #endpoint=True, dtype=int))
    #arr[-1] -= 1
    #print(arr)
    #for key,value in popIDict.items():
        #updatedIDict[key] = [value[i] for i in arr]

    return updatedIDict

def getColor(value):
    if value <= .333:
        return "co"
    elif value <= .667:
        return "yo"
    else:
        return "ro"

def updatedVisualize(network, IDic, position, month):
    print("-- Starting to Visualize [", position+1, "] --")

    #updatedIDic = updateIDic(network)
    #print("UPDATED DIC", updatedIDic)

    map = Basemap(
        projection='merc',
        ellps='WGS84',
        llcrnrlon=-160,urcrnrlon=-60,llcrnrlat=10,urcrnrlat=80,
        resolution="l"
        )
    #map.drawmapboundary("aqua")
    #map.fillcontinents('#555555')
    #map.drawlsmask(land_color='green',ocean_color='aqua',lakes=True)
    map.bluemarble()

    pos = dict()

    for pos_node in network.nodes():
        # Normalize the lat and lon values
        x,y = map(float(network.node[pos_node]['lon']),
                float(network.node[pos_node]['lat']))
        #print("x,y", float(network.node[pos_node]['lon']),float(network.node[pos_node]['lat']))
        pos[pos_node] = [x,y]

    #print("POS", network.nodes())
    # First pass - edges
    nx.draw_networkx_edges(network,pos,edgelist=network.edges(),
            width=1,
            edge_color="gray",
            alpha=0.5,
            arrows=False)

    #Node/vertices colors
    nx.draw_networkx_nodes(network,
            pos,
            linewidths=1,
            node_size=30,
            with_labels=False,
            node_color = "white")


    for pos_node in network.nodes():
        # Normalize the lat and lon values
        x,y = map(float(network.node[pos_node]['lon']),
                float(network.node[pos_node]['lat']))
        #msize = math.ceil(network.node[pos_node]['pop'] * .00001)
        mcolor = getColor(IDic[network.node[pos_node]['IATA']][position])
        print("AIRPORT ATTRIBUTES:", network.node[pos_node]['IATA'], IDic[network.node[pos_node]['IATA']], mcolor)
        map.plot(x, y, mcolor, markersize=7)
        #print(updatedIDic[network.node[pos_node]['IATA']][2])


    #Adjust the plot limits
    cut = 1.05
    xmax = cut * max(xx for xx,yy in pos.values())
    xmin =  min(xx for xx,yy in pos.values())
    xmin = xmin - (cut * xmin)


    ymax = cut * max(yy for xx,yy in pos.values())
    ymin = (cut) * min(yy for xx,yy in pos.values())
    ymin = ymin - (cut * ymin)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    title_string = "Zika Infection for the Month of " + MONTH_MAP[month]
    savefigStr = "infection" + str(position) + ".png"
    print(str(savefigStr))
    #plt.savefig(savefigStr)
    plt.title(title_string)
    plt.axis('off')
    plt.show()
    plt.close()




def visualize(network):
    print("-- Starting to Visualize --")

    m = Basemap(
        projection='merc',
        ellps='WGS84',
        llcrnrlon=-160,urcrnrlon=-60,llcrnrlat=10,urcrnrlat=80,
        resolution="l"
        )

    pos = dict()
    labels = list()
    for pos_node in network.nodes():
        # Normalize the lat and lon values
        x,y = m(float(network.node[pos_node]['lon']),
                float(network.node[pos_node]['lat']))

        pos[pos_node] = [x,y]

    #m.drawmapboundary("aqua")
    #m.fillcontinents('#555555')
    #m.drawlsmask(land_color='green',ocean_color='aqua',lakes=True)
    m.bluemarble()

    # First pass - Green lines
    nx.draw_networkx_edges(network,pos,edgelist=network.edges(),
            width=1,
            edge_color="orange",
            alpha=0.5,
            arrows=False)

    nx.draw_networkx_nodes(network,
            pos,
            linewidths=1,
            node_size=40,
            with_labels=False,
            node_color = "white")

    # nx.draw_networkx_labels(network,pos,labels)


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



if __name__ == "__main__":
    main()
