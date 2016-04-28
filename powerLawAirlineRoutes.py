import sys
import re
import matplotlib.pyplot as plt
import powerlaw
import numpy as np


# python powerLawAirlineRoutes.py ./Data/airports.dat
# ./Data/182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv


def main():

    args = sys.argv

    airportsDataFile = args[1] # airports.dat
    routeDataFile    = args[2] # 182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv

    airportData = open(airportsDataFile, 'r')
    routeData   = open(routeDataFile, 'r')
    outputFile = open('./Data/powerLawData.csv','w')

    domesticAirports = dict()

    # Iterate through each line in original airport file to find all domestic
    # airports (nodes)
    for line in airportData:
        segmentedLine = re.split('"?,?"?',line) # split on the separators
        # Save IATA code for all domestic airports
        if segmentedLine[3] == "United States":
            if segmentedLine[3] not in domesticAirports:
                # Store number of routes, number of passengers flowing out,
                # number of connected airports
                domesticAirports[segmentedLine[4]] = [0,0,dict()]


    routeData.readline() # throw away header
    for line in routeData:
        segmentedLine = re.split('"?,?"?',line) # split on the separators
        # Look for IATA code in approvedAirports list

        # Check to make sure the from and to airports are both domestic and that
        # at least one person traveled between the cities
        if segmentedLine[4] in domesticAirports and segmentedLine[6] in \
                domesticAirports and float(segmentedLine[0]) > 0:

            if segmentedLine[6] not in domesticAirports[segmentedLine[4]][2]:
                domesticAirports[segmentedLine[4]][2][segmentedLine[6]] = 1
                domesticAirports[segmentedLine[4]][0] += 1
            domesticAirports[segmentedLine[4]][1] += float(segmentedLine[0])

    # Find max num
    maxNum = 0
    for key, values in domesticAirports.items():
        if values[0] > maxNum:
            maxNum = values[0]

    domesticAirportsList = [0] * (maxNum + 1)
    # domesticAirportsList2 = dict()
    for key, values in domesticAirports.items():
        domesticAirportsList[values[0]] += 1
    #     for i in range(0,1000000000,1000000):
    #         if (values[1] <= i):
    #             if i not in domesticAirportsList2:
    #                 domesticAirportsList2[i] = 0
    #             domesticAirportsList2[i] += 1
    #             break



    #print (domesticAirports)
    # print (domesticAirportsList)
    #print (domesticAirportsList2)

    # keys = domesticAirportsList2.keys()
    # values = domesticAirportsList2.values()

    #plt.loglog(list(keys), list(values), '*', basex=10)
    #plt.show()
    #print (keys, values)

    xData = list()
    yData = list()
    for i in range(1,len(domesticAirportsList)):
        if domesticAirportsList[i] > 0:
            xData.append(i)
            yData.append(domesticAirportsList[i])

    lines = []
    for i in range(len(xData)):
        lines.append(str(xData[i]) + ", " + str(yData[i]) + "\n")
    outputFile.writelines(lines)

    # xTravelData = list()
    # yTravelData = list()

    #for i in range(1,len(domesticAirports))



    # for i in range(1,maxNum+1):
    #     xData.append(i)

    print ("xdata",xData)
    print ("ydata",yData)
    plt.loglog(xData, yData, '*', basex=10)
    plt.xlabel("Route Degree")
    plt.ylabel("Number of Routes")
    plt.grid(True)
    #plt.plot(xData,domesticAirportsList[1:])
    plt.show()

    # plt.plot(xData, np.poly1d(np.polyfit(xData, yData, 1))(xData))
    # plt.show()




    #fit = powerlaw.Fit(domesticAirportsList)
    #print (fit.power_law.alpha)
    #print (fit.power_law.xmin)

    #fit.power_law.plot_pdf( color= 'b',linestyle='--',label='fit ccdf')
    #fit.plot_pdf( color= 'b')









    # Write to file
    # count = 1
    # for item in foundAirports:
    #     item.insert(0,str(count))
    #     count += 1
    #     outputFile.write(",".join(item))

    airportData.close()
    routeData.close()
    outputFile.close()


if __name__ == '__main__':
    main()