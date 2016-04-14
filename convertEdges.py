import sys
import re

def main():

    args = sys.argv

    dataFile = args[1]
    dataFile1 = args[2]

    rankDict = {}
    data1 = open(dataFile1, 'r')
    for line in data1:
        segmentedLine = line.split(',') # split on the separators
        print(segmentedLine)
        rankDict[segmentedLine[2]] = segmentedLine[0]
    #print(len(rankDict))
    data1.close()

    data = open(dataFile, 'r')
    outputFile = open('./Data/airportsEdges.csv','w')

    approvedAirports = ['ATL','LAX','ORD','DFW','JFK','SFO','MIA','CLT','LAS',
                        'PHX','IAH','MCO','EWR','MSP','BOS','PHL','LGA','FLL',
                        'BWI','IAD','MDW','DCA','HNL','SAN','TPA']
    #count = 0
    # Iterate through each line in original file
    outputList = set()
    for line in data:
        segmentedLine = line.split(',') # split on the separators

        #print(segmentedLine)
        # Look for IATA code in approvedAirports list
        if segmentedLine[2] in approvedAirports and segmentedLine[4] in approvedAirports:
            #print(segmentedLine[2], ":", segmentedLine[4])
            #Save Airport Name, IATA code, lat, long
            outputList.add((segmentedLine[2],segmentedLine[4]))

    outputList = list(outputList)
    for item in outputList:
        if (item[1], item[0]) in outputList:
            outputList.remove((item[1], item[0]))
    #print(outputList)
    #print(len(outputList))

    for item in outputList:
        outputFile.write(str(item[0]) + ',' + str(rankDict[item[0]]) + ',' + str(item[1]) + ',' + str(rankDict[item[1]]) + "\n")

    data.close()

    outputFile.close()


if __name__ == '__main__':
    main()