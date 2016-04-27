import sys
import re

# python convertAirlinePassenger.py
# ./Data/182980864_T_T100D_SEGMENT_ALL_CARRIER2.csv ./Data/airportsMin.csv


def main():

    args = sys.argv

    dataFile = args[1]
    airportsFile = args[2]

    data = open(dataFile, 'r')
    airports = open(airportsFile, 'r')
    outputFile = open('./Data/airlineRoutesPassengerData.csv','w')
    approvedAirports = dict()

    # load approved Airports
    for line in airports:
        segmentedLine = line.split(",")
        approvedAirports[segmentedLine[2]] = segmentedLine[0]
    airports.close()

    foundRoutes = dict()

    # Iterate through each line in airline flight data file
    data.readline() # throw away header
    for line in data:
        segmentedLine = re.split('"?,?"?',line) # split on the separators
        # Look for IATA code in approvedAirports list
        if segmentedLine[4] in approvedAirports and segmentedLine[6] in \
                approvedAirports and float(segmentedLine[0]) > 0:
            # Org IATA code, Org Airport ID, Des IATA code,
            # Des Airport ID, Year, # Passengers
            if (segmentedLine[4],segmentedLine[6]) not in foundRoutes:
                foundRoutes[(segmentedLine[4],segmentedLine[6])] = \
                float(segmentedLine[0])
            else:
                foundRoutes[(segmentedLine[4],segmentedLine[6])] += \
                float(segmentedLine[0])

    for org, des in foundRoutes:
        if foundRoutes[(org,des)] > 365:
            tempList = [org, approvedAirports[org], des, approvedAirports[des],
                        str(foundRoutes[(org,des)]), "\n"]
            outputFile.write(",".join(tempList))

    data.close()
    outputFile.close()


if __name__ == '__main__':
    main()