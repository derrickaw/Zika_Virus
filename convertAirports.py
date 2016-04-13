import sys
import re

def main():

    args = sys.argv

    dataFile = args[1]

    data = open(dataFile, 'r')
    outputFile = open('./Data/airportsMin.csv','w')

    approvedAirports = {'ATL':456002,'LAX':3928864,'ORD':2084044,
                         'DFW':1281047,'JFK':5618852,'SFO':852469,'MIA':430332,
                         'CLT':809958,'LAS':613599,'PHX':1537058,'IAH':2239558,
                         'MCO':262372,'EWR':280579,'MSP':407207,'BOS':655884,
                         'PHL':1560297,'LGA':2872227,'FLL':176013,'BWI':622793,
                         'IAD':335210,'MDW':638345,'DCA':323683,'HNL':350399,
                         'SAN':1381069,'TPA':358699}

    foundAirports = list()

    # Iterate through each line in original file
    for line in data:
        segmentedLine = re.split('"?,?"?',line) # split on the separators
        # Look for IATA code in approvedAirports list
        for item in segmentedLine:
            if item in approvedAirports:
                # Save ID, Airport Name, IATA code, lat, long
                outputList = [segmentedLine[1],segmentedLine[4],
                              segmentedLine[6],segmentedLine[7],
                              str(approvedAirports[item]),'\n']
                foundAirports.append(outputList)
                #outputFile.write(",".join(outputList))

    # Sort based on population
    foundAirports = sorted(foundAirports, key=lambda airport: int(airport[4]),
                           reverse=True)

    # Write to file
    count = 1
    for item in foundAirports:
        item.insert(0,str(count))
        count += 1
        outputFile.write(",".join(item))

    data.close()
    outputFile.close()


if __name__ == '__main__':
    main()