import sys
import re

def main():

    args = sys.argv

    dataFile = args[1]

    data = open(dataFile, 'r')
    outputFile = open('./Data/airportsMin.csv','w')

    approvedAirports = ['ATL','LAX','ORD','DFW','JFK','SFO','MIA','CLT','LAS',
                        'PHX','IAH','MCO','EWR','MSP','BOS','PHL','LGA','FLL',
                        'BWI','IAD','MDW','DCA','HNL','SAN','TPA']
    count = 0
    # Iterate through each line in original file
    for line in data:
        segmentedLine = re.split('"?,?"?',line) # split on the separators
        # Look for IATA code in approvedAirports list
        for item in segmentedLine:
            if item in approvedAirports:
                # Save Airport Name, IATA code, lat, long
                outputList = [segmentedLine[1],segmentedLine[4],
                              segmentedLine[6],segmentedLine[7],'\n']
                outputFile.write(",".join(outputList))

    data.close()
    outputFile.close()


if __name__ == '__main__':
    main()