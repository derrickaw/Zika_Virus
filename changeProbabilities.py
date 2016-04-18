import sys

def main():

    args = sys.argv

    dataFile = args[1]

    data = open(dataFile, 'r')
    outputFile = open('./Data/cityProbabilities.csv','w')

    outputList = list()
    # Iterate through each line in original file
    for line in data:
        segmentedLine = line.split(',') # split on the separators
        segmentedLine[12], segmentedLine[10] = segmentedLine[10], segmentedLine[12]
        segmentedLine[11], segmentedLine[8] = segmentedLine[8], segmentedLine[11]
        segmentedLine[10], segmentedLine[6] = segmentedLine[6], segmentedLine[10]

        print("SEG LINE", segmentedLine)
        outputList.append(segmentedLine)

    for item in outputList:
        outputFile.write(str(item) + "\n")

    data.close()

    outputFile.close()


if __name__ == '__main__':
    main()