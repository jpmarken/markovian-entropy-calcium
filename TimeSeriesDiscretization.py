# TimeSeriesDiscretization.py
#
# Input: csv file where each row is a calcium activity time series with numerical activity values
# Output: csv file where each row is a discrete time series with states delineated by letters

import csv
import copy
import array
import numpy as np

# Parametrize n
n = 2 # Currently compatible up to 26

# Set up to import calcium activity csv
# Alter filename conventions as necessary
inFilename = 'test_flx.csv'
outFilename = 'test_discrete_n-' + str(n) + '.csv'


def calculateQuantiles(self): 
	breakPoints = []
	for ii in range(0,n-1):
		breakPoints.append(np.percentile(self, (100/n)*(ii+1)))

	return breakPoints


def discretize(self, breakPoints): 
	discreteCellValues = []

	alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

	for i in self:
		done = 0
		for j in range(0,n-1):
			if (done == 0) and (i < breakPoints[j]):
				discreteCellValues.append(alphabet[j])
				done = 1
		if done == 0:
			discreteCellValues.append(alphabet[n-1])

	return discreteCellValues


def discretizeEveryCell(self): 
	allDiscretizedCells = []
	for cell in self:
		temporaryCellValues = copy.copy(cell)
		boundaryValues = calculateQuantiles(temporaryCellValues)
		discretizedCell = discretize(cell, boundaryValues)
		allDiscretizedCells.append(discretizedCell)

	return allDiscretizedCells


def writeOutput(self): 
	with open(outFilename, 'wb') as csvfile:
		discreteWriter = csv.writer(csvfile, delimiter = ',')
		discreteWriter.writerows(self)


def main():
	cellValues = []
	with open(inFilename, 'rU') as csvfile: 
		cellreader = csv.reader(csvfile, delimiter = ' ')
		cells = []
		for row in cellreader:
			cells.append(row)

	for cell in cells:
		for i in cell:
			listOfValues = i.split(',')
			newCell = []

			for j in listOfValues: 
				if j is not '': 
					if j != str(0):
						j = float(j)
						newCell.append(j)
			cellValues.append(newCell)

	output = discretizeEveryCell(cellValues)
	writeOutput(output)

main()