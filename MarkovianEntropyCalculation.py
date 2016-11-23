# MarkovianEntropyoCalculation.py
#
# Input: csv file where each row is a discrete time series, from TimeSeriesDiscretization.py
# Output: csv file where each row is a Markovian Entropy value for the given discrete time series
# Output: csv file where each row is the transition probability matrix computed from the given
#         discrete time series, written in column-major order

import csv
import copy
import math
import numpy as np

# Parametrize n
n = 2 # Currently compatible up to 26
k = 1

# Set up to import calcium activity csv
# Alter filename conventions as necessary
inFilename = 'test_discrete_n-' + str(n) + '.csv'
outFilename = 'test'

matrixList = []

def writeMatrices(self, outFilename):
	with open(outFilename, 'wb') as csvfile:
		mWriter = csv.writer(csvfile, delimiter = ',')
		for z in range(0, len(self)):
			mWriter.writerow(self[z])


def calculateEntropy(self): # input is one cell
	m = np.zeros([n]*(k + 1))

	m = list()
	for i in range(0, n):
		m.append([0] * n**k)
	# Note that the m matrix is indexed [column][row]!

	alphabet = dict([('A',0), ('B',1), ('C',2), ('D',3), ('E',4), ('F',5), ('G',6), ('H',7), \
	 ('I',8), ('J',9), ('K',10), ('L',11), ('M',12), ('N',13), ('O',14), ('P',15), \
	 ('Q',16), ('R',17), ('S',18), ('T',19), ('U',20), ('V',21), ('W',22), ('X',23), \
	 ('Y',24), ('Z',25)]) 

	for t in range(k, len(self)):

		temp = list()
		for j in range(0, k + 1):
			temp.append(alphabet[self[t-j]])

		col = temp.pop(0)
		temp.reverse()

		ndx = 0
		for z in range(0, len(temp)):
			ndx += (n**(len(temp)-z-1))*temp[z]

		m[col][ndx] += 1

	for row in range(0, n**k):
		s = 0
		for col in range(0, n):
			s += m[col][row]
		
		if s != 0:
			for col in range(0, n):
				c = copy.deepcopy(m[col][row])
				p = round(c * s**(-1), 4) 
				m[col][row] = p

	e = list()
	for row in range(0, n**k):

		temp = 0
		for col in range(0, n):
			if m[col][row] != 0:
				temp += m[col][row] * math.log(m[col][row], 2)
		e.append(-1 * temp)

	mFlat = [x for sublist in m for x in sublist] # column major order!

	matrixList.append(mFlat)

	return (sum(e)/((n**k)*math.log(n,2)))


def writeEntropies(self, outFilename):
	with open(outFilename, 'wb') as csvfile:
		eWriter = csv.writer(csvfile, delimiter = ',')
		for z in range(0, len(self)):
			eWriter.writerow([self[z]])


def main():
	cellValues = []
	with open(inFilename, 'rU') as csvfile: 
		cellreader = csv.reader(csvfile, delimiter = ' ')
		cells = []
		for row in cellreader:
			cells.append(row)

	cellData = []
	for cell in cells:
		newCell = cell[0].split(',')
		for i in cell:
			listOfValues = i.split(',')
			newCell = []

			for j in listOfValues: 
				newCell.append(j)
			cellValues.append(newCell)
		entropy = calculateEntropy(newCell)
		cellData.append(entropy)

	# write output
	outFilename_entropies = outFilename + '_entropies_n-' + str(n) + '_k-' + str(k) + '.csv'

	outFilename_matrices = outFilename + '_matrices_n-' + str(n) + '_k-' + str(k) + '.csv'

	writeEntropies(cellData, outFilename_entropies)

	writeMatrices(matrixList, outFilename_matrices)

main()


