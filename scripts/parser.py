#!/usr/bin/python
import os
from Bio import SeqIO

data_dir = '../data/'
files = os.listdir(data_dir)

def parse(file):
	print file
	dictionary = {}
	for data in SeqIO.parse(data_dir + str(file), 'fasta'):
		# print data.id
		# print data.seq
		dictionary[data.id] = len(data)
		# lis.append(len(data))
	print max(dictionary.values())
	print dictionary.keys()[int(dictionary.values().index(max(dictionary.values())))]


if __name__ == '__main__':
	for file in files:
		if file.endswith('.fasta'):
			parse(file)