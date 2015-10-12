#!/usr/bin/python
import os
from Bio import SeqIO
import numpy
import matplotlib.pyplot as plt

data_dir = '../data/'
output_dir = '../output/'
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
	plot_bar(dictionary)


# Plot a bar graph for the number of each amino acid in the proteome sequence.
def plot_bar(dictionary, location=None):
	if location == None:
		figs = str(output_dir) + '/plots/'
	else:
		figs = str(location) + '/'
	path_to_dir(figs)
	filename = str(figs) + 'Distribution.png'
	plt.figure().canvas.set_window_title('Distribution')
	print(plt.hist(numpy.asarray(dictionary.values()), bins=100, log=True))
	plt.savefig(filename)
	plt.close()


def path_to_dir(out):
	# Create the specified folder if it does not already exist.
	if not os.path.exists(out) and not out == '':
		os.makedirs(out)

	# If no directory is specified to store the data, store it in the current directory of the user.
	if out == '':
		home = os.path.expanduser('~')
		out = './' + str(uniprot_data.species) + '/'
		os.makedirs(out)


if __name__ == '__main__':
	for file in files:
		if file.endswith('.fasta'):
			parse(file)