#!/usr/bin/python
from __future__ import division
import os, json, re
import numpy
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import bipartite
from Bio import SeqIO
from bitmap import BitMap

data_dir = '../data/'
output_dir = '../output/'
files = os.listdir(data_dir)

def parse(file):
	protein_gene_map = {}
	print 'Processing file: ' + str(file)
	dictionary = {}
	seq_len = {}
	for data in SeqIO.parse(data_dir + str(file), 'fasta'):
		# print data
		# print data.description
		# print data.id
		# print data.seq
		if 'GN=' in data.description:
			gene_name = (data.description).split('GN=')[1].split()[0]
			# print 'Protein: ' + str(data.id)
			# print 'Gene: ' + gene_name
		# gene_name = (data.description).split('=')[2].split('PE')[0].split()[0]
		# print re.split(r'[GN= PE]', data.description)
		protein_gene_map[data.id] = gene_name
		seq_len[data.id] = len(data.seq)
		dictionary[data.id] = {}
		for i in range(0, len(data.seq)):
			if data.seq[i] == 'X':
				pass
			else:
				if data.seq[i] in dictionary[data.id].keys():
					count += 1
					dictionary[data.id][str(data.seq[i])] = float(count/len(data.seq))
				else:
					count = 1
					dictionary[data.id][str(data.seq[i])] = float(1/len(data.seq))
		# print dictionary
		# raw_input('Enter')
	print 'Longest seq: ' + str(max(seq_len.values()))
	print seq_len.keys()[int(seq_len.values().index(max(seq_len.values())))]
	print 'Shortest seq: ' + str(min(seq_len.values()))
	print seq_len.keys()[int(seq_len.values().index(min(seq_len.values())))]
	print 'Average seq length: ' + str(numpy.mean(seq_len.values()))
	jsonify(dictionary, './../output/amino_acid_count.py', 'amino_acid_count')
	jsonify(protein_gene_map, './../output/protein_gene_map.py', 'protein_gene_map')
	# plot_bar(dictionary)
	# export_gexf(protein_gene_map)


def export_gexf(protein_gene_map):
	G = nx.Graph()
	count = 0
	for protein, gene in protein_gene_map.items():
		G.add_node(protein)
		G.add_nodes_from(gene)
		G.add_edge(protein, gene)
		count += 1

		c = bipartite.color(G)
		nx.set_node_attributes(G, 'bipartite', c)
		nx.write_gexf(G,"protein_gene_map.gexf")


def protein_bitmap():
	# bm = BitMap(32)
	pass

# Plot a bar graph for the number of each amino acid in the proteome sequence.
def plot_bar(dictionary, location=None):
	if location == None:
		figs = str(output_dir) + '/plots/'
	else:
		figs = str(location) + '/'
	path_to_dir(figs)
	filename = str(figs) + 'Distribution.png'
	plt.figure().canvas.set_window_title('Distribution')
	plt.hist(numpy.asarray(dictionary.values()), bins=100, log=True)
	plt.savefig(filename)
	plt.close()


def jsonify(dictionary, filename, text='None'):
	a = json.dumps(dictionary, sort_keys=True, indent=4, separators=(',', ': '))
	with open(str(filename), 'w') as outfile:
		if text == 'None':
			outfile.write(a)
		else:
			outfile.write(text + ' = ')
			outfile.write(a)


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
		if file.endswith('.fasta') and '5640' in file:
			parse(file)