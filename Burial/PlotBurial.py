import matplotlib.pyplot as plt
import seaborn as sns
import sys
import numpy as np
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)
from Bio.PDB import NeighborSearch
from Bio.PDB import Selection


def plotHistogram(metric, positions, labels, ylabel, colors, markers = ['o','o','o','o','o','o','o','o','o','o','o',]):  # ['X','o','s','o','s','o','s','o','s','o','s','o','s']):
	sns.set()
	width = 2*len(positions) + 2
	fig = plt.figure(figsize=(width,6))
	ax = fig.add_subplot(111)
	ax.set_ylabel(ylabel, fontsize = 16)
	plt.tick_params(axis='both', which='major', labelsize=14)

	counter = 0
	values = []
	for these_positions in positions:
		label = labels[counter]
		color = colors[counter]
		marker = markers[counter]
		ys = []
		for pos in these_positions:
			ys.append(np.amin(metric[pos]))
		values.append(ys)
		counter +=1
		xs = 0.2 * np.random.random(len(ys)) - 0.1 + counter
		plt.scatter(xs, ys, c=color, marker=marker)      

	plt.boxplot(values,labels=labels)
	plt.xticks(rotation=-40)

	figname =  str("Category_%s_boxplot.pdf") %ylabel
	fig.savefig(figname, bbox_inches='tight')
	plt.tight_layout()
	plt.show()


def readSASAFile(filename, value="percent"):
	data = {}
	counter = 0
	infile = open(filename, 'r')
	for line in infile:
		if line[1] == '-':
			break
	
		counter +=1

		if counter <= 2:
			continue

		fields = line.split()	
		pos = fields[1]

		if value == "percent":
			data[pos] = float(fields[6])

		elif value == "classification":
			if len(fields) == 7:
				data[pos] = 0
			elif len(fields) == 8:
				if fields[-1] == 'i':
					data[pos] = -1
				elif fields[-1] == 'o':
					data[pos] = 1

		elif value == "total":
			data[pos] = float(fields[2])

		elif value == "sidechain":
			data[pos] = float(fields[5])
	
	return data
	

def calculateNeighbors(filename, radius):
	data = {}
	structure = parser.get_structure(filename.split(".pdb")[0], filename)
	atom_list = Selection.unfold_entities(structure, 'A') # A for atoms
	residue_list = Selection.unfold_entities(structure, 'R') # R for residues
	neighbor_search = NeighborSearch(atom_list)
	
	for residue in residue_list:
		resid = str(residue.get_id()[1])
		contacts = []
		for atom in residue.get_list():
			contacts.extend(neighbor_search.search(atom.get_coord(), radius, level = "A"))
		burial = len(contacts)/len(residue.get_list())
		data[resid] = burial
	
	return data


def condenseData(all_data, func=np.mean):
	output = {}	

	for ident in all_data[0].keys():
		vals = []
		for dataset in all_data:
			vals.append(dataset[ident])
		output[ident] = func(vals)
		
	return output
	

radius = float(sys.argv[4])

neighbors = []
for filename in sys.argv[2].split(','):
	neighbors.append(calculateNeighbors(filename, radius))
neighbor_count = condenseData(neighbors, func=np.amin)

SASA_classes = []
total_SASAs = []
sidechain_SASAs = []
percent_SASAs = []
for filename in sys.argv[3].split(','):
	SASA_classes.append(readSASAFile(filename, value='classification'))
	total_SASAs.append(readSASAFile(filename,value='total'))
	sidechain_SASAs.append(readSASAFile(filename,value='sidechain'))
	percent_SASAs.append(readSASAFile(filename,value='percent'))	

sidechain_SASA = condenseData(sidechain_SASAs, func=np.amax)

datafile = open(sys.argv[1])
all_positions = []
all_colors = []
all_labels = []
for line in datafile:
	line = line.strip()
	name, color, mapfilename, positions = line.split()
	positions = sorted(positions.split(','))
	position_list = []
	color_list = []
	label_list = []
	position_list.append(positions)
	all_positions.append(positions)
	label_list.append(name)
	all_labels.append(name)
	color_list.append(color)
	all_colors.append(color)
datafile.close()

plotHistogram(neighbor_count, all_positions, all_labels, "Min(Neighbors)", all_colors)
plotHistogram(sidechain_SASA, all_positions, all_labels, "Max(Sidechain_SASA)", all_colors)

