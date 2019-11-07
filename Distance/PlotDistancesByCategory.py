from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import copy
from numpy import random

groupfilename = sys.argv[1]
structurename = sys.argv[2]

# create parser
parser = PDBParser()
# read structure from file
structure = parser.get_structure('Closed', structurename)
# store key locations in the DHFR structure
model = structure[0]
chain = model['A']
m20 = chain[20]['CA']
sheet = chain[112]['CA']
globular = chain[41]['CA']
ligand = model['X']
for residue in ligand:
	if "N" in residue.get_id()[0]:
		hydride = residue['H4']
		adenosine = residue['C18']


hydride_distances = {}
hydride_distance_list = []
adenosine_distances = {}
adenosine_distance_list = []
m20_distances = {}
m20_distance_list = []
sheet_distances = {}
sheet_distance_list = []
globular_distances = {}
globular_distance_list = []

# record distances between key locations and all alpha carbons
for residue in chain:
	hydride_distance = hydride - residue['CA']
	adenosine_distance = adenosine - residue['CA']
	m20_distance = m20 - residue['CA']
	sheet_distance = sheet - residue['CA']
	globular_distance = globular - residue['CA']
	name = residue.get_id()
	hydride_distances[hydride_distance] = str(name[1])
	hydride_distance_list.append(hydride_distance)
	adenosine_distances[adenosine_distance] = str(name[1])
	adenosine_distance_list.append(adenosine_distance)
	m20_distances[m20_distance] = str(name[1])
	m20_distance_list.append(m20_distance)
	sheet_distances[sheet_distance] = str(name[1])
	sheet_distance_list.append(sheet_distance)
	globular_distances[globular_distance] = str(name[1])
	globular_distance_list.append(globular_distance)


distance_sets = {}
distance_sets["Hydride Transfer"] = [hydride_distances, hydride_distance_list]
distance_sets["Adenosine Ring"] = [adenosine_distances, adenosine_distance_list]
distance_sets["M20 Loop"] = [m20_distances, m20_distance_list]
distance_sets["Sheet domain"] = [sheet_distances, sheet_distance_list]
distance_sets["Globular domain"] = [globular_distances, globular_distance_list]

name_order = []
origin_groups = {}
groupfile = open(groupfilename, 'r')
for line in groupfile:
	name, color, marker, residues = line.split()
	residues = residues.split(',')
	origin_groups[name] = {'residues':residues, 'size':len(residues), 'color':color, 'counts':[], 'distances':[], 'xs':[]}
	name_order.append(name)
groupfile.close()


for setname in distance_sets.keys():
	distances, distance_list = distance_sets[setname]
	groups = copy.deepcopy(origin_groups)
	all = groups['All']['residues']
	steps = []
	step = -0.0
	step_size = 0.25  #0.25 for raw  #2.5 for difference
	while len(all) > 0:
		for distance in distance_list:
			if (distance <= step) and (distance > step - step_size):
				resnum = distances[distance]
				counter = 0
				for name in groups.keys():
					counter += 1
					if resnum in groups[name]['residues']:
						groups[name]['residues'].remove(resnum)
						groups[name]['distances'].append(distance)
						groups[name]['xs'].append(counter -0.1 + 0.2*random.random())
		steps.append(step)
		for name in groups.keys():
			groups[name]['counts'].append(1.0 - len(groups[name]['residues'])/groups[name]['size'])		
		step = step + step_size
	shortname = setname.split()[0]

	#Plot 1 -- 
	sns.set()
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_ylabel("Percent of Residues in Shell", fontsize = 16)
	ax.set_xlabel("Distance from "+ setname, fontsize = 16)
	plt.tick_params(axis='both', which='major', labelsize=14)
	for name in name_order:
		plt.plot(steps, groups[name]['counts'], alpha = 1.0, label=name, linestyle = "-", c=groups[name]['color'], linewidth=4)
	leg = ax.legend(loc = 'lower center', bbox_to_anchor = (0.5, -0.75), fontsize = 16)
	figname =  "CategoryDistanceFrom" + shortname + "_trace.pdf"
	fig.savefig(figname, bbox_inches='tight')
	plt.show()

	#Plot 2 -- 
	width = 0.5*len(name_order) + 1
	fig = plt.figure(figsize=(width,6))
	ax = fig.add_subplot(111)
	ax.set_ylabel("Distance from "+ setname, fontsize = 16)
	plt.tick_params(axis='both', which='major', labelsize=14)
	colors = []
	distances = []
	xs = []
	for name in name_order:
		distances.append(groups[name]['distances'])
		colors.append(groups[name]['color'])
		xs.append(groups[name]['xs'])
		plt.scatter(groups[name]['xs'], groups[name]['distances'], c=groups[name]['color'])	
	plt.boxplot(distances,labels=name_order)
	plt.xticks(rotation=60)
	figname =  "CategoryDistanceFrom" + shortname + "_boxplot.pdf"
	fig.savefig(figname, bbox_inches='tight')
	plt.tight_layout()
	plt.show()
