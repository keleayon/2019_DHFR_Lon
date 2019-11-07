import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math
import numpy as np
import random
from operator import itemgetter
from itertools import combinations


def parseMap(filename, high_limit):
	labels = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*', '','avg']
	selection_map = pd.read_csv(filename, index_col = 'Position  ')
	selection_map.columns = labels
	selection_map.clip_upper(high_limit)
	muts = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*']
	selection_map = selection_map[muts]
	
	return selection_map


def getVectorFromMap(map, position, include_nan=False):
	""""""
	if isinstance(position, float): #Make this an assert statement instead???
		position = int(position)
	vector = np.array(map.iloc[position, :].values).flatten()
	if include_nan:
		return vector
	elif not include_nan:
		return vector[~np.isnan(vector)]


def calculateVectorVectorDistance(vector1, vector2):
	""""""
	length = np.amin([len(vector1),len(vector2)])
	vector1 = -np.sort(-vector1)
	vector2 = -np.sort(-vector2)
	#implement random value removal numpy delete function
	vector1 = vector1[:length]
	vector2 = vector2[:length]
	
	return np.mean(np.abs(vector1 -vector2))


def evaluateDistances(distances, method="Top3"):
	""""""
	if method=="Top3" or method=="top3":
		distances = np.sort(distances)
		return np.mean(distances[:3])	
	elif method=="Top2" or method=="top2":
		distances = np.sort(distances)
		return np.mean(distances[:2])	
	elif method =="mean" or method =="Mean":
		return np.mean(distances)	
	elif method =="median" or method =="Median":
		return np.median(distances)		
	elif method =="mixed" or method =="Mixed":
		return np.mean(distances) + np.mean(distances[:2])	
	else:
		print(f"{method} is not an acceptable evaluation method.")
		exit()


def calculateVectorClusterDistance(vector1, cluster_members, map):
	"""The vector is a iterable object of fitness values for one position.\ 
	The cluster is an iterable object of interger values for positions in the fitness"""

	distances = []
	for position in cluster_members:
		vector = getVectorFromMap(map, position)
		distances.append(calculateVectorVectorDistance(vector1, vector))
	
	return evaluateDistances(distances)


def calculateClusterClusterDistance(cluster1_members, cluster2_members, map):
	
	distances = []
	for position in cluster1_members:
		vector = getVectorFromMap(map, position)
		distances.append(calculateVectorClusterDistance(vector, cluster2_members, map))
		
	for position in cluster2_members:
		vector = getVectorFromMap(map, position)
		distances.append(calculateVectorClusterDistance(vector, cluster1_members, map))

	return evaluateDistances(distances, method="mean")


def assignCluster(map, position, clusters):
	""""""
	vector = getVectorFromMap(map, position)
	rankings = []
	for cluster in clusters.keys():
		rankings.append((cluster, calculateVectorClusterDistance(vector, clusters[cluster], map)))
	rankings.sort(key=itemgetter(1))
	return rankings[0][0] #Identity of the closest cluster
	

def clustering(map, clusters, total, rounds):
	""""""
	#print(clusters)
	for round in range(rounds):
		new_clusters = {}
		for position in range(total):
			closest = assignCluster(map, position, clusters)
			try:
				new_clusters[closest].append(position)
			except KeyError:
				new_clusters[closest] = [position]
		clusters = new_clusters

	return clusters


def generateRandomClusters(map, total, num_clusters):
	""""""
	positions = random.sample(range(total), num_clusters)
	clusters = {}
	for position in positions:
		print(position)
		clusters[position] = [position]
	
	return clusters


def initializeRandomStart(map, total, num_clusters):
	""""""
	clusters = generateRandomClusters(map, total, num_clusters)
	return clustering(map, clusters, total, 1)

		
def generateCategories(length):
	""""""
	categories = {}
	categories['Beneficial1'] = np.linspace(high, 0.0, length)
	categories['Beneficial2'] = np.linspace(high/2, 0.0, length)
	categories['Beneficial3'] = np.linspace(high/3, 0.0, length)
	categories['Beneficial4'] = np.linspace(high/2, low/8, length)
	categories['Tolerant1'] = np.linspace(0.0, 0.0, length)
	categories['Mixed1'] = np.linspace(high, low, length)
	categories['Mixed2'] = np.linspace(high/2, low/2, length)
	categories['Mixed3'] = np.linspace(high/2, low/3, length)
	categories['Mixed4'] = np.linspace(high/3, low/3, length)
	#categories['Mixed5'] = np.linspace(high/4, low/4, length)
	#Add more mixed categories?
	categories['Deleterious1'] = np.linspace(0, low/2, length)
	categories['Deleterious2'] = np.linspace(0, low, length)
	try:
		categories['Deleterious3'] = np.linspace(low, low, length)
		categories['Deleterious3'][0] = 0.0 #Wild-type
	except IndexError:
		pass
	
	try:
		categories['Deleterious4'] = np.linspace(low, dead, length)
		categories['Deleterious4'][0] = 0.0 #Wild-type
	except IndexError:
		pass
		
	try:
		Intolerant = np.linspace(dead, dead, length)
		Intolerant[0] = 0.0 # Wild-type
		categories['Intolerant1'] = Intolerant
	except IndexError:
		pass
	
	try:
		Intolerant[1] = 0.0 # Wild-type + one acceptable mutation
	except IndexError:
		pass
		
	categories['Intolerant2'] = Intolerant 
	try:	
		Intolerant[2] = 0.0 # Wild-type + two acceptable mutations
	except IndexError:
		pass
		
		#categories['Intolerant3'] = Intolerant 
	
	return categories


def initializeOnCategories(map, total):
	""""""
	clusters = {}

	for position in range(total):
		rankings = []
		vector = getVectorFromMap(map, position)
		categories = generateCategories(len(vector))
		for category in categories.keys():
			rankings.append((category, calculateVectorVectorDistance(vector, categories[category])))
		rankings.sort(key=itemgetter(1))
		closest = rankings[0][0]
		try:
			clusters[closest].append(position)
		except KeyError:
			clusters[closest] = [position]		
		
	return clustering(map, clusters, total, 1)


def condenseRelatedClusters(clusters):
	""""""
	condensed_clusters = {}
	for cluster in clusters.keys():
		if isinstance(cluster, int):
			new_cluster = cluster #for random clusters with int names
		elif isinstance(cluster, str):
			new_cluster = cluster[:-1] # for categories with string names
	
		try:
			condensed_clusters[new_cluster].extend(clusters[cluster])
		except KeyError:
			condensed_clusters[new_cluster] = clusters[cluster]
		
	return condensed_clusters


def condenseClustersToNum(clusters, num_clusters, map):
	while len(clusters.keys()) > num_clusters:
		new_clusters = {}
		pairs = list(combinations(clusters.keys(), 2))
		rankings = []
		for pair in pairs:
			rankings.append((pair, calculateClusterClusterDistance(clusters[pair[0]], clusters[pair[1]], map)))
		rankings.sort(key=itemgetter(1))
		closest = rankings[0][0]
		labels = list(clusters.keys())
		labels.remove(closest[0])
		labels.remove(closest[1])
		new_label = f"{closest[0]}-{closest[1]}"
		data = []
		data.extend(clusters[closest[0]])
		data.extend(clusters[closest[1]])
		new_clusters[new_label] = data
		
		for cluster in labels:
			new_clusters[cluster] = clusters[cluster]
		
		clusters = new_clusters
		
	return clusters

#make all of this argparse
mapfile = sys.argv[1]
num_residues = int(sys.argv[2])
high_limit = float(sys.argv[3])
high = float(sys.argv[4])
low = float(sys.argv[5])
dead = float(sys.argv[6])
cluster_type = str(sys.argv[7])
if cluster_type == "Random" or cluster_type == "random":
	num_clusters = int(sys.argv[8])
	num_final_clusters = int(sys.argv[9])


selection_map = parseMap(mapfile, high_limit)

	
if cluster_type == "Categories" or cluster_type == "categories":
	seed_clusters = initializeOnCategories(selection_map, num_residues)
elif cluster_type == "Random" or cluster_type == "random":
	seed_clusters = initializeRandomStart(selection_map, num_residues, num_clusters)

first_clusters = clustering(selection_map, seed_clusters, num_residues, 10)

if cluster_type == "Categories" or cluster_type == "categories":
	condensed_clusters = condenseRelatedClusters(first_clusters)
elif cluster_type == "Random" or cluster_type == "random":
	condensed_clusters = condenseClustersToNum(first_clusters, num_final_clusters, selection_map)


second_clusters = clustering(selection_map, condensed_clusters, num_residues, 10)

for cluster in second_clusters:
	print(cluster)
	outcluster = []
	for pos in second_clusters[cluster]:
		outcluster.append(str(pos+1))
	print(','.join(outcluster))
	print('')
	
		
