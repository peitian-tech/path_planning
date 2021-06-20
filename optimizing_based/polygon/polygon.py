# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import copy

gridEnvList = np.array([[1,1], [1,0], [1,-1], [0,1],[0,0], [0,-1], [-1,1], [-1,0], [-1,-1]])

def float_to_10(a):
	if a>0.00001:
		return 1
	elif a<-0.00001:
		return -1
	else:
		return 0

class Edge():
	def __init__(self, a, b):
		self.a = a
		self.b = b
		ab = b-a
		ao = -a
		proj = ab.dot(ao)/np.linalg.norm(ab)**2
		per_p = a + proj*ab # the perperdicular point from the origin to the line ab
		per_p_norm = np.linalg.norm(per_p)
		if per_p_norm > 0.00001:
			self.distance = per_p_norm
			self.normal = per_p/per_p_norm
		else:
			ab = ab/np.linalg.norm(ab)
			self.normal = np.array([-ab[1], ab[0]])
			self.distance = 0
		self.index = -1


class Polygon():
	def __init__(self, printinfo = True):
		self.vertex = []
		self.vert_concise = False # three neighbouring vertices are not colinear
		self.anticlockwise = False
		self.convex = False
		self.fromA = [] # used in the GJK-EPA
		self.fromB = [] # used in the GJK-EPA
		self.edges = [] # the edges in this polygon
		self.printinfo = printinfo

	def set_vertexes(self, vers, check_flag = True):
		self.vertex = vers
		if check_flag:
			self.check_order()

	def update_edges(self):
		self.edges = []
		for i in range(len(self.vertex)):
			edgei = Edge(self.get_vertex(i), self.get_vertex(i+1))
			edgei.index = i
			self.edges.append(edgei)

	def update_edges_index(self):
		for i in range(len(self.edges)):
			self.edges[i].index = i
			
	def find_closest_edge(self):
		edge_id = 0
		dist = self.edges[0].distance
		for i in range(1, len(self.edges)):
			disttemp = self.edges[i].distance
			if disttemp < dist:
				dist = disttemp
				edge_id = i
		return self.edges[edge_id], edge_id

	def check_order(self):
		# TODO: should be one simply connected domain 
		vertices = np.array(self.vertex)
		xmin_ids = np.where(vertices[:,0] == np.min(vertices[:,0]))[0]
		xmax_ids = np.where(vertices[:,0] == np.max(vertices[:,0]))[0]
		ymin_ids = np.where(vertices[:,1] == np.min(vertices[:,1]))[0]
		ymax_ids = np.where(vertices[:,1] == np.max(vertices[:,1]))[0]
		xyminmax_ids = np.concatenate((xmin_ids, xmax_ids, ymin_ids, ymax_ids))
		cross_values = np.zeros(len(xyminmax_ids))
		for i in range(len(xyminmax_ids)):
			idd = xyminmax_ids[i]
			vec1 = self.get_vertex(idd)-self.get_vertex(idd-1)
			vec2 = self.get_vertex(idd+1)-self.get_vertex(idd)
			if np.cross(vec1, vec2)>=0:
				cross_values[i] = 1
			else:
				cross_values[i] = -1
		if np.sum(cross_values) == len(xyminmax_ids):
			if self.printinfo:
				print("All marginal vertexes are convex vertexes.")
		elif np.sum(cross_values) == -1*len(xyminmax_ids):
			if self.printinfo:
				print("All marginal vertexes have cross value -1. Inverse the order.")
			self.inverse_order()
		else:
			if self.printinfo:
				print("\033[0;31m WARNING: THE VERTEXES ARE IN WRONG ORDER!! USE THE MINIMAL CONVEX HULL INSTEAD! \033[0m")
			self.find_min_convex()
			# os._exit(1)

	def find_min_convex(self):
		'''
		find the mininal convex hull using the Graham scan algorithm
		'''
		vertexes = np.array(self.vertex)
		min_y_id = np.where(vertexes[:,1]==np.min(vertexes[:,1]))[0][0]  # the index of vertex with mininal y
		cos_len_value = np.zeros([len(self.vertex), 2]) # costheta, distance
		cos_len_value[min_y_id, 0] = 2 # any value that is bigger than 1
		cos_len_value[min_y_id, 1] = 0
		for i in range(len(self.vertex)):
			ymin_to_p_vec = vertexes[i] - vertexes[min_y_id]
			ymin_to_p_vec_len = np.linalg.norm(ymin_to_p_vec)
			if ymin_to_p_vec_len > 0.00001:
				ymin_to_p_vec = ymin_to_p_vec/ymin_to_p_vec_len
				cos_len_value[i, 0] = ymin_to_p_vec[0] #ymin_to_p_vec.dot([1,0])
			cos_len_value[i, 1] = ymin_to_p_vec_len
		# sort cos_len_value according to the cos_len_value[:, 0] then cos_len_value[i, 1] in descending order
		sort_order = np.lexsort((-1*cos_len_value[:,1],-1*cos_len_value[:,0]))
		new_vertex = [self.vertex[sort_order[0]]]
		for i in range(1, len(self.vertex)):
			# in when cos_len_value[j,0] = cos_len_value[k,0], they are colinear with the first point, select the one far from the first point
			if cos_len_value[sort_order[i],0] != cos_len_value[sort_order[i-1],0]:
				to_add = self.vertex[sort_order[i]]
				if len(new_vertex)>=2:
					ver1 = np.array(new_vertex[-1]) - np.array(new_vertex[-2]) 
					ver2 = np.array(to_add) - np.array(new_vertex[-1])
					if np.cross(ver1, ver2) > 0:
						new_vertex.append(to_add)
					else:
						new_vertex[-1] = to_add
				else:
					new_vertex.append(to_add)
		if self.printinfo:
			print("There are %d/%d vertexes in the minimal convex hull."%(len(new_vertex), len(self.vertex)))
		self.vertex = new_vertex
		self.anticlockwise = True
		self.convex = True
		self.vert_concise = True

	def add_vertex(self, ver):
		self.vertex.append(ver)

	def remove_vertex(self, ver):
		self.vertex.remove(ver)

	def remove_vertex_byi(self, index):
		i = index%len(self.vertex)
		self.vertex.pop(i)

	def remove_ab_byi(self, index):
		i = index%len(self.vertex)
		self.fromA.pop(i)
		self.fromB.pop(i)

	def get_vertex(self, index):
		i = index%len(self.vertex)
		return np.array(self.vertex[i])

	def get_fromA(self, index):
		i = index%len(self.vertex)
		return np.array(self.fromA[i])

	def get_fromB(self, index):
		i = index%len(self.vertex)
		return np.array(self.fromB[i])

	def inverse_order(self):
		self.vertex.reverse()

	def get_far_proj_vertex(self, direction):
		projmax = np.dot(self.get_vertex(0), direction)
		projmax_p = self.get_vertex(0)
		for i in range(1, len(self.vertex)):
			projtemp = np.dot(self.get_vertex(i), direction)
			if projtemp > projmax:
				projmax = projtemp
				projmax_p = self.get_vertex(i)
		return projmax_p

	def contains(self, p):
		p = np.array(p)
		# check if p is within polygon
		# Currently only support polygon with 3 vertices
		if len(self.vertex) != 3:
			print("\033[0;31m ERROR: CURRENTLY FUNCTION CONTAINS() ONLY SUPPORTS POLYGON WITH 3 VERTICES!! \033[0m")
			os._exit(1)
		p0 = self.get_vertex(0)
		p1 = self.get_vertex(1)
		p2 = self.get_vertex(2)
		# clockwise of anticlockwise
		p01 = p1-p0
		p12 = p2-p1
		plg_dir = float_to_10(np.cross(p01, p12)) # 1, -1, 0
		for i in range(3):
			v1 = self.get_vertex(i+1)-self.get_vertex(i)
			v2 = p-self.get_vertex(i)
			point_dir = float_to_10(np.cross(v1, v2))
			if point_dir==0:
				return True
			elif plg_dir != point_dir:
				return False
		return True

	def plot(self):
		plg_vers = np.array(self.vertex+[self.vertex[0]])
		plt.plot(plg_vers[0,0], plg_vers[0,1],'ro')
		plt.plot(plg_vers[:,0], plg_vers[:,1],'.-')
		plt.plot(plg_vers[-2,0], plg_vers[-2,1],'k*')
		plt.show()

	def plot_setlim(self, xmin, xmax, ymin, ymax):
		plg_vers = np.array(self.vertex+[self.vertex[0]])
		plt.plot(plg_vers[0,0], plg_vers[0,1],'ro')
		plt.plot(plg_vers[:,0], plg_vers[:,1],'.-')
		plt.plot(plg_vers[-2,0], plg_vers[-2,1],'k*')
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.show()

	def plot_in_ax(self, ax):
		plg_vers = np.array(self.vertex+[self.vertex[0]])
		ax.plot(plg_vers[0,0], plg_vers[0,1],'ro')
		ax.plot(plg_vers[:,0], plg_vers[:,1],'.-')
		ax.plot(plg_vers[-2,0], plg_vers[-2,1],'k*')

def plg_remove_colin(plg):
	'''
	remove the vertex that is colinear with other two neighbouring vertex

	@param plg: an instance of class Polygon possibly with at least three neighbouring vertices that are colinear
	@return plg: an instance of class Polygon, represents the same polygon as the input instance but with no colinear vertex
	@Example: plg1 = plg_remove_colin(plg)
	'''
	colin = True
	while colin:
		colin = False
		for i in range(len(plg.vertex)):
			vec1 = plg.get_vertex(i+1) - plg.get_vertex(i)
			vec2 = plg.get_vertex(i+2) - plg.get_vertex(i)
			if np.cross(vec1, vec2)==0:
				colin = True
				plg.remove_vertex_byi(i+1)
				# plg.remove_vertex(list(plg.get_vertex(i+1)))
				break
	plg.vert_concise = True
	return plg

def plg_anticlockwise(plg):
	'''
	reorder the vertices in plg in anticlockwise direction

	@param plg: an instance of class Polygon 
	@return plg: an instance of class Polygon, represents the same polygon as the input instance, the vertices are in anticlockwise order
	@Example: plg1 = plg_remove_colin(plg)
	'''
	vertices = np.array(plg.vertex)
	ylist = list(vertices[:,1])
	min_y_id = ylist.index(min(ylist)) # the index of vertex with mininal y
	vec1 = plg.get_vertex(min_y_id)-plg.get_vertex(min_y_id-1)
	vec2 = plg.get_vertex(min_y_id+1)-plg.get_vertex(min_y_id)

	plg_new_vers = []
	if np.cross(vec1, vec2)>0:
		# plg is already anticlockwise
		for i in range(len(plg.vertex)):
			plg_new_vers.append(list(plg.get_vertex(min_y_id+i)))
	elif np.cross(vec1, vec2)<0:
		for i in range(len(plg.vertex)):
			plg_new_vers.append(list(plg.get_vertex(min_y_id-i)))
	else:
		print("\033[0;31m ERROR: REMOVE THE COLINEAR VERTEX FIRST!! \033[0m")
		os._exit(1)
	plg.set_vertexes(plg_new_vers)
	plg.anticlockwise = True
	return plg

def move_along_dir_by_dis(plg, direction, distance):
	plg1 = copy.deepcopy(plg)
	dx = direction[0]*distance
	dy = direction[1]*distance
	for i in range(len(plg1.vertex)):
		plg1.vertex[i][0] = plg1.vertex[i][0] + dx
		plg1.vertex[i][1] = plg1.vertex[i][1] + dy
	return plg1

def plg_to_convexplg(plg, convex_plg_list, printinfo = True):
	'''
	check if a polygon is convex or concave, if it is concave, decompose it into several convex parts

	@param plg: an instance of class Polygon, with no colinear vertices, anticlockwise order
	@param convex_plg_list: an empty list, store the new-generated convex parts
	@return: no return ,the results will be stored in convex_plg_list
	@Example: convex_plg_list = []; plg_to_convexplg(plg, convex_plg_list)
	'''
	plg_cross_val = np.zeros(len(plg.vertex))
	for i in range(len(plg.vertex)):
		vec1 = plg.get_vertex(i) - plg.get_vertex(i-1)
		vec2 = plg.get_vertex(i+1) - plg.get_vertex(i)
		if np.cross(vec1, vec2)>0:
			plg_cross_val[i] = 1
		elif np.cross(vec1, vec2)<0:
			plg_cross_val[i] = -1
		else:
			print("\033[0;31m ERROR: REMOVE THE COLINEAR VERTEX FIRST!! \033[0m")
			os._exit(1)

	cross_sum = np.sum(plg_cross_val)
	ver_num = len(plg.vertex)
	if cross_sum == ver_num:
		if printinfo:
			print("The input polygon is convex. sum=%d, num=%d"%(np.sum(plg_cross_val), len(plg.vertex)))
		plg.convex = True
		convex_plg_list.append(plg)
		return

	if printinfo:
		print("The input polygon is concave with %d convex vertices and %d concave vertices!!"%((ver_num+cross_sum)/2, (ver_num-cross_sum)/2))
	concave_id_list = np.where(plg_cross_val==-1)[0]
	concave_i = 0
	find_cut_line = False
	while (not find_cut_line) and concave_i < len(concave_id_list):
		concave_vert_id = concave_id_list[concave_i]
		# extend line: vertex[concave_vert_id-1] -- vertex[concave_vert_id]
		# line param: a, b, c, ax+by+c=0
		p_before = plg.get_vertex(concave_vert_id-1)
		p_concave = plg.get_vertex(concave_vert_id)
		if p_before[0] == p_concave[0]:
			a = 1; b = 0; c =-p_before[0];
		else:
			a = (p_before[1]-p_concave[1])/(p_before[0]-p_concave[0]);
			b = -1;
			c = p_before[1] - a*p_before[0];

		for i in range(concave_vert_id+1, concave_vert_id+1+ver_num):
			p1 = plg.get_vertex(i)
			p2 = plg.get_vertex(i+1)
			if (a*p1[0]+b*p1[1]+c)*(a*p2[0]+b*p2[1]+c) <= 0:
				cut_id = i+1
				break;
		# check if vertex concave_vert_id can see vertex cut_id
		cansee = True
		# vector from this concave point to the cut point
		vec1 = plg.get_vertex(cut_id)-plg.get_vertex(concave_vert_id)
		for i in range(cut_id+1, len(plg.vertex)):
			vec2 = plg.get_vertex(i)-plg.get_vertex(concave_vert_id)
			if np.cross(vec1, vec2)<0:
				cansee = False
				break
		if cansee:
			find_cut_line = True
			break;
		else:
			# the connection between vertex concave_vert_id and vertex cut_id passes through one edge of the polygon
			# discard this connection and use the next concave vertex
			concave_i = concave_i+1 
	if not find_cut_line:
		print("\033[0;31m ERROR: CANNOT FIND A CONNECTION TO CUT THE POLYGON!! THIS ALGORITHM IS SHIT! \033[0m")
		os._exit(1)

	# cut the polygon by connecting vertex(concave_vert_id) and vertex(cut_id)
	plg1 = Polygon(printinfo)
	plg2 = Polygon(printinfo)
	for i in range(concave_vert_id, cut_id+1):
		plg1.add_vertex(list(plg.get_vertex(i)))
	for i in range(cut_id, ver_num+concave_vert_id+1):
		plg2.add_vertex(list(plg.get_vertex(i)))
	# remove the possible colinear and neighbouring vertices
	plg1 = plg_remove_colin(plg1)
	plg2 = plg_remove_colin(plg2)
	# recursively decomposite plg1 and plg2 if they are concave
	plg1 = plg_to_convexplg(plg1, convex_plg_list)
	plg2 = plg_to_convexplg(plg2, convex_plg_list)

def map_to_polygons(map):
	edge_p = []
	d_check = 0.05
	for i in range(6, 95):
		for j in range(6, 95):
			status_arr=[]
			for k in range(9):
				if map.if_on_obstacle(i+d_check*gridEnvList[k][0], j+d_check*gridEnvList[k][1]):
					status_arr.append(1)
				else:
					status_arr.append(0)
				if (0 in status_arr) and (1 in status_arr):
					edge_p.append([i,j]); break;

	for i in range(6, 95):
		if map.if_on_obstacle(5+d_check, i-d_check) or map.if_on_obstacle(5+d_check, i+d_check):
			edge_p.append([5,i])
		if map.if_on_obstacle(95-d_check, i-d_check) or map.if_on_obstacle(95-d_check, i+d_check):
			edge_p.append([95,i])
		if map.if_on_obstacle(i-d_check, 5+d_check) or map.if_on_obstacle(i+d_check, 5+d_check):
			edge_p.append([i,5])
		if map.if_on_obstacle(i-d_check, 95-d_check) or map.if_on_obstacle(i+d_check, 95-d_check):
			edge_p.append([i,95])

	edge_p_copy = edge_p[:]
	# extract polygons from edge points
	polygons = []
	# add the walls as polygons
	plg = Polygon(False)
	plg.add_vertex([0,0]);	plg.add_vertex([100,0]);	plg.add_vertex([100,5]);	plg.add_vertex([0,5])
	polygons.append(plg)

	plg = Polygon(False)
	plg.add_vertex([0,5]);	plg.add_vertex([5,5]);	plg.add_vertex([5,95]);	plg.add_vertex([0,95])
	polygons.append(plg)	

	plg = Polygon(False)
	plg.add_vertex([0,95]);	plg.add_vertex([100,95]);	plg.add_vertex([100,100]);	plg.add_vertex([0,100])
	polygons.append(plg)	

	plg = Polygon(False)
	plg.add_vertex([95,5]);	plg.add_vertex([100,5]);	plg.add_vertex([100,95]);	plg.add_vertex([95,95])
	polygons.append(plg)

	while len(edge_p)>0:
		plg = Polygon()
		ver_p = edge_p[0]
		# print("first vertex: ", ver_p)
		plg.add_vertex( ver_p )
		edge_p.remove( ver_p )
		findNewVex = True
		while findNewVex:
			findNewVex = False
			xp_dir = False
			xn_dir = False
			yp_dir = False
			yn_dir = False

			if [ver_p[0]+1, ver_p[1]] in edge_p:
				xp_dir = True
			if [ver_p[0]-1, ver_p[1]] in edge_p:
				xn_dir = True
			if [ver_p[0], ver_p[1]+1] in edge_p:
				yp_dir = True
			if [ver_p[0], ver_p[1]-1] in edge_p:
				yn_dir = True
			
			i = 1
			while (not findNewVex) and xp_dir:
				if [ver_p[0]+i, ver_p[1]] in edge_p:
					edge_p.remove([ver_p[0]+i, ver_p[1]])
					i=i+1
				else:
					plg.add_vertex([ver_p[0]+i-1, ver_p[1]])
					findNewVex = True
					ver_p = [ver_p[0]+i-1, ver_p[1]]
					xp_dir = False
			i = 1
			while (not findNewVex) and xn_dir:
				if [ver_p[0]-i, ver_p[1]] in edge_p:
					edge_p.remove([ver_p[0]-i, ver_p[1]])
					i=i+1
				else:
					plg.add_vertex([ver_p[0]-i+1, ver_p[1]])
					findNewVex = True
					ver_p = [ver_p[0]-i+1, ver_p[1]]
					xp_dir = False
			i = 1
			while (not findNewVex) and yp_dir:
				if [ver_p[0], ver_p[1]+i] in edge_p:
					edge_p.remove([ver_p[0], ver_p[1]+i])
					i=i+1
				else:
					plg.add_vertex([ver_p[0], ver_p[1]+i-1])
					findNewVex = True
					ver_p = [ver_p[0], ver_p[1]+i-1]
					yp_dir = False
			i = 1
			while (not findNewVex) and yn_dir:
				if [ver_p[0], ver_p[1]-i] in edge_p:
					edge_p.remove([ver_p[0], ver_p[1]-i])
					i=i+1
				else:
					plg.add_vertex([ver_p[0], ver_p[1]-i+1])
					findNewVex = True
					ver_p = [ver_p[0], ver_p[1]-i+1]
					yn_dir = False
		polygons.append(plg)

	print("There are %d polygons in this map."%len(polygons))

	vertexes = []
	for i in range(len(polygons)):
		for j in range(len(polygons[i].vertex)):
			vertexes.append(polygons[i].vertex[j])

	print("Decomposite the concave polygons into convex subpolygons")
	polygons_convex = []
	for i in range(len(polygons)):
		polygons[i] = plg_remove_colin(polygons[i])
		polygons[i] = plg_anticlockwise(polygons[i])
		plg_list=[]
		plg_to_convexplg(polygons[i], plg_list)
		for j in range(len(plg_list)):
			polygons_convex.append(plg_list[j])
	print("There are %d convex polygons in this map."%len(polygons_convex))

	vertexes = []
	for i in range(len(polygons_convex)):
		for j in range(len(polygons_convex[i].vertex)):
			vertexes.append(polygons_convex[i].vertex[j])

	for i in range(len(polygons_convex)):
		print("%dth convex polygon has %d vertexes"%(i,len(polygons_convex[i].vertex)))

	# return edge_p_copy, vertexes, polygons_convex
	return polygons_convex