# -*- coding: utf-8 -*-
'''
This class to calculate the distance between two convex polygons using the GJK and EPA algorithm.
The code is largely rewritten from C# version in https://github.com/youlanhai/learn-physics/tree/master/Assets/05-gjk-epa
Rewritten is permitted by the author.
Reference:
https://github.com/youlanhai/learn-physics
https://blog.csdn.net/you_lan_hai/article/details/108293780
https://dyn4j.org/2010/04/gjk-gilbert-johnson-keerthi/
https://dyn4j.org/2010/04/gjk-distance-closest-points/
https://dyn4j.org/2010/05/epa-expanding-polytope-algorithm/
'''

import os
import numpy as np
import matplotlib.pyplot as plt
from polygon import *

class GJK_EPA():
	def __init__(self, plg1, plg2):
		if not (plg1.convex and plg2.convex):
			raise ValueError("ERROR: THE INPUT POLYGONS SHOULD BE CONVEX!!")
		self.epsilon = 0.00001
		self.maxIterCount = 50
		self.plg1 = plg1
		self.plg2 = plg2
		self.simplex = Polygon()
		self.isCollision = False

	def printSimplexP(self):
		print("simplex_p:")
		for i in range(len(self.simplex.vertex)):
			print(self.simplex.vertex[i])

	def cal_distance(self):
		direction = self.findFirstDirection();
		self.get_support(direction, True)
		self.get_support(-direction, True)
		direction = -self.getClosestPointToOrigin(self.simplex.get_vertex(0), self.simplex.get_vertex(1))
		for i in range(self.maxIterCount):
			if np.linalg.norm(direction) < self.epsilon:
				self.isCollision = True
				break
			p = self.get_support(direction, False)[0]
			if np.linalg.norm(p-self.simplex.get_vertex(0)) < self.epsilon or np.linalg.norm(p-self.simplex.get_vertex(1)) < self.epsilon:
				self.isCollision = False
				break
			self.get_support(direction, True)
			if self.simplex.contains([0,0]):
				self.isCollision = True
				break
			direction = self.findNextDirection()


		if not self.isCollision:
			return self.ComputeClosestPoint()
		else:
			# EPA
			self.simplex.update_edges()
			for i in range(self.maxIterCount):
				ce, ce_id = self.simplex.find_closest_edge() # currentEPAEdge
				penetrationVector = ce.normal * ce.distance
				point = self.get_support(ce.normal, False)[0]
				dist = point.dot(ce.normal)
				if np.abs(dist-ce.distance) < self.epsilon:
					penetrationVector = dist * ce.normal
					p0 = self.simplex.get_vertex(ce_id)
					p1 = self.simplex.get_vertex(ce_id+1)
					p01 = p1-p0
					if np.linalg.norm(p01) < self.epsilon:
						closestOnA = p0
						closestOnB = p0
					else:
						r2 = -p01.dot(p0)/np.linalg.norm(p01)**2
						if r2<0:
							r2=0
						elif r2>1:
							r2=1
						r1 = 1.0-r2
						closestOnA = self.simplex.get_fromA(ce_id) * r1 + self.simplex.get_fromA(ce_id+1)  * r2
						closestOnB = self.simplex.get_fromB(ce_id) * r1 + self.simplex.get_fromB(ce_id+1)  * r2
					return penetrationVector, -np.linalg.norm(penetrationVector), closestOnA, closestOnB
				
				point, frA, frB = self.get_support(ce.normal, False)
				self.simplex.vertex.insert(ce.index+1, list(point))
				self.simplex.fromA.insert(ce.index+1, list(frA))
				self.simplex.fromB.insert(ce.index+1, list(frB))
				self.simplex.update_edges()

				# new_edge0 = Edge(ce.a, point)
				# self.simplex.edges[ce.index] = new_edge0
				# new_edge1 = Edge(point, ce.b)
				# self.simplex.edges.insert(ce.index+1, new_edge1)
				# self.simplex.update_edges_index()
			fig, ax = plt.subplots()
			self.plg1.plot_in_ax(ax)
			self.plg2.plot_in_ax(ax)
			plt.axis("Equal")
			plt.show()
			raise ValueError("ERROR: NO RESULT WHEN EPA, TRY TO INCREASE ITERATION TIMES!!")

	def findFirstDirection(self):
		direction =  self.plg1.get_vertex(0)-self.plg2.get_vertex(0)
		if np.linalg.norm(direction) < self.epsilon:
			direction = self.plg1.get_vertex(1) - self.plg2.get_vertex(0)
		return direction

	def get_support(self, direction, addsimplex):
		proj1max_p = self.plg1.get_far_proj_vertex(direction)
		proj2max_p = self.plg2.get_far_proj_vertex(-direction)
		if addsimplex:
			self.simplex.vertex.append(list(proj1max_p-proj2max_p))
			self.simplex.fromA.append(list(proj1max_p))
			self.simplex.fromB.append(list(proj2max_p))
		return proj1max_p-proj2max_p, proj1max_p, proj2max_p

	def getClosestPointToOrigin(self, p1, p2):
		p12 = p2-p1
		p10 = -p1
		if np.linalg.norm(p12)<self.epsilon:
			return p1
		proj = p12.dot(p10)/np.linalg.norm(p12)**2
		if proj < 0:
			return p1
		elif proj > 1:
			return p2
		else:
			return p1+p12*proj

	def findNextDirection(self):
		if len(self.simplex.vertex) == 2:
			return -self.getClosestPointToOrigin(self.simplex.get_vertex(0), self.simplex.get_vertex(1))
		elif len(self.simplex.vertex) == 3:
			p1 = self.getClosestPointToOrigin(self.simplex.get_vertex(2), self.simplex.get_vertex(0))
			p2 = self.getClosestPointToOrigin(self.simplex.get_vertex(2), self.simplex.get_vertex(1))
			if np.linalg.norm(p1) < np.linalg.norm(p2):
				self.simplex.remove_vertex_byi(1)
				self.simplex.remove_ab_byi(1)
				return -p1
			else:
				self.simplex.remove_vertex_byi(0)
				self.simplex.remove_ab_byi(0)
				return -p2
		else:
			print("\033[0;31m ERROR: THE SIMPLEX DOESNOT HAVE 3 VERTICES!! \033[0m")
			return np.zeros(2)

	def ComputeClosestPoint(self):
		p0 = self.simplex.get_vertex(0)
		p1 = self.simplex.get_vertex(1)
		p01 = p1-p0
		if np.linalg.norm(p01) < self.epsilon:
			closestOnA = p0
			closestOnB = p0
		else:
			r2 = -p01.dot(p0)/np.linalg.norm(p01)**2
			if r2<0:
				r2=0
			elif r2>1:
				r2=1
			r1 = 1.0-r2
			closestOnA = self.simplex.get_fromA(0) * r1 + self.simplex.get_fromA(1)  * r2
			closestOnB = self.simplex.get_fromB(0) * r1 + self.simplex.get_fromB(1)  * r2
		return closestOnA-closestOnB, np.linalg.norm(closestOnA-closestOnB), closestOnA, closestOnB