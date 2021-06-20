# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import time
from polygon import *
from gjk_epa import *
from cvxopt import solvers, matrix

sys.path.append("../") 
import my_map as mymap

def test_polygon_self():
	'''
	test item: the functions of polygon class
	'''
	plg = Polygon()
	# vers=[[0,0],[1,1],[2,0],[3,3],[0,4],[-1,3.5], [0,2]]
	# vers = [[0,0],[1,1],[2,0],[3,4],[2,1.5],[0,3]]
	vers = [[1,2], [2,2], [2.5,3], [1.5, 4], [0.5,3], [4,1], [5,1], [5.5,2], [4.5,3], [3.5,2]]
	plg.set_vertexes(vers)
	plg.find_min_convex()
	plg1 = move_along_dir_by_dis(plg, [1,0], 5)
	fig, ax = plt.subplots()
	plg.plot_in_ax(ax) # plot this plg
	plg1.plot_in_ax(ax) # plot this plg
	plt.show()

def test_decomposition():
	'''
	test item: the decomposition of a concave polygon
	'''
	plg = Polygon()
	# vers=[[0,0],[1,1],[2,0],[3,3],[0,4],[-1,3.5], [0,2]]
	vers=[[0,0],[1,1],[2,0],[3,4],[2,1.5],[0,3]]
	plg.set_vertexes(vers)
	plg.plot() # plot this plg

	# decomposition the polygon
	convex_plg_list=[]
	plg_to_convexplg(plg, convex_plg_list)

	# plot every convex polygon
	# for i in range(len(convex_plg_list)):
	# 	convex_plg_list[i].plot()

	# plot the polygon the every convex subpolygon
	plt.ion()
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plg.plot_in_ax(ax)
	plt.title("Decompositon of a Concave Polygon")

	subplg, = ax.plot([], [], 'r--')
	fig.canvas.draw()
	fig.canvas.flush_events()
	for i in range(len(convex_plg_list)):
		showtimes = 1
		if i== len(convex_plg_list)-1:
			showtimes = 10000000000
		j = 0
		while j<showtimes:
			time.sleep(1)
			vertices = convex_plg_list[i].vertex
			vertices = np.array(vertices+[vertices[0]])
			subplg.set_xdata(vertices[:,0])
			subplg.set_ydata(vertices[:,1])
			fig.canvas.draw()
			fig.canvas.flush_events()
			j=j+1

	plt.ioff()

def test_polygon_distance():
	'''
	test item: calculate the distance between two polygons
	'''
	plg1 = Polygon(False); plg2 = Polygon(False)
	plg1.convex = True; plg2.convex = True

	# pentagons
	# plg1.set_vertexes([[1,2], [2,2], [2.5,3], [1.5, 4], [0.5,3]])
	# plg2.set_vertexes([[4,1], [5,1], [5.5,2], [4.5,3], [3.5,2]])
	# plg2.set_vertexes([[2,1], [3,1], [3.5,2], [2.5,3], [1.5,2]])


	# plg1.set_vertexes([[0,5], [5,5], [5,95],[0,95]])
	# plg2.set_vertexes([[7.25, 92.25], [7.75, 92.25], [7.75, 92.75], [7.25, 92.75]])

	# plg2 in plg1
	# plg1.set_vertexes([[5, 90], [10,90], [10,95], [5,95]])
	# plg2.set_vertexes([[7.25, 92.25], [7.75, 92.25], [7.75, 92.75], [7.25, 92.75]])
	
	# wrong epa distance, dis=0 epa dis=-10
	# plg1.set_vertexes([[10,15], [30,15], [30,20], [10,20]])
	# dx = 2
	# plg2.set_vertexes([[30-dx,10],[34-dx,10],[34-dx,25],[30-dx,25]])
	
	# plg1.set_vertexes([[2,2],[2,5],[4,5],[4,2]])
	# plg2.set_vertexes([[3,3],[5,4],[5,2]])
	
	# px = 48.5844; py = 50.1788; a = 2.9; ha=a/2.0;
	px = 56.27; py = 49.5514; a = 2.9; ha=a/2.0;
	plg1.set_vertexes([[px-ha,py-ha], [px+ha, py-ha], [px+ha, py+ha], [px-ha, py+ha]])
	plg2.set_vertexes([[50,30], [55,25], [55,70], [50,70]])

	gjkepa = GJK_EPA(plg1, plg2)

	res = gjkepa.cal_distance()
	print(res)

	fig, ax = plt.subplots()
	plg1.plot_in_ax(ax)
	plg2.plot_in_ax(ax)
	plt.title("distance = %f"%res[1])
	plt.axis("Equal")
	plt.show()

def test_gjkepa_on_map():
	map = mymap.Map([[]])
	map.load_map_file('../../common/map1.txt')
	(start_node , target_node) = map.generate_start_target()

	ob_polygons = map_to_polygons(map)

	fig, ax = plt.subplots()
	for ob_plg in ob_polygons:
		ob_plg.plot_in_ax(ax)
	plt.axis("Equal")
	plt.show()



def test_cvxopt():
	N = 9
	path = np.zeros([N, 2])
	path[0,:] = np.array([7.5, 92.5])
	path[-1,:] = np.array([87.5, 12.5])
	stvec = path[-1,:] - path[0,:]
	len_st = np.linalg.norm(stvec)
	stvec = stvec/len_st
	cos_theta = np.dot(np.array([1,0]), stvec)
	sin_theta = np.cross(np.array([1,0]), stvec)
	tmatrix = np.array([[cos_theta, -sin_theta],[sin_theta, cos_theta]])
	pvec = path[0,:]
	period_num = (int)(N/5)
	omega = 2*np.pi*period_num/len_st
	for i in range(N):
		px = i*len_st/(float)(N-1); py = 2*np.sin(omega*px)
		path[i,:]=np.dot(tmatrix, np.array([px,py])) + pvec

	pathlin = np.zeros([N, 2])
	dx = (path[-1,0] - path[0,0])/(N-1.0)
	dy = (path[-1,1] - path[0,1])/(N-1.0)
	for i in range(N):
		pathlin[i, 0] = path[0,0] + i*dx
		pathlin[i, 1] = path[0,1] + i*dy

	# Q, p, G, h
	Q_opt = np.diag(np.ones(2*(N-2))*4)
	for i in range(np.size(Q_opt, 0)-2):
		Q_opt[i,i+2] = -2
		Q_opt[i+2, i] = -2

	print(Q_opt)

	p_opt = np.zeros(2*(N-2))
	p_opt[0] = -2*path[0,0]
	p_opt[1] = -2*path[0,1]
	p_opt[-2] = -2*path[N-1,0]
	p_opt[-1] = -2*path[N-1,1]	

	G_opt1 = np.diag(np.ones(2*(N-2)))
	G_opt2 = np.diag(-1*np.ones(2*(N-2)))
	G_opt = np.vstack((G_opt1, G_opt2))

	h_opt1 = path[1:-1].reshape(-1,1) + 5
	h_opt2 = path[1:-1].reshape(-1,1)*(-1)+5
	h_opt = np.vstack((h_opt1, h_opt2))

	sol=solvers.qp(matrix(Q_opt), matrix(p_opt), matrix(G_opt), matrix(h_opt))
	pathnew = path.copy()
	pathnew[1:-1,:] = np.array(sol['x']).reshape(-1,2)

	pathsinnorm = 0
	pathlinnorm = 0
	pathnewnorm = 0
	for i in range(N-1):
		pathsinnorm = pathsinnorm + np.linalg.norm(path[i+1,:]-path[i,:])**2
		pathlinnorm = pathlinnorm + np.linalg.norm(pathlin[i+1,:]-pathlin[i,:])**2
		pathnewnorm = pathnewnorm + np.linalg.norm(pathnew[i+1,:]-pathnew[i,:])**2

	pathmid = path[1:-1,:].reshape(-1,1)[:,0]
	pathlinmid = pathlin[1:-1,:].reshape(-1,1)[:,0]
	pathnewmid = pathnew[1:-1,:].reshape(-1,1)[:,0]
	path_sin_norm_cal = pathmid.dot(Q_opt).dot(pathmid)/2.0 + p_opt.dot(pathmid) + np.linalg.norm(path[0,:])**2 + np.linalg.norm(path[-1,:])**2
	path_lin_norm_cal= pathlinmid.dot(Q_opt).dot(pathlinmid)/2.0 + p_opt.dot(pathlinmid) + np.linalg.norm(path[0,:])**2 + np.linalg.norm(path[-1,:])**2
	path_new_norm_cal= pathnewmid.dot(Q_opt).dot(pathnewmid)/2.0 + p_opt.dot(pathnewmid) + np.linalg.norm(path[0,:])**2 + np.linalg.norm(path[-1,:])**2

	print("path_sin_norm=", pathsinnorm, "path_sin_norm_cal=", path_sin_norm_cal)
	print("pathlinnorm=", pathlinnorm, "path_lin_norm_cal=", path_lin_norm_cal)
	print("pathnewnorm=", pathnewnorm, "path_new_norm_cal=", path_new_norm_cal)

	plt.plot(path[:,0], path[:,1],'.-', label="original")
	plt.plot(pathlin[:,0], pathlin[:,1],'.-', label="pathlin")
	plt.plot(pathnew[:,0], pathnew[:,1],'.-', label="new")
	plt.legend(loc=2)
	plt.show()

def test_cvxopt1():
	curpath = np.loadtxt("curpath.txt")
	temppath = np.loadtxt("temppath.txt")
	Q_opt = np.loadtxt("Q_opt.txt")
	p_opt = np.loadtxt("p_opt.txt")
	G_opt = np.loadtxt("G_opt.txt")
	h_opt = np.loadtxt("h_opt.txt")
	sol=solvers.qp(matrix(Q_opt), matrix(p_opt), matrix(G_opt), matrix(h_opt))

	curmid = curpath[1:-1,:].reshape(-1,1)[:,0]
	solmid = temppath[1:-1,:].reshape(-1,1)[:,0]
	aaa = np.array(sol['x']).reshape(-1,1)[:,0]
	curmmm = curmid.copy()
	# curmmm[67] = solmid[67]
	print(curmid)
	plt.plot(curmid.reshape(-1,2)[:,0], curmid.reshape(-1,2)[:,1], label="curpath")
	plt.plot(aaa.reshape(-1,2)[:,0], aaa.reshape(-1,2)[:,1], label="solpath")
	plt.plot(solmid.reshape(-1,2)[:,0], solmid.reshape(-1,2)[:,1], label="solpath1")
	plt.legend()
	plt.show()
	#print(np.array(sol['x']).reshape(-1,2))
	print(curmid.dot(Q_opt).dot(curmid)/2.0 + p_opt.dot(curmid))
	print(solmid.dot(Q_opt).dot(solmid)/2.0 + p_opt.dot(solmid))
	print(aaa.dot(Q_opt).dot(aaa)/2.0 + p_opt.dot(aaa))
	print(curmmm.dot(Q_opt).dot(curmmm)/2.0 + p_opt.dot(curmmm))
	print(curmid.dot(Q_opt).dot(curmid)/2.0 + p_opt.dot(curmid) - solmid.dot(Q_opt).dot(solmid)/2.0 - p_opt.dot(solmid))


if __name__ == '__main__':

	# test_polygon_self()
	# test_decomposition()
	# test_polygon_distance()
	# test_cvxopt()
	# test_cvxopt1()
	test_gjkepa_on_map()
