# -*- coding: utf-8 -*-

import os
import sys
sys.path.append("polygon/") 
import numpy as np
import matplotlib.pyplot as plt
from polygon import *
from gjk_epa import *
from cvxopt import solvers, matrix
import copy
import my_map as mymap

printplgmsg = False

class TrajOpt():
	def __init__(self, env_map, start_node, target_node, side_length, N, init_traj, dsafe, CCD_flag, trick1_flag, trick2_flag, trick3_flag, mu0, s0, c, tau_plus, tau_minus, k, ftol, xtol, ctol):
		self.env_map = env_map
		self.start_node = start_node
		self.target_node = target_node;
		self.N = N # number of trajectory points
		self.a = side_length # the robot is represented by a square with side length a
		self.init_traj = init_traj
		self.ob_polygons = map_to_polygons(env_map) # the convex representation of obstacles in environment
		self.CCD_flag = CCD_flag # continuous collision detection 

		self.trick1 = trick1_flag #check if find a collision-free path without calculating the cost, if true, exit
		self.trick2 = trick2_flag #if this path leads to collision in other path points, update these points in p_opt
		self.trick3 = trick3_flag #if allow safe distance to be 0. (dsafe>0 when optimizing, and dsafe=0 when final checking) will allow faster convergence

		# construct the init trajectory to be optimazied
		path = np.zeros([N, 2])
		if init_traj=='lin':
			dx = (target_node.x - start_node.x)/(N-1.0)
			dy = (target_node.y - start_node.y)/(N-1.0)
			for i in range(N):
				path[i, 0] = start_node.x + i*dx
				path[i, 1] = start_node.y + i*dy
		elif init_traj=='sin':
			stvec = np.array([target_node.x - start_node.x, target_node.y - start_node.y])
			len_st = np.linalg.norm(stvec)
			stvec = stvec/len_st
			cos_theta = np.dot(np.array([1,0]), stvec)
			sin_theta = np.cross(np.array([1,0]), stvec)
			tmatrix = np.array([[cos_theta, -sin_theta],[sin_theta, cos_theta]])
			pvec = np.array([start_node.x, start_node.y])
			period_num = (int)(N/5)
			omega = 2*np.pi*period_num/len_st
			for i in range(N):
				px = i*len_st/(float)(N-1); py = 2*np.sin(omega*px)
				path[i,:]=np.dot(tmatrix, np.array([px,py])) + pvec
		else:
			print("ERROR: unsupported initial trajectory type!")

		self.path = path
		self.rob_polygons = []
		self.rob_swept_polygons = []

		# parameter for TrajOpt algorithm
		self.mu = mu0 #initial penalty coefficient
		self.mutol = 10000000.0
		self.s = s0 #initial trust region size
		self.s0 = s0 #initial trust region size
		self.c = c #step acceptance parameter
		self.tau_plus = tau_plus #trust region expansion factor
		self.tau_minus = tau_minus #trust region shrinkage factor
		self.k = k #penalty scaling factor
		self.ftol = ftol #convergence thresholds for merit
		self.xtol = xtol #convergence thresholds for x
		self.ctol = ctol #constraint satisfaction threshold
		self.dsafe = dsafe # safe distance to the obstacle
		self.dcheck = dsafe + 1 #check distance in every iteration
		self.manual_epa_dirs = [[[]]*len(self.ob_polygons) for ii in range(self.N-1)] # set the epa direction manually in order to overcome the obstacle
		self.rob_ob_infos = self.update_gjkepa_info(self.path) # the distance array from robs to obstacles

		# Q, p, G, h for the QP solver
		Q_opt = np.diag(np.ones(2*(N-2))*4)
		for i in range(np.size(Q_opt, 0)-2):
			Q_opt[i,i+2] = -2
			Q_opt[i+2,i] = -2
		self.Q_opt = Q_opt
		p_opt = np.zeros(2*(N-2))
		p_opt[0] = -2*path[0,0]
		p_opt[1] = -2*path[0,1]
		p_opt[-2] = -2*path[N-1,0]
		p_opt[-1] = -2*path[N-1,1]
		self.p_opt = p_opt
		self.last_real_p_opt = p_opt
		G_opt1 = np.diag(np.ones(2*(N-2)))
		G_opt2 = np.diag(-1*np.ones(2*(N-2)))
		self.G_opt = np.vstack((G_opt1, G_opt2))
		h_opt1 = self.path[1:-1].reshape(-1,1) + self.s
		h_opt2 = self.path[1:-1].reshape(-1,1)*(-1)+self.s
		self.h_opt = np.vstack((h_opt1, h_opt2))

	def update_gjkepa_info(self, path):
		if self.CCD_flag:
			rob_plgs = self.update_rob_swept_polygons(path)
		else:
			rob_plgs = self.update_rob_polygons(path)
		
		m = len(rob_plgs)
		n = len(self.ob_polygons)
		rob_ob_infos = [[0]*n for i in range(m)]
		for i in range(m):
			for j in range(n):
				if self.CCD_flag and len(self.manual_epa_dirs[i][j])>0:
						# if continuous-time detection and the epa_dir is manually set
						# do not use gjkepa to calculate the epa_dir and distance
						epa_vec = self.manual_epa_dirs[i][j][0]
						p_on_rob = rob_plgs[i].get_far_proj_vertex(epa_vec)
						p_on_ob = self.ob_polygons[j].get_far_proj_vertex(-epa_vec)
						epadis = -(p_on_rob-p_on_ob).dot(epa_vec)/np.linalg.norm(epa_vec)
						# print("check: this should be negative: ", epadis)
						rob_ob_infos[i][j] = (epa_vec, epadis, p_on_rob, p_on_ob)
				else:
					gjkepa = GJK_EPA(rob_plgs[i], self.ob_polygons[j])
					rob_ob_infos[i][j] = gjkepa.cal_distance() # vector from ob to rob and norm of this vec
		return rob_ob_infos

	def update_rob_polygons(self, path):
		rob_polygons = []
		ha = self.a/2.0
		for i in range(np.size(path, 0)):
			plg = Polygon(printplgmsg)
			x = path[i,0]; y = path[i,1];
			vers = [[x-ha,y-ha], [x+ha, y-ha], [x+ha,y+ha], [x-ha,y+ha]]
			plg.set_vertexes(vers)
			plg.convex = True
			rob_polygons.append(plg)
		return rob_polygons

	def update_rob_swept_polygons(self, path):
		rob_swept_polygons = []
		ha = self.a/2.0
		for i in range(np.size(path, 0)-1):
			plg = Polygon(printplgmsg)
			xi = path[i,0]; yi = path[i,1];
			xi1 = path[i+1,0]; yi1 = path[i+1,1];
			vers = [[xi-ha,yi-ha], [xi+ha, yi-ha], [xi+ha,yi+ha], [xi-ha,yi+ha], [xi1-ha,yi1-ha], [xi1+ha, yi1-ha], [xi1+ha,yi1+ha], [xi1-ha,yi1+ha]]
			plg.set_vertexes(vers, False)
			plg.find_min_convex()
			rob_swept_polygons.append(plg)
		return rob_swept_polygons

	def compute_ConstrainCost(self, path, safe0, selfpath):
		'''
		compute the constrain cost (obstacle etc.) of the given path
		@param path: the path, Nx2 array
		@param safe0: if allow safe distance to be 0. (dsafe>0 when optimizing, and dsafe=0 when final checking) will allow faster convergence
		@param selfpath: if the input path is current path. 
						 if true, use current rob_ob_information; 
						 if not, recalculate the rob_ob_information of the input path
		@return obstacle_sum: the constrainst cost
		'''
		obstacle_sum = 0
		if selfpath:
			rob_ob_infos = self.rob_ob_infos
		else:
			rob_ob_infos = self.update_gjkepa_info(path)

		m = len(rob_ob_infos)
		n = len(rob_ob_infos[0])
		for i in range(m):
			for j in range(n):
				if safe0:
					obstacle_sum = obstacle_sum + np.max([0-rob_ob_infos[i][j][1], 0])
				else:
					obstacle_sum = obstacle_sum + np.max([self.dsafe-rob_ob_infos[i][j][1], 0])
		return obstacle_sum

	def compute_PathCost(self, path, selfpath):
		'''
		compute the path cost = smooth_cost + constraint_const
		@param path: the path, Nx2 array
		@param selfpath: if the input path is current path. 
						 if true, use current rob_ob_information; 
						 if not, recalculate the rob_ob_information of the input path
		@return the path cost
		'''
		aa=(path[0:-1,:]-path[1:,:]).reshape(-1,1)[:,0]
		square_sum = aa.dot(aa)

		ccost = self.compute_ConstrainCost(path, False, selfpath)
		return self.mu * ccost + square_sum, square_sum, ccost

	def compute_ModelConstrainCost(self, path):
		# 暂时不用
		# from pathold to path
		if self.CCD_flag:
			rob_polygons = self.update_rob_polygons(self.path)
			for i in range(0, self.N-1):
				for j in range(len(self.ob_polygons)):
					diffpi = path[i,:] - self.path[i,:]
					diffpi1 = path[i+1,:] - self.path[i+1,:]
					vec = self.rob_ob_infos[i][j][0]
					dis = self.rob_ob_infos[i][j][1]
					p0 = rob_polygons[i].get_far_proj_vertex(vec) # TODO vec or -vec
					p1 = rob_polygons[i].get_far_proj_vertex(vec)
			return 0
		else:
			obstacle_sum = 0
			for i in range(1, self.N-1):
				for j in range(len(self.ob_polygons)):
					diff_p = path[i,:] - self.path[i,:]
					vec = self.rob_ob_infos[i][j][0]
					dis = self.rob_ob_infos[i][j][1]
					sdnew = dis + diff_p.dot(vec)/dis
					obstacle_sum = obstacle_sum + np.max([self.dsafe - sdnew, 0])
			return obstacle_sum

	def compute_ModelCost(self, path):
		aa=(path[0:-1,:]-path[1:,:]).reshape(-1,1)[:,0]
		square_sum = aa.dot(aa)
		# return self.mu * self.compute_ModelConstrainCost(path) + square_sum
		curmid = path[1:-1,:].reshape(-1,1)[:,0]
		return curmid.dot(self.Q_opt).dot(curmid)/2.0 + self.p_opt.dot(curmid) + square_sum
		
	def qp_opt(self, Q, p, G, h):
		# use cvxopt QP solvers
		solvers.options['show_progress'] = False
		sol = solvers.qp(matrix(Q), matrix(p), matrix(G), matrix(h))
		return np.array(sol['x']).reshape(-1, 2)

	def get_collision_rob(self, path, selfpath):
		# 获取当前路径中与障碍碰撞的rob_plg序号和穿透方向（障碍指向rob）
		obstacle_sum = 0
		if selfpath:
			rob_ob_infos = self.rob_ob_infos
		else:
			rob_ob_infos = self.update_gjkepa_info(path)

		collision_robs= []
		m = len(rob_ob_infos)
		n = len(rob_ob_infos[0])
		for i in range(m):
			for j in range(n):
				if rob_ob_infos[i][j][1]<self.dsafe:
					collision_robs.append([i, rob_ob_infos[i][j][0]/float(rob_ob_infos[i][j][1])])

		return collision_robs

	def ConvexifyProblem(self):
		# update p. Q, G are not changed
		# 根据rob_ob_infos中的信息，更新p_opt，该值决定优化过程中碰撞的rob_plg向哪个方向运动可以规避碰撞（减小穿透距离）
		p_opt = np.zeros(2*(self.N-2))
		p_opt_r_vecs = [[] for ii in range(self.N-2)] # 存储N-2个中间路径点穿透方向，应沿反方向移动，最后把所有方向叠加
		p_opt_r_vecs_ob_id = [[] for ii in range(self.N-2)] # 存储每个中间路径点与之相碰撞的障碍物序号
		if self.CCD_flag:
			# TODO
			rob_polygons = self.update_rob_polygons(self.path)
			rob_swept_polygons = self.update_rob_swept_polygons(self.path) 
			for i in range(0, self.N-1):
				for j in range(len(self.ob_polygons)):
					vec = self.rob_ob_infos[i][j][0]
					dis = self.rob_ob_infos[i][j][1]
					if dis < dsafe:
						p_on_rob = self.rob_ob_infos[i][j][2] # 穿透向量在rob_swept上的点
						p_on_ob = self.rob_ob_infos[i][j][3] # 穿透向量在障碍上的点
						p0 = rob_polygons[i].get_far_proj_vertex(vec) # rob_swept前的路径点多边形上的点
						p1 = rob_polygons[i+1].get_far_proj_vertex(vec) #rob_swept后的路径点多边形上的点
						p0p_dis = np.linalg.norm(p0-p_on_rob)
						p1p_dis = np.linalg.norm(p1-p_on_rob)
						alpha = p1p_dis/(p0p_dis+p1p_dis) # p_on_rob = a*p0+(1-a)*p1
						if i > 0: # 同时更新前一个点和后一个点的穿透向量，不更新路径起点 
							p_opt_r_vecs[i-1].append( -alpha*vec/float(dis) ) 
							p_opt_r_vecs_ob_id[i-1].append(j)
						if i < self.N-2: # 同时更新前一个点和后一个点的穿透向量，不更新路径终点
							p_opt_r_vecs[i].append( -(1-alpha)*vec/float(dis) )
							p_opt_r_vecs_ob_id[i].append(j)
			# 若某个中间路径点有多个移动方向，且方向相反，合方向较小，会导致该点没法移动。当点附近路径穿过障碍时出现该情况
			# 方案为尝试沿该移动方向的垂直方向正反移动，沿正向移动move_dis1后脱离该障碍，沿反向移动move_dis2后脱离该障碍
			# 检查正反向脱离后是否又与其他障碍物碰撞，选择未于其他障碍碰撞的方向；若都未碰撞，选择移动距离较小一方
			# 该方案相当于改变了穿透方向和穿透距离，在优化中计算成本改善时因需要用到穿透距离，因此需要修改rob_ob_infos中的相应穿透距离
			# 需要修改的rob_ob_infos中的索引记录在manual_epa_dirs中：记录内容为修改的单位化穿透向量
			# 注意：重新调用该函数ConvexifyProblem()时需用原本的rob_ob_infos，重置manual_epa_dirs为空，再update_gjkepa_info()
			for i in range(self.N-2):
				if len(p_opt_r_vecs[i])>0:
					p_opt_r_vecs_i_comb = np.sum(np.array(p_opt_r_vecs[i]), axis = 0)
					# 判断是否有多个移动方向，且方向相反，合方向较小
					if np.linalg.norm(p_opt_r_vecs_i_comb) < 0.1 and len(p_opt_r_vecs[i])>1: #has at least two vecs with opposite directions
						print("update the coefficient of %dth path point: "%(i+1), self.path[i+1,:])
						# 选择最大的移动方向
						vecinorm = [np.linalg.norm(veci) for veci in p_opt_r_vecs[i]]
						max_id = vecinorm.index(max(vecinorm))
						maxnorm_vec = p_opt_r_vecs[i][max_id]
						# try to move along the direction perpendicular to maxnorm_vec
						vec_vertical = np.array([-maxnorm_vec[1], maxnorm_vec[0]]) # 与最大移动方向垂直的单位方向
						vec_vertical = vec_vertical/np.linalg.norm(vec_vertical) # 与最大移动方向垂直的单位方向
						# calculate two ends along vec_vertical of the obstacle
						this_ob_plg = self.ob_polygons[p_opt_r_vecs_ob_id[i][max_id]] #最大移动方向对应障碍物
						rob_swept_plg = rob_swept_polygons[i+max_id] #最大移动方向对应的rob_swept_plg
						p_on_ob = self.rob_ob_infos[i+max_id][p_opt_r_vecs_ob_id[i][max_id]][3] #contact point on obstacle
						ob_plg_end1 = this_ob_plg.get_far_proj_vertex(vec_vertical) # 沿正垂直方向障碍物的远端点
						ob_plg_end2 = this_ob_plg.get_far_proj_vertex(-vec_vertical) # 沿负垂直方向障碍物的远端点
						# 两个远端点到当前路径点的距离
						dis_to_end1 = ((ob_plg_end1-self.path[i+1,:]).dot(vec_vertical)) #should be positive without abs
						dis_to_end2 = ((ob_plg_end2-self.path[i+1,:]).dot(-vec_vertical)) #should be positive without abs
						print("this two should be positive: ", dis_to_end1, dis_to_end2)
						# 正向移动直至与this_ob_plg不再碰撞
						collision_on_end1 = True
						print("move the %dth rob_swept_plg along direction:"%(i+max_id), vec_vertical)
						move_dis1 = dis_to_end1/2.0
						while collision_on_end1:
							rob_swept_polygons_temp1 = move_along_dir_by_dis(rob_swept_plg, vec_vertical, move_dis1)
							gjkepa = GJK_EPA(rob_swept_polygons_temp1, this_ob_plg)
							if gjkepa.cal_distance()[1]<0: # vector from ob to rob and norm of this vec
								move_dis1 = move_dis1 + dis_to_end1/10.0
							else:
								collision_on_end1 = False
						print("after moving move_dis1=%f, the %dth rob_swept_plg is not in collision with this_ob_plg"%(move_dis1, i+max_id))
						# 反向移动直至与this_ob_plg不再碰撞
						print("move the %dth rob_swept_plg along direction:"%(i+max_id), -vec_vertical)
						collision_on_end2 = True
						move_dis2 = dis_to_end2/2.0
						while collision_on_end2:
							rob_swept_polygons_temp2 = move_along_dir_by_dis(rob_swept_plg, -vec_vertical, move_dis2)
							gjkepa = GJK_EPA(rob_swept_polygons_temp2, this_ob_plg)
							if gjkepa.cal_distance()[1]<0: 
								move_dis2 = move_dis2 + dis_to_end2/10.0
							else:
								collision_on_end2 = False
						print("after moving move_dis2=%f, the %dth rob_swept_plg is not in collision with this_ob_plg"%(move_dis2, i+max_id))
						# 检查移动后是否与其他障碍碰撞
						# check if rob_swept_polygons_temp1 and rob_swept_polygons_temp2 is in collision with other obstacles,
						# select the one that is collision-free with all other obstacles
						rob_temp1_collision_free = True
						rob_temp2_collision_free = True
						for ob_plg in self.ob_polygons:
							if rob_temp1_collision_free:
								gjkepa = GJK_EPA(rob_swept_polygons_temp1, ob_plg)
								if gjkepa.cal_distance()[1]<0:
									rob_temp1_collision_free = False
							if rob_temp2_collision_free:
								gjkepa = GJK_EPA(rob_swept_polygons_temp2, ob_plg)
								if gjkepa.cal_distance()[1]<0:
									rob_temp2_collision_free = False
						# 选择一个方向作为新的穿透方向
						# both are collision free, select the direction with smaller dis
						if rob_temp1_collision_free and rob_temp2_collision_free:
							if move_dis1 < move_dis2:
								p_opt_r_vecs_i_comb = -vec_vertical*move_dis1
								self.manual_epa_dirs[i][p_opt_r_vecs_ob_id[i][max_id]] = [-vec_vertical]
								self.manual_epa_dirs[i+1][p_opt_r_vecs_ob_id[i][max_id]] = [-vec_vertical]
							else:
								p_opt_r_vecs_i_comb = vec_vertical*move_dis2
								self.manual_epa_dirs[i][p_opt_r_vecs_ob_id[i][max_id]] = [vec_vertical]
								self.manual_epa_dirs[i+1][p_opt_r_vecs_ob_id[i][max_id]] = [vec_vertical]
						# temp2, if move along -vec_vertical, collision with other objects
						elif rob_temp1_collision_free and (not rob_temp2_collision_free):
							p_opt_r_vecs_i_comb = -vec_vertical*move_dis1
							self.manual_epa_dirs[i][p_opt_r_vecs_ob_id[i][max_id]] = [-vec_vertical]
							self.manual_epa_dirs[i+1][p_opt_r_vecs_ob_id[i][max_id]] = [-vec_vertical]
						elif (not rob_temp1_collision_free) and rob_temp2_collision_free:
							p_opt_r_vecs_i_comb = vec_vertical*move_dis2
							self.manual_epa_dirs[i][p_opt_r_vecs_ob_id[i][max_id]] = [vec_vertical]
							self.manual_epa_dirs[i+1][p_opt_r_vecs_ob_id[i][max_id]] = [vec_vertical]
						elif (not rob_temp1_collision_free) and (not rob_temp2_collision_free):
							# TODO： 双向移动都不行！！死胡同，算法找不到路径，强制退出算了
							print("WARNING！！！！ cannot find the right optimization direction for the  %dth point:"%(i+1), self.path[i+1,:] )
						print("modified coefficient is: ", p_opt_r_vecs_i_comb)

					# update p_opt
					p_opt[2*i] = self.mu*p_opt_r_vecs_i_comb[0]
					p_opt[2*i+1] = self.mu*p_opt_r_vecs_i_comb[1]
			p_opt[0] = p_opt[0] - 2*self.path[0,0]
			p_opt[1] = p_opt[1] - 2*self.path[0,1]
			p_opt[-2] = p_opt[-2] - 2*self.path[-1,0]
			p_opt[-1] = p_opt[-1] - 2*self.path[-1,1]
		else:	
			for i in range(1, self.N-1):
				x_coeff = 0
				y_coeff = 0
				for j in range(len(self.ob_polygons)):
					vec = self.rob_ob_infos[i][j][0] # vector from ob to rob
					dis = self.rob_ob_infos[i][j][1] #
					# TODO dis>dsafe; dis<dsafe dis=0??
					if dis < dsafe:
						x_coeff = x_coeff - vec[0]/np.sign(dis)
						y_coeff = y_coeff - vec[1]/np.sign(dis)
						# print("%dth points has negative dis %f, vec is "%(i, dis), vec, "point is: ", self.path[i,:])
				p_opt[2*(i-1)] = self.mu * x_coeff
				p_opt[2*(i-1)+1] = self.mu * y_coeff
			p_opt[0] = p_opt[0] - 2*self.path[0,0]
			p_opt[1] = p_opt[1] - 2*self.path[0,1]
			p_opt[-2] = p_opt[-2] - 2*self.path[-1,0]
			p_opt[-1] = p_opt[-1] - 2*self.path[-1,1]
		
		self.p_opt = p_opt
		
	def update_cvxopt_h(self):
		# update h
		h_opt1 = self.path[1:-1].reshape(-1,1) + self.s
		h_opt2 = self.path[1:-1].reshape(-1,1)*(-1) + self.s
		self.h_opt = np.vstack((h_opt1, h_opt2))

	def start_plan(self, line, fig):
		self.rob_ob_infos = self.update_gjkepa_info(self.path)
		total_ittimes = 0
		find_solution = False
		while self.mu < self.mutol:
			convexify_i = 0
			s_too_low = False
			findfreepath = False
			while (not s_too_low) and (not findfreepath) and True:
				# print("convexify_i = ", convexify_i)
				convexify_i = convexify_i+1
				self.manual_epa_dirs = [[[]]*len(self.ob_polygons) for ii in range(self.N-1)]
				self.rob_ob_infos = self.update_gjkepa_info(self.path)
				self.ConvexifyProblem() # update Q, p, G, h 

				# manual_epa_dirs maybe changed after this, update rob_ob_infos using latest manual_epa_dirs
				self.rob_ob_infos = self.update_gjkepa_info(self.path)

				cur_collision_robs = self.get_collision_rob(self.path, True)
				cur_collision_robs_ids = [ccur[0] for ccur in cur_collision_robs]
				print("current path has collision on points: ", cur_collision_robs_ids)
				trust_i = 0
				last_true_cost, last_smcost, last_ccost = self.compute_PathCost(self.path, True)
				last_model_cost = self.compute_ModelCost(self.path)
				while True:
					total_ittimes = total_ittimes + 1
					# print("trust_i = ", trust_i)
					trust_i = trust_i+1
					self.update_cvxopt_h() # update Q, p, G, h
					path = self.qp_opt(self.Q_opt, self.p_opt, self.G_opt, self.h_opt)
					# np.savetxt("Q_opt.txt", self.Q_opt, fmt="%.6f")
					# np.savetxt("p_opt.txt", self.p_opt, fmt="%.6f")
					# np.savetxt("G_opt.txt", self.G_opt, fmt="%.6f")
					# np.savetxt("h_opt.txt", self.h_opt, fmt="%.6f")
					# np.savetxt("curpath.txt", self.path, fmt="%.6f")
					temppath = self.path.copy() # deep copy
					temppath[1:-1,:] = path
					# np.savetxt("temppath.txt", temppath, fmt="%.6f")

					sol_collision_robs = self.get_collision_rob(temppath, False);
					sol_collision_robs_ids = [ssol[0] for ssol in sol_collision_robs]
					print("optimized path has collision on points: ", sol_collision_robs_ids)
					################## trick one: check if find a collision-free path without calculating the cost, if true, exit
					if len(sol_collision_robs) == 0 and self.trick1:
						print("no collision on this optimized path")
						self.path[1:-1,:] = path
						self.manual_epa_dirs = [[[]]*len(self.ob_polygons) for ii in range(self.N-1)]
						self.rob_ob_infos = self.update_gjkepa_info(self.path)
						findfreepath = True
						break
					########
					################## trick two: # if this path leads to collision in other path points, update these points in p_opt
					if len(sol_collision_robs) != 0 and self.trick2:
						print("update p_opt basing on the path generated by cvxopt", len(sol_collision_robs))
						cur_collision_robs_ids = [ccur[0] for ccur in cur_collision_robs]
						for solcr in sol_collision_robs:
							robi = solcr[0]
							if robi not in cur_collision_robs_ids:
								cur_collision_robs.append(solcr)
								self.p_opt[2*(robi-1)] = self.p_opt[2*(robi-1)] - self.mu*solcr[1][0]
								self.p_opt[2*(robi-1)+1] = self.p_opt[2*(robi-1)+1] - self.mu*solcr[1][0]
					########

					# update
					cur_true_cost, cur_smcost, cur_ccost = self.compute_PathCost(temppath, False)
					cur_model_cost = self.compute_ModelCost(temppath)
					true_improve = last_true_cost - cur_true_cost
					model_improve = last_model_cost - cur_model_cost

					if true_improve > 0:# and true_improve/model_improve > self.c:
						print("improved, true_improve=%f, model_improve=%f, s=%f"%(true_improve,model_improve, self.s))
						self.s = self.tau_plus * self.s
						self.path[1:-1,:] = path
						self.manual_epa_dirs = [[[]]*len(self.ob_polygons) for ii in range(self.N-1)]
						self.rob_ob_infos = self.update_gjkepa_info(self.path)
						last_true_cost = cur_true_cost
						last_model_cost = cur_model_cost
						break
					else:
						# print("no improvement, true_improve=%f, model_improve=%f, s=%f"%(true_improve,model_improve, self.s))
						self.s = self.tau_minus * self.s

					if self.s < self.xtol:
						s_too_low = True
						break

				line.set_xdata(self.path[:,0])
				line.set_ydata(self.path[:,1])
				fig.canvas.draw()
				fig.canvas.flush_events()

				if (true_improve>0) and (true_improve < self.ftol):
					break

			s_too_low = False
			self.manual_epa_dirs = [[[]]*len(self.ob_polygons) for ii in range(self.N-1)]
			self.rob_ob_infos = self.update_gjkepa_info(self.path)
			true_cost, smcost, ccost = self.compute_PathCost(self.path, True)
			print("total cost(%f) = smcost(%f) + %f*ccost(%f)"%(true_cost, smcost, self.mu, ccost))
			################## trick three: use dsafe=0 to check the collision states
			ccost1 = self.compute_ConstrainCost(self.path, self.trick3, True)
			#######
			if self.trick3:
				print("ccost with dsafe=0 is: ", ccost1)

			if ccost1 < self.ctol:
				find_solution = True
				break
			else:
				print("enlarge mu from %f to %f"%(self.mu, self.k*self.mu))
				self.mu = self.k * self.mu
				self.s = self.s0
		if find_solution:
			print("!!!! SUCCEED !!!!")
		else:
			print("!!!! FAIL !!!!")
		print("trick1=", self.trick1, ", trick2=", self.trick2, ", trick3=", self.trick3)
		print("total iteration times: ", total_ittimes)
		# keep the figure on
		while True:
			line.set_xdata(self.path[:,0])
			line.set_ydata(self.path[:,1])
			fig.canvas.draw()
			fig.canvas.flush_events()

if __name__ == '__main__':
	map = mymap.Map([[]])
	map.load_map_file('../common/map1.txt')

	(start_node , target_node) = map.generate_start_target()
	
	'''parameters for TrajOpt
	@param side_length: each path point is represented by a square with side_length
	@param N: total number of path point
	@param init_traj: how to set the initial trajectory: 'lin' or 'sin'
	@param dsafe: safe distance when checking collision
	@param CCD_flag: if enable continuous-time collision detection
	@param mu0: initial penalty coefficient
	@param s0: initial trust region size
	@param c: step acceptance parameter
	@param tau_plus: trust region expansion factor
	@param tau_minus: trust region shrinkage factor
	@param k: penalty scaling factor
	@param ftol: convergence thresholds for merit
	@param xtol: convergence thresholds for x
	@param ctol: constraint satisfaction threshold
	'''
	side_length = 1.9;  N = 40; init_traj = "sin"; dsafe = 0.5; CCD_flag = True; 
	trick1_flag = True; trick2_flag = False; trick3_flag = True;
	mu0 = 50; s0 = 2; c = 0.1; 
	tau_plus = 1.2; tau_minus = 0.8; k = 2; ftol = 0.00001; xtol = 0.5; ctol = 0.01;
	TrajOptPlanner = TrajOpt(map, start_node, target_node, side_length, N, init_traj, dsafe, CCD_flag, trick1_flag, trick2_flag, trick3_flag, mu0, s0, c, tau_plus, tau_minus, k, ftol, xtol, ctol)

	######### start the planner
	print("start the planner")
	plt.ion()
	fig, ax = plt.subplots()
	map.draw_map(fig, ax, 'TrajOpt Demo')

	line, = ax.plot(TrajOptPlanner.path[:,0], TrajOptPlanner.path[:,1], 'r.-', linewidth=2)
	fig.canvas.draw()
	fig.canvas.flush_events()
	TrajOptPlanner.start_plan(line, fig)

	######### for test
	# fig, ax = plt.subplots()
	# rob_polygons = TrajOptPlanner.update_rob_polygons(TrajOptPlanner.path)
	# rob_swept_polygons = TrajOptPlanner.update_rob_swept_polygons(TrajOptPlanner.path)
	# TrajOptPlanner.rob_ob_infos = TrajOptPlanner.update_gjkepa_info(TrajOptPlanner.path)
	# TrajOptPlanner.ConvexifyProblem()

	# obi=5
	# TrajOptPlanner.ob_polygons[obi].plot_in_ax(ax)
	# robi = 20
	# print("p_id=", robi+1, TrajOptPlanner.p_opt[2*(robi-1)+2], TrajOptPlanner.p_opt[2*(robi-1)+3])
	# print("p_id=", robi+2, TrajOptPlanner.p_opt[2*(robi-1)+4], TrajOptPlanner.p_opt[2*(robi-1)+5])

	# for i in range(3):
	# 	rob_swept_polygons[robi+i].plot_in_ax(ax)
	# 	plt.plot(TrajOptPlanner.path[robi+i, 0], TrajOptPlanner.path[robi+i, 1], 'ko')
	# 	plt.plot(TrajOptPlanner.path[robi+i+1, 0], TrajOptPlanner.path[robi+i+1, 1], 'ko')
	# 	gjkepa = GJK_EPA(rob_swept_polygons[robi+i], TrajOptPlanner.ob_polygons[obi])
	# 	print(gjkepa.cal_distance()) 

	# plt.axis("Equal")
	# plt.show()


	


