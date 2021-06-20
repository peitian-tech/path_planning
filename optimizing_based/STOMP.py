# -*- coding: utf-8 -*-
import sys
sys.path.append("polygon/") 
import my_map as mymap
import numpy as np
import matplotlib.pyplot as plt

def signed_edt(point, env_map):
	# if this point is occupied 
	p_status = env_map.if_on_obstacle(point[0], point[1]) 
	if not p_status:
		return 0
	# calculate the Manhattan Distance to the nearest point that has different type of this point
	m_dist = 0
	find_p = False
	nearest_p = np.zeros(2)
	while (not find_p):
		m_dist = m_dist + 1
		for i in range(0, m_dist+1):
			nearest_p[0] = point[0] + i
			nearest_p[1] = point[1] + m_dist-i
			if nearest_p[0]<0 or nearest_p[0]>=100 or nearest_p[1]<0 or nearest_p[1]>=100:
				npstatus = True
			else:
				npstatus = env_map.if_on_obstacle(nearest_p[0], nearest_p[1])
			if npstatus != p_status:
				find_p = True
				break
			nearest_p[0] = point[0] - i
			nearest_p[1] = point[1] + m_dist-i
			if nearest_p[0]<0 or nearest_p[0]>=100 or nearest_p[1]<0 or nearest_p[1]>=100:
				npstatus = True
			else:
				npstatus = env_map.if_on_obstacle(nearest_p[0], nearest_p[1])
			if npstatus != p_status:
				find_p = True
				break
			nearest_p[0] = point[0] + i
			nearest_p[1] = point[1] - (m_dist-i)
			if nearest_p[0]<0 or nearest_p[0]>=100 or nearest_p[1]<0 or nearest_p[1]>=100:
				npstatus = True
			else:
				npstatus = env_map.if_on_obstacle(nearest_p[0], nearest_p[1])
			if npstatus != p_status:
				find_p = True
				break
			nearest_p[0] = point[0] - i
			nearest_p[1] = point[1] - (m_dist-i)
			if nearest_p[0]<0 or nearest_p[0]>=100 or nearest_p[1]<0 or nearest_p[1]>=100:
				npstatus = True
			else:
				npstatus = env_map.if_on_obstacle(nearest_p[0], nearest_p[1])
			if npstatus != p_status:
				find_p = True
				break

	# return the Euclidean distance
	if p_status:
		return -1*np.linalg.norm(point-nearest_p)
	else:
		return np.linalg.norm(point-nearest_p)

def compute_PathPointCost(point, env_map, dsafe, rb, vel):
	return (dsafe + rb -signed_edt(point, env_map))*vel

def compute_ConstrainCost(path, env_map, dsafe, rb):
	edtsum = 0
	for i in range(1, np.size(path, 0)):
		vel = np.linalg.norm(path[i,:]-path[i-1,:])
		edtsum = edtsum + compute_PathPointCost(path[i,:], env_map, dsafe, rb, vel)
	return edtsum

def compute_PathCost(path, R, env_map, dsafe, rb):
	return compute_ConstrainCost(path, env_map, dsafe, rb)# + 0.001*0.5*np.sum(path.transpose().dot(R).dot(path))

def generate_PathNoise(Rinv):
	N = np.size(Rinv, 0)
	gaussout=np.random.multivariate_normal(np.zeros(N),Rinv,size=(1,2))[0].transpose()
	# gaussout=np.random.multivariate_normal(np.zeros(N),np.diag(np.ones(N)*0.036),size=(1,2))[0].transpose()
	return gaussout*100

class STOMP():
	def __init__(self, env_map, start_node, target_node, N, K, init_traj, dsafe, rb, convThr):
		self.env_map = env_map
		self.start_node = start_node
		self.target_node = target_node
		self.N = N
		self.K = K
		self.init_traj = init_traj
		self.dsafe = dsafe
		self.rb = rb
		self.convThr = convThr
		# construct the init trajectory to be optimazied
		path = np.zeros([N, 2])
		if init_traj=='lin':
			dx = (target_node.x - start_node.x)/(N-1.0)
			dy = (target_node.y - start_node.y)/(N-1.0)
			for i in range(N):
				path[i, 0] = start_node.x + i*dx
				path[i, 1] = start_node.y + i*dy
		else:
			print("ERROR: unsupported initial trajectory type!")

		self.path = path

		# precompute
		A = np.zeros([N+2, N])
		A[0, 0] = 1; A[1, 0] = -2; A[1, 1] = 1;
		A[N+1, N-1] = 1; A[N, N-2] = 1; A[N, N-1] = -2;
		for i in range(2, N):
			A[i, i-2] = 1
			A[i, i-1] = -2
			A[i, i] = 1
		R = np.dot(A.transpose(), A)
		Rinv = np.linalg.inv(R)
		M = 1.0/N * np.dot(Rinv, np.diag(1.0/np.max(Rinv, 0)))
		Rinv = Rinv/np.sum(Rinv)
		print(max(np.diag(Rinv)))

		self.A = A
		self.R = R
		self.Rinv = Rinv
		self.M = M

	def start_plan(self, line, fig):
		self.Qcost = compute_PathCost(self.path, self.R, self.env_map, self.dsafe, self.rb)
		print("Qcost: %f"%self.Qcost)
		self.Qcostold = -999
		ittimes = 0
		while self.Qcost > self.convThr: #abs(Qcost - Qcostold) > convThr:
			ittimes = ittimes + 1
			self.Qcostold = self.Qcost

			path_noise = list(range(self.K))

			############### weight * noisy point
			# for k in range(self.K):
			# 	path_noise[k] = generate_PathNoise(self.Rinv)
			# dpath = np.zeros(self.path.shape)
			# for i in range(1, self.N-1):
			# 	noise_pathpoint_cost = np.zeros(self.K)
			# 	noisy_pathpoint_weight = np.zeros(self.K)
			# 	vel = np.linalg.norm(self.path[i,:] - self.path[i-1,:])
			# 	for k in range(self.K):
			# 		noise_pathpoint_cost[k] = compute_PathPointCost(self.path[i,:]+path_noise[k][i,:], self.env_map, self.dsafe, self.rb, vel)
			# 	minS = min(noise_pathpoint_cost)
			# 	maxS = max(noise_pathpoint_cost)
			# 	for k in range(self.K):
			# 		if maxS==minS:
			# 			noisy_pathpoint_weight[k] = 1
			# 		else:
			# 			noisy_pathpoint_weight[k] = np.exp(-10*(noise_pathpoint_cost[k]-minS)/(maxS-minS))
			# 	noisy_pathpoint_weight = noisy_pathpoint_weight/np.sum(noisy_pathpoint_weight)

			# 	for k in range(self.K):
			# 		dpath[i,:] = dpath[i,:] + noisy_pathpoint_weight[k]*path_noise[k][i,:]
			#########

			############### weight*noisy path
			noisy_path_cost = np.zeros(self.K)
			noisy_path_weight = np.zeros(self.K)
			for i in range(self.K):
				path_noise[i] = generate_PathNoise(self.Rinv)
				noisy_path_cost[i] = compute_PathCost(self.path + path_noise[i], self.R, self.env_map, self.dsafe, self.rb)
			minS = min(noisy_path_cost)
			maxS = max(noisy_path_cost)
			for i in range(self.K):
				noisy_path_weight[i] = np.exp(-2*(noisy_path_cost[i]-minS)/(maxS-minS))
			noisy_path_weight = noisy_path_weight/sum(noisy_path_weight)
			dpath = np.zeros(self.path.shape)
			for i in range(self.K):
				dpath = dpath + noisy_path_weight[i]*path_noise[i]
			#######

			dpath = self.M.dot(dpath)*min(self.Qcost, 20)
			# print("dpathmax: %f"%np.max(dpath))
			self.path = self.path+dpath
			self.Qcost = compute_PathCost(self.path, self.R, self.env_map, self.dsafe, self.rb)
			print("Qcost: %f"%self.Qcost)
			line.set_xdata(self.path[:,0])
			line.set_ydata(self.path[:,1])
			fig.canvas.draw()
			fig.canvas.flush_events()

		print("iteration %d times."%ittimes)
		while True:
			line.set_xdata(self.path[:,0])
			line.set_ydata(self.path[:,1])
			fig.canvas.draw()
			fig.canvas.flush_events()

def update_fig(i):
	line.set_xdata(stompplanner.path[:,0])
	line.set_ydata(stompplanner.path[:,1])
	return line, ax

if __name__ == '__main__':
	map = mymap.Map([[]])
	map.load_map_file('map1.txt')

	(start_node , target_node) = map.generate_start_target()
	stompplanner = STOMP(map , start_node , target_node , 100 , 50 , 'lin', 0, 0, 0.0001)

	print("plot")
	plt.ion()
	fig, ax = plt.subplots()
	map.draw_map(fig, ax, 'STOMP demo')

	line, = ax.plot(stompplanner.path[:,0], stompplanner.path[:,1], 'r.-', linewidth=2)
	fig.canvas.draw()
	stompplanner.start_plan(line, fig)




