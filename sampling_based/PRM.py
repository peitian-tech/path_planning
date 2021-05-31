import myAstar as astar
import my_map as mymap
import matplotlib.pyplot as plt
import heapq as pq
import numpy as np
import os
def PRM(map , start_node , target_node):
    #随机生成的点数
    point_num = 100
    #选距离最近的k个点连边
    k = 4
    #碰撞检测步长
    collision_check_step = 2
    #固定距离策略用的固定最大距离
    max_distance = 15
    graph=[start_node , target_node]
    # graph_map = np.loadtxt("graph_map.txt");
    # for i in graph_map:
    #     node = astar.Node(0,i[0],i[1])
    #     plt.plot(node.x , node.y ,'yo')
    #     graph.append(node)
    for i in range(point_num):
        node = map.generate_random_valid_node()
        # with open("graph_map.txt","a") as f:
        #     f.write(str(node.x) + " " + str(node.y) + "\n")
        plt.plot(node.x , node.y ,'yo')
        graph.append(node)
    for cnode in graph:
        inode_list = []
        #计算cnode与其他每个节点的距离，选出最小的k个
        for inode in graph:
            if inode is not cnode:
                #元组数组，元素为（代价，节点）
                inode_list.append((cnode.distance(inode) , inode))
        pq.heapify(inode_list)
        #print("inode_list:",inode_list)
        # 策略1，选距离cnode一定距离之内的点
        for node in inode_list:
            if cnode.distance(node[1]) > max_distance:
                continue
            if map.if_cross_obstacle(cnode , node[1] , collision_check_step):
                continue
            cnode.add_adjacent_node(node[1], node[0])
            if cnode not in node[1].adjacent_node:
                node[1].add_adjacent_node(cnode, node[0])
            plt.plot([cnode.x, node[1].x], [cnode.y, node[1].y], 'y')
            plt.text((cnode.x + node[1].x) / 2, (cnode.y + node[1].y) / 2, str(format(node[0], '.1f')))

        # # 策略2，选除了检测到碰撞之外的k个最小距离节点
        # smallest_num = 0
        # while len(inode_list) != 0:
        #     node = pq.heappop(inode_list)
        #     if map.if_cross_obstacle(cnode , node[1] , collision_check_step):
        #         continue
        #     cnode.add_adjacent_node(node[1] ,node[0])
        #     if cnode not in node[1].adjacent_node:
        #         node[1].add_adjacent_node(cnode ,node[0] )
        #     plt.plot([cnode.x , node[1].x] , [cnode.y , node[1].y] , 'y')
        #     plt.text((cnode.x + node[1].x)/2, (cnode.y + node[1].y)/2, str(format(node[0], '.1f')))
        #     smallest_num = smallest_num + 1
        #     if smallest_num == k:
        #         break

        # 策略3，选距离最小的k个节点，并去掉检测到碰撞的节点
        # smallest_node = pq.nsmallest(k,inode_list)
        # #print("smallest_node:",smallest_node)
        # #把这k个作为cnode的临节点，如果node的临节点中没有cnode，也需要把cnode加到node的临节点中
        # for node in smallest_node:
        #     if map.if_cross_obstacle(cnode , node[1] , collision_check_step):
        #         continue
        #     cnode.add_adjacent_node(node[1] ,node[0])
        #     if cnode not in node[1].adjacent_node:
        #         node[1].add_adjacent_node(cnode ,node[0] )
        #     plt.plot([cnode.x , node[1].x] , [cnode.y , node[1].y] , 'y')
        #     plt.text((cnode.x + node[1].x)/2, (cnode.y + node[1].y)/2, str(format(node[0], '.1f')))
    return graph
if __name__ == '__main__':
    map = mymap.Map([[]])
    map.load_map_file('map2.txt')
    map.draw_map()
    (start_node , target_node) = map.generate_start_target()
    graph = PRM(map , start_node , target_node)
    astar =  astar.Astar()
    path = astar.search_path_of_graph(graph , start_node , target_node)
    print("path len:",len(path))
    map.draw_path(path)
    plt.show()