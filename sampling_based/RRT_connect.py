import math
import random
import matplotlib.pyplot as plt
import RRT as rrt
import common.my_map as mymap

def connect(map , tree , node , step_length ,collision_check_step):
    while True:
        print("connect extend")
        (new_node , result) = rrt.extend(map,tree,node,step_length,collision_check_step)
        if result == 0:
            return new_node , 0    #trap
        else:
            if new_node.distance(node) < 3:
                #node.parent_node = new_node
                #tree.append(node)
                #plt.plot([node.x, node.x], [new_node.y, new_node.y], 'b')
                return new_node , 1    #reach


def RRT_connect(map , start_node , target_node , max_point_num , step_length ,collision_check_step):
    tree_a = [start_node]
    tree_b = [target_node]

    for i in range(max_point_num):
        random_node = map.generate_random_node()
        (new_node, result) = rrt.extend(map, tree_a, random_node, step_length, collision_check_step)
        if result == 0:
            continue
        (reach_node, result) = connect(map, tree_b, new_node, step_length, collision_check_step)
        if result == 1:
            path = []
            node = reach_node
            path.append(node)
            while node != tree_b[0] and not node.parent_node is None:
                node = node.parent_node
                path.append(node)
            path = path[::-1]

            node = new_node
            while node != tree_a[0] and not node.parent_node is None:
                node = node.parent_node
                path.append(node)
            return path
        tree_a , tree_b = tree_b , tree_a


if __name__ == '__main__':
    map = mymap.Map([[]])
    map.load_map_file('../common/map2.txt')
    map.draw_map()
    (start_node , target_node) = map.generate_start_target()
    path = RRT_connect(map , start_node , target_node , 10000 , 5 , 0.5)
    print("path len:",len(path))
    map.draw_path(path)
    plt.ioff()
    plt.show()