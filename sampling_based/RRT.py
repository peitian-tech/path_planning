import math
import random
import matplotlib.pyplot as plt
import common.my_map as mymap
#查找距离cnode最近的node
def search_nearest_node(tree, cnode):
    nearest_node = tree[0]
    nearest_distance = cnode.distance(tree[0])
    for node in tree:
        if math.fabs(node.x - cnode.x) >= nearest_distance:
            continue
        if math.fabs(node.y - cnode.y) >= nearest_distance:
            continue
        dis = node.distance(cnode)
        if dis < nearest_distance:
            nearest_node = node
            nearest_distance = dis
    return nearest_node

#查找距离cnode小于search_range的所有node
def search_range_node(tree , cnode , search_range):
    nodes = []
    for node in tree:
        if math.fabs(node.x - cnode.x) >= search_range:
            continue
        if math.fabs(node.y - cnode.y) >= search_range:
            continue
        dis = node.distance(cnode)
        if dis < search_range:
            nodes.append(node)
    return nodes

def optimize_tree(map , tree , new_node , search_range ,collision_check_step):
    range_nodes = search_range_node(tree , new_node , search_range)
    #重新选择父节点
    for node in range_nodes:
        if node is not new_node.parent_node:
            new_gvalue = node.gvalue + node.distance(new_node)
            if  new_gvalue < new_node.gvalue and not map.if_cross_obstacle(node , new_node , collision_check_step):
                new_node.parent_node = node
                new_node.gvalue = new_gvalue
    #重新布线随机树
    for node in range_nodes:
        if node is not new_node.parent_node:
            new_gvalue = new_node.gvalue + node.distance(new_node)
            if  new_gvalue < node.gvalue and not map.if_cross_obstacle(node , new_node , collision_check_step):
                node.parent_node = new_node
                node.gvalue = new_gvalue

def extend(map , tree , node , step_length , collision_check_step):
    # 在树中查找离该随机点最近的节点
    near_node = search_nearest_node(tree, node)
    #print("near node:", near_node.x, "  ", near_node.y)
    new_node = map.find_inter_node(near_node, node, step_length)
    #print("new node:", new_node.x, "  ", new_node.y)
    if map.if_cross_obstacle(near_node, new_node, collision_check_step):
        return new_node , 0
    new_node.parent_node = near_node
    new_node.gvalue = near_node.gvalue + near_node.distance(new_node)
    tree.append(new_node)
    optimize_tree(map, tree, new_node, 10, collision_check_step)
    plt.plot([near_node.x, new_node.x], [near_node.y, new_node.y], 'b')
    return new_node , 1

def RRT(map , start_node , target_node , way_point_cache , max_point_num , step_length ,collision_check_step):
    tree = [start_node]
    new_node = start_node
    for i in range(max_point_num):
        # goal bias,随机以目标点作为随机生成点
        rand_value = random.random()
        if rand_value < 0.1:
            random_node = target_node
        elif len(way_point_cache) > 0:
            if rand_value < 0.7:
                rand_index = random.randint(0 , len(way_point_cache)-1)
                random_node = way_point_cache[rand_index]
            else:
                random_node = map.generate_random_node()
        else:
            random_node = map.generate_random_node()
        print(i , "random node:" , random_node.x , "  ",random_node.y)
        (new_node , result) = extend(map, tree, random_node, step_length, collision_check_step)
        #if result is not 0:
        #    plt.plot([near_node.x, new_node.x], [near_node.y, new_node.y], 'b')
        if result is not 0 and new_node.distance(target_node) < 3:
            target_node.parent_node = new_node
            break
        # 加上下面两行后是RRT*
        if result is 1:
           optimize_tree(map, tree, new_node, 10, collision_check_step)

    path = []
    node = target_node
    path.append(node)
    while node != start_node and not node.parent_node is None:
        node = node.parent_node
        path.append(node)
    return path

def start_plan(event):
    print("start_plan")
    (start_node , target_node) = map.generate_start_target()
    path = RRT(map , start_node , target_node , [] , 10000 , 5 , 0.5)
    print("path len:",len(path))
    map.draw_path(path)
   # plt.show()

if __name__ == '__main__':
    map = mymap.Map([[]])
    map.load_map_file('../common/map2.txt')
    map.draw_map()
    map.set_start_action(start_plan)
    (start_node , target_node) = map.generate_start_target()
    path = RRT(map , start_node , target_node , [], 10000 , 5 , 0.5)
    print("path len:",len(path))
    map.draw_path(path)
    plt.ioff()
    plt.show()
