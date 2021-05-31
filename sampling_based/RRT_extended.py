import matplotlib.pyplot as plt
import RRT as rrt
import common.my_map as mymap

def start_plan(event):
    print("start_plan")
    (start_node , target_node) = map.generate_start_target()
    print("way_point_cache:",way_point_cache)
    path = rrt.RRT(map , start_node , target_node , way_point_cache , 10000 , 5 , 0.5)
    print("path len:",len(path))
    map.draw_path(path)
    add_path_to_cache(path, way_point_cache)
   # plt.show()

def add_path_to_cache(path , cache):
    for node in path:
        cache.append(node)
        if len(cache) > 30:
            del cache[0]

if __name__ == '__main__':
    map = mymap.Map([[]])
    map.load_map_file('../common/map1.txt')
    map.draw_map()
    map.set_start_action(start_plan)
    (start_node , target_node) = map.generate_start_target()
    path = rrt.RRT(map , start_node , target_node , [], 10000 , 5 , 0.5)
    print("path len:",len(path))
    map.draw_path(path)

    way_point_cache = []
    add_path_to_cache(path , way_point_cache)
    plt.ioff()
    plt.show()