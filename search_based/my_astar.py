import numpy as np
import matplotlib.pyplot as plt
import random
import heapq as pq
import common.my_map as mymap

class Astar:
    open_list = []
    close_list = []
    def search_path_of_graph(self, graph , start_node , target_node):
        #起点加入openlist
        pq.heappush(self.open_list , start_node)
        while True:
            current_node = pq.heappop(self.open_list)
            #print("current_node fvalue:", current_node.fvalue)
            #从openlist中取出代价f最小的节点，加到closelist中
            self.close_list.append(current_node)
            #print("adjacnet node num:" , len(current_node.adjacent_node))
            see_target = False
            for i in range(len(current_node.adjacent_node)):
                adja_node = current_node.adjacent_node[i]
                #忽略在closelist中的node
                if adja_node in self.close_list:
                    continue
                #如果它不在 open list 中，把它加入 open list ，并且把当前方格设置为它的父亲，记录该方格的 F ， G 和 H 值
                if adja_node not in self.open_list:
                    #print("add to open list")
                    adja_node.parent_node = current_node
                    adja_node.gvalue = current_node.gvalue + current_node.adjacent_cost[i]
                    adja_node.hvalue = adja_node.distance(target_node)
                    adja_node.update_fvalue()
                    pq.heappush(self.open_list, adja_node)
                    if adja_node is target_node:
                        print("see target.")
                        see_target = True
                        break
                #如果它已经在 open list 中，检查这条路径 ( 即经由当前方格到达它那里 ) 是否更好，用 G 值作参考。更小的 G 值表示这是更好的路径。
                #如果是这样，把它的父亲设置为当前方格，并重新计算它的 G 和 F 值。如果你的 open list 是按 F 值排序的话，改变后你可能需要重新排序
                else:
                    new_gvalue = current_node.gvalue + current_node.adjacent_cost[i]
                    if new_gvalue < adja_node.gvalue:
                        adja_node.parent_node = current_node
                        adja_node.gvalue = new_gvalue
                        adja_node.update_fvalue()
                        pq.heapify(self.open_list)

            if len(self.open_list)==0:
                print("fail.")
                break
            if see_target == True:
                print("success.")
                break
        path = []
        node = target_node
        path.append(node)
        while node != start_node and not node.parent_node is None:
            node = node.parent_node
            path.append(node)
        return path

    def search_path_of_map(self, map):
        (graph , start_node , target_node) = map.convert_to_graph()
        return self.search_path_of_graph(graph , start_node , target_node)


if __name__ == '__main__':
    map = mymap.Map([[]])
    map.load_map_file('../common/map1.txt')
    map.draw_map()
    astar = Astar()
    path = astar.search_path_of_map(map)
    map.draw_path(path)
    plt.ioff()
    plt.show()
