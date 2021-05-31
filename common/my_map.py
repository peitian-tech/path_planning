import matplotlib.pyplot as plt
from pylab import *
mpl.rcParams['font.sans-serif'] = ['SimHei']

class Node:
    # 构造函数
    def __init__(self, type, x, y):
        # 节点类型
        self.type = type
        # 节点坐标
        self.x = x
        self.y = y
        # 相邻节点
        self.adjacent_node = []
        # 相邻节点代价值
        self.adjacent_cost = []
        # 父节点索引
        self.parent_node = None
        # f = g + h
        self.gvalue = 0
        self.hvalue = 0
        self.fvalue = 0
    def __lt__(self, other):
        return self.fvalue < other.fvalue
    # 添加一个相邻节点
    def add_adjacent_node(self, node, cost):
        self.adjacent_node.append(node)
        self.adjacent_cost.append(cost)
    def update_fvalue(self):
        self.fvalue = self.gvalue + self.hvalue
    def distance(self , other):
        return sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

class Map:
    # 用二维列表表示地图，0：可通行 1：障碍 2：起点 3：目标点
    # 构造函数
    def __init__(self, map_grid):
        self.init_para(map_grid)

    def init_para(self, map_grid):
        self.map_grid = map_grid
        self.row_num = len(self.map_grid)
        self.column_num = len(self.map_grid[0])
        self.grid_width = 100 / max(self.column_num, self.row_num)

    def load_map_file(self, map_file):
        self.init_para(np.loadtxt(map_file, dtype=int))

    # 画地图
    def draw_map(self):
        plt.ion()
        fig = plt.figure()
        fig.canvas.mpl_connect('button_press_event', self.OnClick)
        self.btn_start = Button(plt.axes([0.3, 0.01, 0.07, 0.05]), 'start', color='khaki', hovercolor='yellow')
        ax = fig.add_subplot(111)
        plt.title('路径规划算法演示')
        ax.set_aspect('equal')
        ax.set_xlim([0, 100])
        ax.set_ylim([0, 100])
        for i in range(self.row_num):
            for j in range(self.column_num):
                color = 'black'
                if self.map_grid[i][j] == 0:
                    color = 'green'
                elif self.map_grid[i][j] == 1:
                    color = 'black'
                elif self.map_grid[i][j] == 2:
                    color = 'red'
                elif self.map_grid[i][j] == 3:
                    color = 'blue'
                else:
                    print("error!")
                ax.add_patch(
                    plt.Rectangle((self.grid_width * j, self.grid_width * (self.row_num - 1 - i)), self.grid_width,
                                  self.grid_width, linewidth=1,
                                  edgecolor='gray', facecolor=color))
        #plt.show(block=False)
        #plt.show()
        fig.canvas.draw()

    # 画路径
    def draw_path(self, path):
        for i in range(len(path)):
            if i < len(path) - 1:
                plt.plot([path[i].x, path[i + 1].x], [path[i].y, path[i + 1].y],'r')
        #plt.show(block=False)
        plt.gcf().canvas.draw()

    # 将地图转为图数据结构
    def convert_to_graph(self):
        graph_node_list = []
        start_node = None
        target_node = None
        for i in range(self.row_num):
            graph_row_list = []
            for j in range(self.column_num):
                type = self.map_grid[i][j]
                node = Node(type, self.grid_width * j + self.grid_width / 2,
                            self.grid_width * (self.row_num - 1 - i) + self.grid_width / 2)
                graph_row_list.append(node)
            graph_node_list.append(graph_row_list)
        graph = []
        for i in range(self.row_num):
            for j in range(self.column_num):
                node = graph_node_list[i][j]
                type = node.type
                if type == 1:
                    continue
                if i > 0:
                    if self.map_grid[i - 1][j] != 1:
                        node.add_adjacent_node(graph_node_list[i - 1][j], self.grid_width)
                    if j > 0 and self.map_grid[i - 1][j - 1] != 1:
                        node.add_adjacent_node(graph_node_list[i - 1][j - 1], sqrt(2) * self.grid_width)
                    if j < self.column_num - 1 and self.map_grid[i - 1][j + 1] != 1:
                        node.add_adjacent_node(graph_node_list[i - 1][j + 1], sqrt(2) * self.grid_width)
                if j > 0 and self.map_grid[i][j - 1] != 1:
                    node.add_adjacent_node(graph_node_list[i][j - 1], self.grid_width)
                if j < self.column_num - 1 and self.map_grid[i][j + 1] != 1:
                    node.add_adjacent_node(graph_node_list[i][j + 1], self.grid_width)
                if i < self.row_num - 1:
                    if self.map_grid[i + 1][j] != 1:
                        node.add_adjacent_node(graph_node_list[i + 1][j], self.grid_width)
                    if j > 0 and self.map_grid[i + 1][j - 1] != 1:
                        node.add_adjacent_node(graph_node_list[i + 1][j - 1], sqrt(2) * self.grid_width)
                    if j < self.column_num - 1 and self.map_grid[i + 1][j + 1] != 1:
                        node.add_adjacent_node(graph_node_list[i + 1][j + 1], sqrt(2) * self.grid_width)
                if type == 2:
                    start_node = node
                elif type == 3:
                    target_node = node
                graph.append(node)
        print("graph size:", len(graph))
        return graph , start_node , target_node
    def generate_random_node(self):
        x = self.row_num * self.grid_width * random()
        y = self.column_num * self.grid_width * random()
        node = Node(0, x, y)
        return node
    #地图内随机生成一个不在障碍物上的点，并构建node返回
    def generate_random_valid_node(self):
        while True:
            x = self.row_num * self.grid_width * random()
            y = self.column_num * self.grid_width * random()
            if not self.if_on_obstacle(x , y):
                break
        node = Node(0,x,y)
        return node

    def find_inter_node(self , start, target, step_length):
        dis = start.distance(target)
        if dis == 0:
            return start
        x = start.x + (target.x - start.x) / dis * step_length
        y = start.y + (target.y - start.y) / dis * step_length
        return Node(0, x, y)

    #生成起点和目标点的Node
    def generate_start_target(self):
        start_node = Node(2,0,0)
        target_node = Node(3,0,0)
        for i in range(self.row_num):
            for j in range(self.column_num):
                type = self.map_grid[i][j]
                if type == 2:
                    start_node = Node(type, self.grid_width * j + self.grid_width / 2,
                            self.grid_width * (self.row_num - 1 - i) + self.grid_width / 2)
                elif type == 3:
                    target_node = Node(type, self.grid_width * j + self.grid_width / 2,
                            self.grid_width * (self.row_num - 1 - i) + self.grid_width / 2)
        return start_node , target_node

    #判断两个节点之间是否经过障碍物
    def if_cross_obstacle(self , a , b , step):
        dis= a.distance(b)
        if dis == 0:
            return False
        step_x = (b.x - a.x) / dis * step
        step_y = (b.y - a.y) / dis * step
        point_num = int(dis / step)
        for i in range(point_num):
            if self.if_on_obstacle(a.x + step_x * i , a.y + step_y * i):
                return True
        return False

    #判断一个点是否在障碍物上
    def if_on_obstacle(self , x , y):
        column = int(x / self.grid_width)
        row = self.row_num - 1 - int(y / self.grid_width)
        return self.map_grid[row][column] == 1

    def set_start_action(self , action):
        self.btn_start.on_clicked(action)

    def OnClick(self ,event):
        column = int(event.xdata / self.grid_width)
        row = self.row_num - 1 - int(event.ydata / self.grid_width)
        color=''
        if event.button == 1: #mouse left button,add obstacle
            print("left:", event.xdata, event.ydata)
            self.map_grid[row][column] = 1
            color='black'
        if event.button == 3: #mouse right button,delete obstacle
            print("right:", event.xdata, event.ydata)
            self.map_grid[row][column] = 0
            color='green'

        plt.gca().add_patch(
            plt.Rectangle((self.grid_width * column, self.grid_width * (self.row_num - 1 - row)), self.grid_width,
                          self.grid_width, linewidth=1, edgecolor='gray', facecolor=color))
        plt.gcf().canvas.draw()