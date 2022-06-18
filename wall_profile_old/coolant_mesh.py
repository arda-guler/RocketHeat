from coolant_meshing import *
from material import SS304L, CuCrZr
import matplotlib.pyplot as plt

SS304L = SS304L()
CuCrZr = CuCrZr()

# Section 1
#   |---|
#    ______________________________________________________  
#   |   
#   |                        SECTION 3
# C |   
# o |    __________________________________________________
# m | S |
# b | E |
# u | C |
# s | T |
# t | I |
# i | O |                Coolant Channel
# o | N |
# n |   |
#   | 1 |
#   |   |
#   |   |__________________________________________________ _
#   |                                                       |
#   |               SECTION 2                               | Section 2
#   |                                                       |
#   |______________________________________________________ -

T_init = 273 + 25

def create_coolant_channel_mesh(L_x, L_y, cell_L_x, cell_L_y, channel_L_x, channel_L_y):
    n_x = int(L_x/cell_L_x)
    n_y = int(L_y/cell_L_y)

    n_channel_x_start = int((L_x - channel_L_x)/cell_L_x)
    n_channel_y_start = int(((L_y - channel_L_y)/cell_L_y) / 2)
    n_channel_y_end = n_channel_y_start + int(channel_L_y/cell_L_y) - 1
    
    cells = [ [] for i in range(n_y) ]
    cells_x = []
    cells_y = []
    cells_flag = []

    # create Section 1 cells
    for x_i in range(n_channel_x_start):
        for y_i in range(n_y):
            new_cell_size = [cell_L_x, cell_L_y, 0.5]
            new_cell_index = [x_i, y_i, 0]
            new_cell_pos = [x_i*cell_L_x + cell_L_x*0.5, y_i*cell_L_y + cell_L_y*0.5, 0]
            if new_cell_index[0] == 0:
                new_cell = cell(new_cell_index, new_cell_pos, new_cell_size, CuCrZr, T_init, "boundary_chamber")
                cells_flag.append("r")
            elif new_cell_index[0] == n_channel_x_start-1 and n_channel_y_start <= new_cell_index[1] <= n_channel_y_end:
                new_cell = cell(new_cell_index, new_cell_pos, new_cell_size, CuCrZr, T_init, "boundary_coolant")
                cells_flag.append("b")
            else:
                new_cell = cell(new_cell_index, new_cell_pos, new_cell_size, CuCrZr, T_init)
                cells_flag.append("k")
            cells[y_i].insert(x_i, new_cell)
            cells_x.append(new_cell_pos[0])
            cells_y.append(new_cell_pos[1])

    # create Section 2 cells
    for x_i in range(n_channel_x_start, n_x, 1):
        for y_i in range(n_channel_y_start):
            new_cell_size = [cell_L_x, cell_L_y, 0.5]
            new_cell_index = [x_i, y_i, 0]
            new_cell_pos = [x_i*cell_L_x + cell_L_x*0.5, y_i*cell_L_y + cell_L_y*0.5, 0]
            if n_channel_x_start <= new_cell_index[0] and new_cell_index[1] == n_channel_y_start-1:
                new_cell = cell(new_cell_index, new_cell_pos, new_cell_size, CuCrZr, T_init, "boundary_coolant")
                cells_flag.append("b")
            else:
                new_cell = cell(new_cell_index, new_cell_pos, new_cell_size, CuCrZr, T_init)
                cells_flag.append("k")
            cells[y_i].insert(x_i, new_cell)
            cells_x.append(new_cell_pos[0])
            cells_y.append(new_cell_pos[1])

    # create Section 3 cells
    for x_i in range(n_channel_x_start, n_x, 1):
        for y_i in range(n_channel_y_end+1, n_y, 1):
            new_cell_size = [cell_L_x, cell_L_y, 0.5]
            new_cell_index = [x_i, y_i, 0]
            new_cell_pos = [x_i*cell_L_x + cell_L_x*0.5, y_i*cell_L_y + cell_L_y*0.5, 0]
            if n_channel_x_start <= new_cell_index[0] and new_cell_index[1] == n_channel_y_end + 1:
                new_cell = cell(new_cell_index, new_cell_pos, new_cell_size, CuCrZr, T_init, "boundary_coolant")
                cells_flag.append("b")
            else:
                new_cell = cell(new_cell_index, new_cell_pos, new_cell_size, CuCrZr, T_init)
                cells_flag.append("k")
            cells[y_i].insert(x_i, new_cell)
            cells_x.append(new_cell_pos[0])
            cells_y.append(new_cell_pos[1])

    print("Close the 'Mesh and Boundaries' plot to continue.")
    plt.scatter(cells_x, cells_y, c=cells_flag, marker="s")
    plt.title("Mesh and Boundaries")
    plt.show()
    print("Proceeding with analysis...")
    # n_x, n_y, n_z, cells
    clt_mesh = mesh(n_x, n_y, 1, cells)
    return clt_mesh
