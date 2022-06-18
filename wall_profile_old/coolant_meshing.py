from material import CuCrZr, SS304L, Jet_A1
CuCrZr = CuCrZr()
SS304L = SS304L()
Jet_A1 = Jet_A1()

class mesh:
    def __init__(self, n_x, n_y, n_z, cells):
        self.n_x = n_x
        self.n_y = n_y
        self.n_z = n_z
        self.cells = cells

    def get_cell(self, index_x, index_y, index_z):
        try:
            return self.cells[index_y][index_x]
        except IndexError:
            return None

# a cell is a 3D cube of material
class cell:
    def __init__(self, index, pos, size, mtl, temp, flag=None):
        self.index = index
        self.pos = pos
        self.mtl = mtl
        self.T = temp
        self.flag = flag

        # size[0] = x -- towards right
        # size[1] = y -- into secreen
        # size[2] = z -- upwards

        self.L_x = size[0]
        self.L_y = size[1]
        self.L_z = size[2]
        
        self.A_x = size[1] * size[2]
        self.A_y = size[0] * size[2]
        self.A_z = size[0] * size[1]

        self.V = size[0] * size[1] * size[2]
        self.m = self.mtl.get_density() * self.V/(1000000000) # this is now in kg
        self.m *= 1000 # now in grams

    def get_index(self):
        return self.index

    def get_flag(self):
        return self.flag

    def get_T(self):
        return self.T

    def get_A_x(self):
        return self.A_x

    def get_A_y(self):
        return self.A_y

    def get_A_z(self):
        return self.A_z

    def get_V(self):
        return self.V

    def get_m(self):
        return self.m

    def get_pos(self):
        return self.pos

    def get_heat_cpc(self):
        return self.mtl.get_specific_heat(self.T) * self.m

    def get_spec_heat(self):
        return self.mtl.get_specific_heat(self.T)

    def get_density(self):
        return self.mtl.get_density()

    def get_thermal_conductivity(self):
        return self.mtl.get_thermal_conductivity(self.T)

def create_equal_cell_mesh(n_x, n_y, n_z, L_x, L_y, L_z, mtl, T_in, T_left, T_right, T_front, T_rear):
    cells = []
    for z_i in range(n_z):
        for y_i in range(n_y):
            for x_i in range(n_x):
                
                new_cell_pos = [L_x * x_i,
                                L_y * y_i,
                                L_z * z_i]

                # boundary conditions for boundary cells
                if x_i == 0:
                    new_cell = cell([x_i, y_i, z_i], new_cell_pos, [L_x, L_y, L_z], mtl, T_left)
                elif x_i == n_x - 1:
                    new_cell = cell([x_i, y_i, z_i], new_cell_pos, [L_x, L_y, L_z], mtl, T_right)
                elif y_i == 0:
                    new_cell = cell([x_i, y_i, z_i], new_cell_pos, [L_x, L_y, L_z], mtl, T_front)
                elif y_i == n_y - 1:
                    new_cell = cell([x_i, y_i, z_i], new_cell_pos, [L_x, L_y, L_z], mtl, T_rear)

                # inner region cell
                else:
                    new_cell = cell([x_i, y_i, z_i], new_cell_pos, [L_x, L_y, L_z], mtl, T_in)
                    
                cells.append(new_cell)

    new_mesh = mesh(n_x, n_y, n_z, cells)
    return new_mesh
