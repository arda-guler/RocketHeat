import math
import matplotlib.pyplot as plt

from coolant_meshing import *
from coolant_mesh import *
from material import SS304L, CuCrZr, Jet_A1

SS304L = SS304L()
CuCrZr = CuCrZr()
Jet_A1 = Jet_A1()

def dist(pos1, pos2):
    return math.sqrt((pos2[0] - pos1[0])**2 + (pos2[1] - pos1[1])**2 + (pos2[2] - pos1[2])**2)

def conduct_differential_heat_on_mesh(mesh, T_chamber, T_coolant, dt):
    
    n_x = mesh.n_x
    n_y = mesh.n_y
    n_z = mesh.n_z

    dTs = [ [] for i in range(n_y) ]

    heat_gen_turn_positive = 0
    heat_gen_turn_negative = 0

    for i in range(len(mesh.cells)):
        for j in range(len(mesh.cells[i])):
            cell = mesh.get_cell(j, i, 0)
            cell_i_x, cell_i_y, cell_i_z = cell.get_index()[0], cell.get_index()[1], cell.get_index()[2]

            cell_left = mesh.get_cell(cell_i_x - 1, cell_i_y, cell_i_z)
            cell_right = mesh.get_cell(cell_i_x + 1, cell_i_y, cell_i_z)
            cell_front = mesh.get_cell(cell_i_x, cell_i_y - 1, cell_i_z)
            cell_rear = mesh.get_cell(cell_i_x, cell_i_y + 1, cell_i_z)

            if cell_left:
                dx = dist(cell.get_pos(), cell_left.get_pos()) * 0.001 # convert from mm to m
            elif cell_right:
                dx = dist(cell.get_pos(), cell_right.get_pos()) * 0.001 # convert from mm to m
            else:
                dx = 10**5 # if there are no left or right adjacent cells, make dx very big to prevent heat flux

            if cell_front:
                dy = dist(cell.get_pos(), cell_front.get_pos()) * 0.001 # convert from mm to m
            elif cell_rear:
                dy = dist(cell.get_pos(), cell_rear.get_pos()) * 0.001 # convert from mm to m
            else:
                dy = 10**5 # if there are no front or rear adjacent cells, make dy very big to prevent heat flux

            if cell_left:
                T_left = cell_left.get_T()
            else:
                T_left = cell.get_T()

            if cell_right:
                T_right = cell_right.get_T()
            else:
                T_right = cell.get_T()

            if cell_front:
                T_front = cell_front.get_T()
            else:
                T_front = cell.get_T()

            if cell_rear:
                T_rear = cell_rear.get_T()
            else:
                T_rear = cell.get_T()

            T_c = cell.get_T()
            
            k = cell.get_thermal_conductivity()
            rho = cell.get_density()
            spec_heat = cell.get_spec_heat()

            alpha = k/(rho * spec_heat) # thermal diffusivity

            # mm2 to m2 = 0.000001
            if cell.get_flag() == "boundary_chamber":
                heat_gen =  * cell.get_A_x() * (T_chamber - T_c)
                heat_gen_turn_positive += heat_gen
            elif cell.get_flag() == "boundary_coolant":
                heat_gen =  * cell.get_A_x() * (T_c - T_coolant)
                heat_gen_turn_negative -= heat_gen
            else:
                heat_gen = 0 # if this is not a boundary cell, no heat influx/outflux other than conduction
            
            dT = alpha * ((((T_left - 2*T_c + T_right)/(dx**2)) + ((T_front - 2*T_c + T_rear)/(dy**2))) + (heat_gen/k)) * dt
            dTs[i].append(dT)

    return dTs, heat_gen_turn_positive, heat_gen_turn_negative

def conduct_heat_in_time_interval(mesh, T_chamber, T_coolant, time_end, dt=0.001):
    time = 0
    calc_steps = int(time_end/dt)
    time_steps = 0

    heat_in = 0
    heat_out = 0
    
    while time < time_end:
        dTs, heat_flux_positive, heat_flux_negative = conduct_differential_heat_on_mesh(mesh, T_chamber, T_coolant, dt)
        for i in range(len(mesh.cells)):
            for j in range(len(mesh.cells[i])):
                cell = mesh.get_cell(j, i, 0)
                cell.T += dTs[i][j]

        heat_in += heat_flux_positive
        heat_out += heat_flux_negative

        time += dt
        time_steps +=1

        if time_steps % 100 == 0:
            print(str((time_steps/calc_steps) * 100) + "%")

    print("Heat influx:", str(heat_in))
    print("Heat outflux:", str(heat_out))

def plot_mesh_T(mesh):
    xs = []
    ys = []
    Ts = []
                
    for x_i in range(mesh.n_x):
        xs.append(mesh.get_cell(x_i, 0, 0).get_pos()[0])

    for y_i in range(mesh.n_y):
        ys.append(mesh.get_cell(0, y_i, 0).get_pos()[1])

    for y_i in range(mesh.n_y):
        Ts.append([])
        for x_i in range(mesh.n_x):
            cell = mesh.get_cell(x_i, y_i, 0)
            if cell:
                T = cell.get_T() - 273
            else:
                T = math.nan
            Ts[y_i].append(T)

    clrplot = plt.contourf(xs, ys, Ts)
    plt.colorbar()
    cplot = plt.contour(xs, ys, Ts, colors="k")
    plt.title("Temperature Gradient (degrees C)")
    plt.clabel(cplot, fontsize=10)
    plt.xlabel("X position (mm)")
    plt.ylabel("Y position (mm)")
    plt.show()

# - - - ENGINE GEOMETRY - - -
L_engine =
D_chm = 
D_thrt = 
D_exit = 

chm_contract_angle = 
ROC_chm =

nzl_expansion_angle = 
nzl_downstream_ROC = 
nzl_upstream_ROC = 

cochan_inner_wall_dist = 
cochan_radial_depth = 
cochan_angular_size = 

# - - - COMBUSTION / CEA - - -
D_star = D_thrt # m
Pc = 
r_c = nzl_downstream_ROC # m
Tc = 
c_star = 
gas_conductivity = 
avg_molecular_mass = 

Tw = 

# - - - COMBUSTION CHAMBER INPUTS - - -
visc_chm = 
gamma_chm = 

# - - - THROAT INPUTS - - -
visc_thrt =  
gamma_thrt = 

# - - - MATERIALS - - -
chm_inner_wall_material =
chm_outer_shell_material =

# - - - COOLANT - - -
coolant_fluid = 
coolant_mass_flow =
T_coolant = 

T_chamber = Tc

# L_x, L_y, cell_L_x, cell_L_y, channel_L_x, channel_L_y
mesh = create_coolant_channel_mesh()
conduct_heat_in_time_interval(mesh, T_chamber, T_coolant, )
plot_mesh_T(mesh)
