# - - - - - - - - - - - - - - - - - - - -
# REGENERATIVE HEAT TRANSFER COMPUTATIONS
# - - - - - - - - - - - - - - - - - - - -
# Program to compute transient heat
# transfers in an LRE.
# - - - - - - - - - - - - - - - - - - - -
# Authors:
# H. Arda Güler
# - - - - - - - - - - - - - - - - - - - -

import math
pi = math.pi
import re

from film_coeff import *
from geometry import *
from material import SS304L, CuCrZr, Jet_A1
from mach import *
from plot import *

SS = SS304L()
CCZ = CuCrZr()
JetA1 = Jet_A1()

materials = [SS, CCZ, JetA1]

def get_material_by_name(line):
    global materials
    
    if "SS" in line:
        return materials[0]
    elif "CCZ" in line:
        return materials[1]
    elif "Jet A1" in line:
        return materials[2]

config_filename = input("Enter filename to import values from: ")
if not config_filename.endswith(".txt"):
    config_filename += ".txt"
config_file = open(config_filename, "r")
config_lines = config_file.readlines()

# - - - ENGINE GEOMETRY - - -
L_engine = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[0])[0]) # m
D_chm = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[1])[0]) # m
D_thrt = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[2])[0]) # m
D_exit = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[3])[0]) # m
A_star = pi * (D_thrt/2)**2 # m2

a_chmContract = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[4])[0]) # deg
ROC_chm = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[5])[0]) # m

a_nzlExp = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[6])[0]) # deg
ROC_thrtDn = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[7])[0]) # m
ROC_thrtUp = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[8])[0]) # m

n_cochan = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[9])[0]) # number of coolant channels
L_cochanInnerWallDist = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[10])[0]) # m
L_cochanSideWall = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[11])[0]) # m
L_cochanDepth = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[12])[0]) # m

# - - - COMBUSTION / CEA - - -
D_star = D_thrt # m
P_c = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[13])[0]) # Pa
r_c = ROC_thrtDn # m
T_c = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[14])[0]) # K
c_star = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[15])[0]) # m/s, CEA
gasConductivity = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[16])[0]) # W m-1 K-1, CEA
avgMolecularMass = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[17])[0]) # g mol-1

T_w = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[18])[0]) # K, wall temp

# - - - COMBUSTION CHAMBER INPUTS - - -
visc_chm = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[19])[0]) # millipoise, CEA
gamma_chm = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[20])[0]) # CEA

# - - - THROAT INPUTS - - -
visc_thrt =  float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[21])[0]) # millipoise, CEA
gamma_thrt = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[22])[0]) # CEA

# - - - MATERIALS - - -
mtl_innerWall = get_material_by_name(config_lines[23])
mtl_outerShell = get_material_by_name(config_lines[24])

# - - - COOLANT - - -
mtl_clt = get_material_by_name(config_lines[25])
mdot_clt = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[26])[0])/n_cochan # kg s-1 (per channel)
T_clt = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[27])[0]) # manifold coolant temp. K

# - - - ANALYSIS - - -
fineness_vertical = int(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[28])[0])
time_end = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[29])[0]) # s
time_step = float(re.findall(r"(?<![a-zA-Z:])[-+]?\d*\.?\d+", config_lines[30])[0]) # s
n_steps = int(time_end/time_step)

# calculate engine geometry
geom_x, geom_y, x_step, engine_lengths = calculate_geometry(L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm,
                                                            a_nzlExp, ROC_thrtDn, ROC_thrtUp, fineness_vertical)

# calculate Mach distribution
subsonic_x, subsonic_M, supersonic_x, supersonic_M = calc_mach_num(L_engine, engine_lengths[4], T_c, gamma_thrt, avgMolecularMass, fineness_vertical,
                                                                   L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)

# generate cylinders
m_engine = 0
cylinders = []
for i in range(fineness_vertical):
    x = i * x_step
    r_in = get_inner_radius_at(x, L_engine, D_chm, D_thrt, D_exit, a_chmContract, ROC_chm, a_nzlExp, ROC_thrtDn, ROC_thrtUp)
    r_clt = r_in + L_cochanInnerWallDist
    r_out = r_clt + L_cochanDepth
    a_clt = (360/n_cochan) - (360 * L_cochanSideWall)/(2 * pi * r_clt)
    Mach = get_mach_num_at(x, subsonic_M, subsonic_x, supersonic_M, supersonic_x, engine_lengths)
    new_cylinder = cylinder(x, r_in, r_out, x_step, n_cochan, r_clt, a_clt, mtl_innerWall, T_w, Mach)
    cylinders.append(new_cylinder)
    m_engine += new_cylinder.get_m()

# calculate Cp and Pr
Cp_chm = (gamma_chm/(gamma_chm-1)) * uni_gas_const / avgMolecularMass # kJ kg-1 K-1, CEA
Pr_chm = (4*gamma_chm) / (9*gamma_chm - 5) # unitless

Cp_thrt = (gamma_thrt/(gamma_thrt-1)) * uni_gas_const / avgMolecularMass # kJ kg-1 K-1, CEA
Pr_thrt = (4*gamma_thrt) / (9*gamma_thrt - 5) # unitless

Q_ins = []
Q_outs = []
xs = []
cylinder_temps = []
coolant_temps = []
Reynolds = []
Nusselts = []
T_gases = []
h_gs = []
h_ls = []
clt_vels = []
Q_in_fulls = []
Q_out_fulls = []
flow_areas = []
wet_perimeters = []
D_hydros = []
    
time = 0
j = 0
for t_step in range(n_steps):

    if t_step % 100 == 0:
        Q_ins.append([])
        Q_outs.append([])
        cylinder_temps.append([])
        coolant_temps.append([])
        Reynolds.append([])
        Nusselts.append([])
        h_gs.append([])
        h_ls.append([])
        clt_vels.append([])
        Q_in_full = 0
        Q_out_full = 0

    T_clt_current = T_clt # revert to manifold temperature
    
    # loop backwards when computing cylinder heat transfer (go from manifold to injector face)
    for cylinder in cylinders[::-1]:
        
        if cylinder.x < engine_lengths[4]:
            vis = visc_chm
            gamma = gamma_chm
            Cp = Cp_chm
            Pr = Pr_chm
        else:
            vis = visc_thrt
            gamma = gamma_thrt
            Cp = Cp_thrt
            Pr = Pr_thrt

        M = cylinder.get_Mach()
        A = pi * cylinder.r_in**2

        # calculate flow temperature at point
        T_gas = (1 + (gamma-1)/2 * M**2)**(-1) * T_c
            
        # calculate heat transfer
        h_g = get_convection_coeff(D_star, vis, Cp, Pr, P_c, c_star, r_c, A_star, A, gamma, M, cylinder.T, T_c)
        Q_in = h_g * (T_gas - cylinder.T) * cylinder.get_A_chm() * time_step
        Q_in_full += Q_in


        # compute Reynold's number
        # https://en.wikipedia.org/wiki/Hydraulic_diameter
        wet_perimeter = (2 * pi * cylinder.r_clt) * (cylinder.a_clt/360) + (2 * pi * cylinder.r_out) * (cylinder.a_clt/360) + 2 * L_cochanDepth
        flow_area = cylinder.A_cochan_flow
        D_hydro = 4 * (flow_area / wet_perimeter)
        Reynolds_num = (mdot_clt * D_hydro) / (mtl_clt.get_viscosity(T_clt_current) * cylinder.A_cochan_flow)
        
        h_l = get_coolant_film_coeff(mtl_clt, T_clt_current, cylinder.get_A_cochan_flow(), mdot_clt, D_hydro)
        Q_out = h_l * (cylinder.T - T_clt_current) * cylinder.get_A_clt() * time_step
        Q_out_full += Q_out

        # increase cylinder temp
        Q_net = Q_in - Q_out
        dT = Q_net/cylinder.get_heat_capacity()
        cylinder.T += dT

        # increase coolant fluid temp.
        m_flow = mdot_clt * time_step
        dT_clt = (Q_out/n_cochan)/(m_flow * mtl_clt.get_specific_heat(T_clt_current))
        T_clt_current += dT_clt

        clt_vel = mdot_clt / (mtl_clt.get_density(T_clt_current) * cylinder.A_cochan_flow)

        # compute Nusselt number (Dittus Boelter)
        Nusselt_num = 0.023 * Reynolds_num**0.8 * Pr**0.3

        # record data for plotting
        if t_step == 0:
            xs.insert(0, cylinder.x)
            T_gases.insert(0, T_gas)
            flow_areas.insert(0, flow_area)
            wet_perimeters.insert(0, wet_perimeter)
            D_hydros.insert(0, D_hydro)
        if t_step % 100 == 0:
            Q_ins[j].insert(0, Q_in/time_step) # convert to W
            Q_outs[j].insert(0, Q_out/time_step) # convert to W
            cylinder_temps[j].insert(0, cylinder.T - 273) # convert to celcius
            coolant_temps[j].insert(0, T_clt_current - 273) # convert to celcius
            Reynolds[j].insert(0, Reynolds_num)
            Nusselts[j].insert(0, Nusselt_num)
            h_gs[j].insert(0, h_g)
            h_ls[j].insert(0, h_l)
            clt_vels[j].insert(0, clt_vel)

    if t_step % 100 == 0:
        Q_in_fulls.append(Q_in_full)
        Q_out_fulls.append(Q_out_full)

    # proceed to next time step
    time += time_step
    if t_step % 100 == 0:
        j+=1
        print(str((t_step/n_steps) * 100), "%")

plot_data(time_step, xs, cylinder_temps, coolant_temps, Q_ins, Q_outs, Reynolds, Nusselts,
          T_gases, h_gs, h_ls, clt_vels, Q_in_fulls, Q_out_fulls, geom_x, geom_y,
          flow_areas, wet_perimeters, D_hydros)
