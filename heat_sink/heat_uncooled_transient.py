# - - - - - - - - - - - - - - - - - - - - -
# HEAT SINK LRE THERMAL COMPUTATIONS
# - - - - - - - - - - - - - - - - - - - - -
# Program to calculate transient heat
# transfer in an uncooled (heat sink)
# LRE.
# - - - - - - - - - - - - - - - - - - - - -
# Authors:
# H. Arda. Güler
# - - - - - - - - - - - - - - - - - - - - -
# References: (1) Rocket Propulsion Elements 7th Ed. by George P. Sutton
# - - - - - - - - - - - - - - - - - - - - -

from datetime import datetime
import matplotlib.pyplot as plt
import os
import math

from material import SS304L, CuCrZr
from film_coeff import *

pi = math.pi
uni_gas_const = 8.314472 # m2 kg s-2 K-1 mol-1

class segment:
    def __init__(self, index, r_in, thickness, height, temp, material):
        self.index = index                                              # unitless
        self.r_in = r_in                                                # m
        self.thickness = thickness                                      # m
        self.height = height                                            # m
        self.temp = temp                                                # K
        self.material = material
        
        self.r_out = self.r_in + self.thickness                         # m
        self.in_area = 2 * pi * self.r_in * self.height                 # m2
        self.out_area = 2 * pi * self.r_out * self.height               # m2
        self.volume = pi * self.height * (self.r_out**2 - self.r_in**2) # m3

    def get_material(self):
        # returns material object
        return self.material

    def get_temp(self):
        # returns average temp of cell in K
        return self.temp

    def get_mass(self):
        # returns mass of cell in kg
        return self.get_material().get_density() * self.volume

    def get_inner_radius(self):
        # return inner radius in m
        return self.r_in

    def get_outer_radius(self):
        # return outer radius in m
        return self.r_out

    def get_mid_radius(self):
        # return mid-point radius in m
        return (self.r_in + self.r_out) * 0.5

    def get_inner_area(self):
        # return inner facing area in m2
        return self.in_area

    def get_outer_area(self):
        # return outer area in m2
        return self.out_area

    def get_thickness(self):
        # return thickness in m
        return self.thickness

# init materials
##SS = SS304L()
##CCZ = CuCrZr()

## - - - - - - - - - - - - - - - - - - - - -
##                  INPUTS
## - - - - - - - - - - - - - - - - - - - - -

# - - - COMMON INPUTS - - -
total_thickness = # wall thickness, m
fineness = # number of wall slices                                  
segment_height = # m
segment_init_temp = + 273 # K                              
segment_material = 

D_star =  # m
Pc =  # Pa
r_c =  # m
chamber_temp =  # K
c_star =  # m/s, CEA
gas_conductivity = # W m-1 K-1, CEA
avg_molecular_mass = # g mol-1



# - - - COMBUSTION CHAMBER INPUTS - - -
chamber_inner_radius = # m
Mach_chamber = 0 # assume no velocity in combustion chamber
viscosity_chamber = # millipoise, CEA
gamma_chamber = # CEA
Cp_chamber = (gamma_chamber/(gamma_chamber-1)) * uni_gas_const / avg_molecular_mass # kJ kg-1 K-1, CEA
Pr_chamber = (4*gamma_chamber) / (9*gamma_chamber - 5) # unitless



# - - - THROAT INPUTS - - -
throat_radius = # m
Mach_throat = 1 # obvious
viscosity_throat = # millipoise, CEA
gamma_throat = # CEA
Cp_throat = (gamma_throat/(gamma_throat-1)) * uni_gas_const / avg_molecular_mass # kJ kg-1 K-1, CEA
Pr_throat = (4*gamma_throat) / (9*gamma_throat - 5) # unitless



# - - - GEOMETRY CALCULATIONS - - -
A_star = pi * throat_radius**2 # m2
A_chamber = pi * chamber_inner_radius**2 # m2
A_throat = A_star

segment_thickness = total_thickness/fineness



# - - - ANALYSIS SETUP - - -
time_incr = # s
engine_shutdown_time = # s
time_end = # s



## - - - - - - - - - - - - - - - - - - - - -
##              CALCULATIONS
## - - - - - - - - - - - - - - - - - - - - -

# Rocket Propulsion Elements, Eqn. (8-24) at p. 315
# Transient Heat Transfer Analysis of an Uncooled Metal Combustion Chamber

def perform_heat_analysis(inner_radius, seg_thickness, seg_height, seg_T0, seg_material,\
                          D_star, vis, Cp, Pr, Pc, c_star, r_c, A_star, A, gamma, M, chamber_temp,\
                          time_incr, time_end, engine_shutdown_time):

    temps_data = []
    h_list = []
    segments = []
    time = 0
    time_run = 0
    total_heat = 0
    
    # generate wall segments
    for i in range(fineness):
        new_r_in = inner_radius + i * seg_thickness
        new_segment = segment(i, new_r_in, seg_thickness, seg_height, seg_T0, seg_material)
        segments.append(new_segment)
    
    while time <= time_end:
        temps_data.append([time, [], []])
        
        for segment_index in range(len(segments)):
            # get properties of current segment
            current_segment = segments[segment_index]
            s_mass = current_segment.get_mass()
            s_inner_area = current_segment.get_inner_area()
            s_outer_area = current_segment.get_outer_area()
            s_thickness = current_segment.get_thickness()
            s_temp = current_segment.get_temp()
            s_spec_heat = current_segment.get_material().get_specific_heat(s_temp)
            s_thermal_conductivity = current_segment.get_material().get_thermal_conductivity(s_temp)

            if segment_index == 0 and time < engine_shutdown_time:
                # add heat to the innermost segment from combustion
                gas_convection_coeff = get_convection_coeff(D_star, vis, Cp, Pr, Pc, c_star, r_c, A_star, A, gamma, M, s_temp, chamber_temp)
                Q_transfer = gas_convection_coeff * (chamber_temp - s_temp) * s_inner_area * time_incr

                h_list.append([gas_convection_coeff, time])
                total_heat += Q_transfer

            else:
                # if this is not the innermost segment, use Q_transfer from inner neighboring layer
                pass
            
            current_segment.temp += Q_transfer / (s_mass * s_spec_heat)

            if not segment_index == len(segments) - 1:
                # this is not the outermost segment, transfer some heat to the next layer
                temp_diff = s_temp - segments[segment_index + 1].get_temp()
                Q_transfer = s_thermal_conductivity * s_outer_area * (temp_diff/s_thickness) * time_incr
                
                # remove Q_transfer amount of heat from current layer, reducing temperature
                current_segment.temp -= Q_transfer/(s_mass * s_spec_heat)

            else:
                # this is the outermost segment, nowhere to transfer heat (uncooled chamber)
                pass

            temps_data[time_run][1].append(current_segment.get_mid_radius() * 1000) # convert to mm for plot
            temps_data[time_run][2].append(current_segment.get_temp() - 273) # convert to degrees C for plot

        time += time_incr
        time_run += 1

    return temps_data, h_list, total_heat

# provide a readback of the input so user doesn't miss anything
print("TRANSIENT HEAT ANALYSIS - UNCOOLED\n")
print("Input Parameters:")
print("Chamber inner radius:", chamber_inner_radius, "m")
print("Segment thickness:", segment_thickness, "m")
print("Segment height:", segment_height, "m")
print("Wall initial temperature:", segment_init_temp, "K")
print("Thrust chamber material:", segment_material.get_name())
print("Pc:", Pc, "Pa")
print("C*:", c_star, "m/s")
print("D*:", D_star, "m")
print("Throat radius of curvature:", r_c, "m")
print("A*:", A_star, "m^2")
print("Tc:", chamber_temp, "K\n")

print("Viscosity (chamber):", viscosity_chamber, "millipoise")
print("Cp (chamber):", Cp_chamber, "kJ kg-1 K-1")
print("Pr (chamber):", Pr_chamber)
print("A (chamber):", A_chamber, "m^2")
print("Gamma (chamber):", gamma_chamber)
print("Mach num (chamber):", Mach_chamber, "\n")

print("Viscosity (throat):", viscosity_throat, "millipoise")
print("Cp (throat):", Cp_throat, "kJ kg-1 K-1")
print("Pr (throat):", Pr_throat)
print("A (throat):", A_throat, "m^2")
print("Gamma (throat):", gamma_throat)
print("Mach num (throat):", Mach_throat, "\n")

print("Engine shutdown time:", engine_shutdown_time, "s")
print("Time increment:", time_incr, "s")
print("Time end:", time_end, "s\n")
# --- readback ends ---

# Instruct script to perform two analyses, one with chamber input and another with throat
print("Running chamber analysis...")
chamber_temps_list, chamber_h_list, chamber_total_heat = perform_heat_analysis(chamber_inner_radius, segment_thickness, segment_height, segment_init_temp,\
                                                                               segment_material, D_star, viscosity_chamber, Cp_chamber, Pr_chamber, Pc, c_star,\
                                                                               r_c, A_star, A_chamber, gamma_chamber, Mach_chamber, chamber_temp,\
                                                                               time_incr, time_end, engine_shutdown_time)

print("Running throat analysis...")
throat_temps_list, throat_h_list, throat_total_heat = perform_heat_analysis(throat_radius, segment_thickness, segment_height, segment_init_temp,\
                                                                            segment_material, D_star, viscosity_throat, Cp_throat, Pr_throat, Pc, c_star,\
                                                                            r_c, A_star, A_throat, gamma_throat, Mach_throat, chamber_temp,\
                                                                            time_incr, time_end, engine_shutdown_time)

## - - - - - - - - - - - - - - - - - - - - -
##              PLOTTING
## - - - - - - - - - - - - - - - - - - - - -

def plot_heat_analysis(temps_data, h_list, display=[1,True,""], style="2d"):

    plt.figure(display[0])
    if style == "3d":
        ax = plt.axes(projection = "3d")
        for i in range(0, len(temps_data), int(len(temps_data)/8)):
            ax.plot3D(temps_data[i][1], temps_data[i][2], temps_data[i][0], 'green')
                
        plt.title("Temp. Gradient Across Wall (Transient)")
        plt.xlabel("Wall Position (mm)")
        plt.ylabel("Temperature (degrees C)")

    else:
        for i in range(0, len(temps_data), int(len(temps_data)/12)):
            red = min(1, max(temps_data[i][2])/segment_material.get_melting_point("C"))
            blue = 1 - red
            plt.plot(temps_data[i][1], temps_data[i][2], color=(red, 0, blue))
            plt.plot(temps_data[i][1], temps_data[i][2], "ro", color=(red, 0, blue), markersize=3)

        plt.grid()
        plt.title("%s\nTemperature Gradient Across %s Wall (Transient)\nuntil T = %.2f s" % (display[2], segment_material.get_name(), time_end))
        plt.xlabel("Wall Position (mm)")
        plt.ylabel("Temperature (degrees C)")

    fig_a = plt.gcf()

    plt.figure(display[0]+1)
    plt.title("%s\nConvective Heat Transfer Coefficient vs. Time" % (display[2]))
    for i in range(0, len(h_list), int(len(h_list)/25)):
        plt.plot(h_list[i][1], h_list[i][0], color=(1, 0.8, 0))
        plt.plot(h_list[i][1], h_list[i][0], "ro", markersize=3)

    plt.grid()
    plt.xlabel("Time (s)")
    plt.ylabel("Convective Heat Transfer Coeff. (W m-2 K-1)")

    fig_b = plt.gcf()

    if display[1]:
        plt.show()

    return fig_a, fig_b

print("\nOutput:")
print("Chamber segment heat (J):", chamber_total_heat)
print("Throat segment heat (J):", throat_total_heat, "\n")

# save the plots because user might want to export them later
fig_chamber, fig_chamber_h = plot_heat_analysis(chamber_temps_list, chamber_h_list, [1, True, "Combustion Chamber"])
fig_throat, fig_throat_h = plot_heat_analysis(throat_temps_list, throat_h_list, [3, True, "Throat"])

## - - - - - - - - - - - - - - - - - - - - -
##              EXPORTING
## - - - - - - - - - - - - - - - - - - - - -

# does user want to export this analysis?
export_go = input("EXPORT HEAT ANALYSIS DATA? (y/N):")

# if answer isn't yes, we can stop the program now
if not export_go.lower() == "y":
    quit()

# user said yes, export the data
# name the analysis with current time
now = datetime.now()
date_time = now.strftime("%Y-%m-%d_%H-%M-%S")
folder_name = "heat_analysis_" + date_time

# create new folder to contain analysis files
print("Creating folder...")
if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# write parameters and text output to file
print("Writing data to analysis.txt...")
with open(folder_name + "/analysis.txt", "w") as f:
    f.write("Chamber inner radius: " + str(chamber_inner_radius) + " m\n")
    f.write("Segment thickness: " + str(segment_thickness) + " m\n")
    f.write("Segment height: " + str(segment_height) + " m\n")
    f.write("Wall initial temperature: " + str(segment_init_temp) + " K\n")
    f.write("Thrust chamber material: " + str(segment_material.get_name()) + "\n")
    f.write("Pc: " + str(Pc) + " Pa\n")
    f.write("C*: " + str(c_star) + " m/s\n")
    f.write("D*: " + str(D_star) + " m\n")
    f.write("Throat radius of curvature: " + str(r_c) + " m\n")
    f.write("A*: " + str(A_star) + " m^2\n")
    f.write("Tc: " + str(chamber_temp) + " K\n\n")

    f.write("Viscosity (chamber): " + str(viscosity_chamber) + " millipoise\n")
    f.write("Cp (chamber): " + str(Cp_chamber) + " kJ kg^-1 K^-1\n")
    f.write("Pr (chamber): " + str(Pr_chamber) + "\n")
    f.write("A (chamber): " + str(A_chamber) + " m^2\n")
    f.write("Gamma (chamber): " + str(gamma_chamber) + "\n")
    f.write("Mach num (chamber): " + str(Mach_chamber) + "\n\n")

    f.write("Viscosity (throat): " + str(viscosity_throat) + " millipoise\n")
    f.write("Cp (throat): " + str(Cp_throat) + " kJ kg-1 K-1\n")
    f.write("Pr (throat): " + str(Pr_throat) + "\n")
    f.write("A (throat): " + str(A_throat) + " m^2\n")
    f.write("Gamma (throat): " + str(gamma_throat) + "\n")
    f.write("Mach num (throat): " + str(Mach_throat) + "\n\n")

    f.write("Engine shutdown time: " + str(engine_shutdown_time) + " s\n")
    f.write("Time increment: " + str(time_incr) + " s\n")
    f.write("Time end: " + str(time_end) + "s\n\n")

    f.write("Heat transfer to analyzed chamber segment for full duration: " + str(chamber_total_heat) + " J\n")
    f.write("Heat transfer to analyzed throat segment for full duration: " + str(throat_total_heat) + " J\n")

# export plots
print("Exporting plots...")

fig_chamber.savefig(folder_name + '/chamber_temp.png')
fig_chamber_h.savefig(folder_name + '/chamber_coeff.png')
fig_throat.savefig(folder_name + '/throat_temp.png')
fig_throat_h.savefig(folder_name + '/throat_coeff.png')
