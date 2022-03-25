import time
import datetime
import os
import matplotlib.pyplot as plt
import shutil

def plot_data(time_step, xs, cylinder_temps, coolant_temps, Q_ins, Q_outs, Reynolds, Nusselts, T_gases,
              h_gs, h_ls, clt_vels, Q_in_fulls, Q_out_fulls, geom_x, geom_y,
              flow_areas, wet_perimeters, D_hydros, filename=None):

    # PRINT TOTAL Q
    Q_in_total = 0
    for Q in Q_in_fulls:
        Q_in_total += Q

    Q_out_total = 0
    for Q in Q_out_fulls:
        Q_out_total += Q

    print("\nTotal Q_in:", Q_in_total * 100)
    print("Total Q_out:", Q_out_total * 100)
    print("Net Q:", (Q_in_total - Q_out_total) * 100)
    
    num_frames = len(cylinder_temps)
    fig, ax = plt.subplots()

    # WALL TEMP. PLOT
    plt.figure(1)

    for i in range(0, num_frames, int(num_frames/10)):
        red = min(1, max(cylinder_temps[i])/600)
        blue = 1 - red
        plt.plot(xs, cylinder_temps[i], color=(red, 0, blue))

    plt.grid()
    plt.title("Wall Temperature")
    plt.xlabel("Position (m)")
    plt.ylabel("Temperature (C)")

    # COOLANT TEMP. PLOT
    plt.figure(2)

    for i in range(0, num_frames, int(num_frames/10)):
        red = min(1, max(coolant_temps[i])/350)
        blue = 1 - red
        plt.plot(xs, coolant_temps[i], color=(red, 0, blue))

    plt.grid()
    plt.title("Coolant Temperature")
    plt.xlabel("Position (m)")
    plt.ylabel("Temperature (C)")

    # HEAT PLOT
    plt.figure(3)

    for i in range(0, num_frames, int(num_frames/10)):
        plt.plot(xs, Q_ins[i], color=(1,0,0))
        plt.plot(xs, Q_outs[i], color=(0,0,1))

    plt.grid()
    plt.title("Heat Transfer")
    plt.xlabel("Position (m)")
    plt.ylabel("Heat Transferred (W)")

    # REYNOLDS NUMBER PLOT
    plt.figure(4)

    for i in range(0, num_frames, int(num_frames/10)):
        plt.plot(xs, Reynolds[i])

    plt.grid()
    plt.title("Reynolds Number")
    plt.xlabel("Position (m)")
    plt.ylabel("Reynolds Number")

    # NUSSELT NUMBER PLOT
    plt.figure(5)

    for i in range(0, num_frames, int(num_frames/10)):
        plt.plot(xs, Nusselts[i])

    plt.grid()
    plt.title("Nusselts Number")
    plt.xlabel("Position (m)")
    plt.ylabel("Nusselts Number")

    # COMBUSTION GAS TEMP. PLOT
    plt.figure(6)

    plt.plot(xs, T_gases)

    plt.grid()
    plt.title("Gas Temperature (K)")
    plt.xlabel("Position (m)")
    plt.ylabel("Gas Temperature")

    # GAS CONVECTION COEFF. PLOT
    plt.figure(7)

    for i in range(0, num_frames, int(num_frames/10)):
        plt.plot(xs, h_gs[i])

    plt.grid()
    plt.title("hg")
    plt.xlabel("Position (m)")
    plt.ylabel("hg")

    # LIQUID FILM COEFF. PLOT
    plt.figure(8)

    for i in range(0, num_frames, int(num_frames/10)):
        plt.plot(xs, h_ls[i])

    plt.grid()
    plt.title("hl")
    plt.xlabel("Position (m)")
    plt.ylabel("hl")

    # COOLANT VELOCITY PLOT
    plt.figure(9)

    for i in range(0, num_frames, int(num_frames/10)):
        plt.plot(xs, clt_vels[i])

    plt.grid()
    plt.title("Coolant Flow Velocity")
    plt.xlabel("Position (m)")
    plt.ylabel("Velocity (m s-1)")

    # TOTAL Q IN PLOT
    plt.figure(10)

    plt.plot(Q_in_fulls)

    plt.grid()
    plt.title("Total Q In")
    plt.xlabel("Time")
    plt.ylabel("Total Q In")

    # TOTAL Q OUT PLOT
    plt.figure(11)

    plt.plot(Q_out_fulls)

    plt.grid()
    plt.title("Total Q Out")
    plt.xlabel("Time")
    plt.ylabel("Total Q Out")

##    # HEAT PLOT (3D)
##    plt.figure(12)
##    ax = plt.axes(projection='3d')
##    
##    times = []
##    for i in range(num_frames):
##        times.append(i*time_step)
##
##    for i in range(0, num_frames, int(num_frames/10)):
##        ax.plot3D(Q_ins[i], xs, times[i], color="red")
##
##    for i in range(0, num_frames, int(num_frames/10)):
##        ax.plot3D(Q_outs[i], xs, times[i], color="blue")
##
##    # HEAT DIFF. PLOT (3D)
##    plt.figure(13)
##    ax = plt.axes(projection='3d')
##
##    Q_nets = []
##    for i in range(len(Q_ins)):
##        Q_nets.append([])
##        for j in range(len(Q_ins[0])):
##            Q_nets[i].append(Q_ins[i][j] - Q_outs[i][j])
##    
##    for i in range(0, num_frames, int(num_frames/10)):
##        ax.plot3D(Q_nets[i], xs, times[i], color="green")

    # ENGINE GEOMETRY
    plt.figure(12)

    plt.axes().set_aspect('equal')
    plt.plot(geom_x, geom_y)
    geom_y_negative = []
    for y in geom_y:
        geom_y_negative.append(-y)
    plt.plot(geom_x, geom_y_negative)

    plt.grid()
    plt.title("Thrust Chamber Geometry")
    plt.xlabel("X")
    plt.ylabel("Y")

    # FLOW AREA
    plt.figure(13)
    plt.plot(xs, flow_areas)
    plt.grid()
    plt.title("Coolant Flow Area")
    plt.xlabel("X")
    plt.ylabel("Area m2")

    # WET PERIMETER
    plt.figure(14)
    plt.plot(xs, wet_perimeters)
    plt.grid()
    plt.title("Wet Perimeter")
    plt.xlabel("X")
    plt.ylabel("Perimeter m")

    # HYDRAULIC DIAMETERS
    plt.figure(15)
    plt.plot(xs, D_hydros)
    plt.grid()
    plt.title("Hydraulic Diameter")
    plt.xlabel("X")
    plt.ylabel("Diameter m")

    if not filename:
        folder_name = "heat_analysis_" + datetime.datetime.now().strftime("%y%m%d%H%M%S")
    else:
        folder_name = filename.split("/")[2].split(".")[0]
        
    print("Exporting figures to folder: " + folder_name)

    try:
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)
    except:
        print("ERROR: Could not create folder. Try saving figures by hand.")
        plt.show()
        return

    try:
        shutil.copy(filename, folder_name)
    except:
        print("WARNING: Could not copy inputs file to analysis folder.")

    try:
        for i in range(1, 16):
            new_fig = plt.figure(i)
            save_str = folder_name + "/figure_" + str(i) + ".png"
            new_fig.savefig(save_str)
    except:
        print("ERROR: Could not save some or all of the figures. Try saving figures by hand.")
        plt.show()
        return

    print("Figures exported successfully!")
    
    print("Clearing figures from memory...")
    for i in range(1, 16):
        new_fig = plt.figure(i)
        new_fig.clear()
        plt.close()

    print("Analysis done!")
    # qc = input("Press Enter to quit...")
