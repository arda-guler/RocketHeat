import time
import matplotlib.pyplot as plt

def plot_data(time_step, xs, cylinder_temps, coolant_temps, Q_ins, Q_outs, Reynolds, Nusselts, T_gases,
              h_gs, h_ls, clt_vels, Q_in_fulls, Q_out_fulls):

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

    # HEAT PLOT (3D)
    plt.figure(12)
    ax = plt.axes(projection='3d')
    
    times = []
    for i in range(num_frames):
        times.append(i*time_step)

    for i in range(0, num_frames, int(num_frames/10)):
        ax.plot3D(Q_ins[i], xs, times[i], color="red")

    for i in range(0, num_frames, int(num_frames/10)):
        ax.plot3D(Q_outs[i], xs, times[i], color="blue")

    # HEAT DIFF. PLOT (3D)
    plt.figure(13)
    ax = plt.axes(projection='3d')

    Q_nets = []
    for i in range(len(Q_ins)):
        Q_nets.append([])
        for j in range(len(Q_ins[0])):
            Q_nets[i].append(Q_ins[i][j] - Q_outs[i][j])
    
    for i in range(0, num_frames, int(num_frames/10)):
        ax.plot3D(Q_nets[i], xs, times[i], color="green")

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
    plt.title("Combustion Gas Convection")
    plt.xlabel("Position (m)")
    plt.ylabel("hg (W m-2 K-1)")

    # LIQUID FILM COEFF. PLOT
    plt.figure(8)

    for i in range(0, num_frames, int(num_frames/10)):
        plt.plot(xs, h_ls[i])

    plt.grid()
    plt.title("Coolant Convection")
    plt.xlabel("Position (m)")
    plt.ylabel("hl (W m-2 K-1)")

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

    plt.show()
