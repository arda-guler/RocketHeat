import matplotlib.pyplot as plt
from material import *

mtl = CuCrZr()

T_wall =  + 273
T_air =  + 273

print("Initial T:", T_wall, "K")

# r2 = inner wall radius + channel size + outer shell
r2 = 
r1 = 
area_factor = r2/r1

D = 2 * r2
A =  * area_factor
print("Area:", A, "m2")

# TO-DO: Learn to actually calculate this
h = 7 # W m-2 K-1

mass = 

t_step = 0.001
t_end = 3600
t = 0
Ts = []
ts = []

while t < t_end:
    spec_heat = mtl.get_specific_heat(T_wall)

    heat_capacity = mass * spec_heat

    Q = (T_wall - T_air) * A * h * t_step
    T_wall -= Q/heat_capacity
    Ts.append(T_wall - 273)
    ts.append(t/60)
    t += t_step

print("Last T:", T_wall, "K")
plt.plot(ts, Ts)
plt.title("Natural Convection Air Cooling")
plt.ylabel("Temperature (C)")
plt.xlabel("Time after shutdown (min)")
plt.show()
