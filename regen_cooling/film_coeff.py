import math
pi = math.pi
from material import SS304L, CuCrZr, Jet_A1

SS304L = SS304L()
CuCrZr = CuCrZr()
Jet_A1 = Jet_A1()

# Bartz Equation

# D. R. Bartz - A Simple Equation for Rapid Estimation of Rocket
#               Nozzle Convective Heat Transfer Coefficients

# NASA - Experimental Studies of the Heat Transfer to RBCC Rocket
#        Nozzles for CFD Application to Design Methodologies

def get_convection_coeff(D_star, vis, Cp, Pr, Pc, c_star, r_c, A_star, A, gamma, M, Tw, T0):

    # takes:
    # D_star = throat diameter in m
    # vis = viscosity in millipoise
    # Cp = specific heat at const. press. kJ kg-1 K-1
    # Pr = prandtl number
    # Pc = stagnation press in Pa
    # r_c = throat radius of curvature in m
    # c_star = characteristic velocity in m/s
    # A_star = throat area in m2
    # A = local area in m2
    # M = local mach number

    # gives:
    # h_g = convection coefficient in W m-2 K-1

    vis = vis * 0.0001 # convert millipoise to Pa.s
    Cp = Cp * 1000 # convert kJ kg-1 K-1 to J kg-1 K-1

    # sigma = correction factor for boundary layer
    sigma = 1 / ((0.5 * (Tw/T0) * (1 + (gamma-1)/2 * M**2) + 0.5)**(0.68)) * ((1 + (gamma-1)/2 * M**2)**(0.12))

    h_g = (0.026 / D_star**0.2) * (vis**0.2 * Cp)/Pr**0.6 * (Pc/c_star)**0.8 * (D_star/r_c)**0.1 * (A_star/A)**0.9 * sigma
    
    return h_g

# Sutton Heat Transfer Analysis (8-25)
def get_coolant_film_coeff(coolant_fluid, coolant_temp, A_cochan_flow, coolant_mass_flow, D_hydro):
    
    # takes:
    # coolant fluid choice (e.g. LOX or Jet A-1)
    # coolant temperature in K
    # coolant flow area in m2
    # coolant mass flow rate in kg s-1

    # gives:
    # h_liq: convection coefficient in W m-2 K-1
    
    coolant_density = coolant_fluid.get_density(coolant_temp)
    coolant_viscosity = coolant_fluid.get_viscosity(coolant_temp)
    coolant_spec_heat = coolant_fluid.get_specific_heat(coolant_temp)
    coolant_conductivity = coolant_fluid.get_thermal_conductivity(coolant_temp)
    
    #D_eq = ((4 * A_cochan_flow) / pi)**(0.5)
    D_eq = 2 * D_hydro
    flow_vel = coolant_mass_flow / (coolant_density * A_cochan_flow)
    
    h_liq = 0.023 * coolant_spec_heat * (coolant_mass_flow/A_cochan_flow)
    h_liq *= ((D_eq * flow_vel * coolant_density)/coolant_viscosity)**(-0.2)
    h_liq *= (coolant_viscosity * coolant_spec_heat * coolant_conductivity)**(-2/3)

    return h_liq

def get_coolant_film_coeff_NEW(mtl_clt, T_clt, mdot_clt, D_hydro, cy):

    Pr = mtl_clt.get_specific_heat(T_clt) * mtl_clt.get_viscosity(T_clt) / mtl_clt.get_thermal_conductivity(T_clt)
    
    # compute Reynold's number
    # https://en.wikipedia.org/wiki/Hydraulic_diameter
    Reynolds_num = (mdot_clt * D_hydro) / (mtl_clt.get_viscosity(T_clt) * cy.A_cochan_flow)

    # compute Nusselt number
    Nusselt_num = 0.021 * (Reynolds_num**(0.8)) * (Pr**(0.4)) * (0.64 + 0.36 * (T_clt/cy.T))

    h_liq = Nusselt_num * mtl_clt.get_thermal_conductivity(T_clt) / D_hydro
    return h_liq

