import numpy as np
from eMach.mach_eval.analyzers.mechanical.rotor_thermal import SPM_RotorThermalProblem,SPM_RotorThermalAnalyzer
# Example Machine Dimensions
r_sh=5E-3 # [m]
d_m=3E-3 # [m]
r_ro=12.5E-3 # [m]
d_ri=r_ro-r_sh - d_m # [m]
d_sl=1E-3 # [m]
l_st=50E-3 # [m]
l_hub=3E-3 # [m]
r_si=r_ro+d_sl+1E-3 # [m]

# Define Material Dictionary
mat_dict= {'shaft_therm_conductivity': 51.9, # W/m-k ,
           'core_therm_conductivity': 28, # W/m-k
           'magnet_therm_conductivity': 8.95, # W/m-k ,
           'sleeve_therm_conductivity': 0.71, # W/m-k,
           'air_therm_conductivity'     :.02624, #W/m-K
           'air_viscosity'              :1.562E-5, #m^2/s
           'air_cp'                     :1, #kJ/kg
           'rotor_hub_therm_conductivity':205.0} #W/m-K}
# Operating Conditions
T_ref=25 # [C]
omega=120E3*2*np.pi/60 # [rad/s]
losses={'rotor_iron_loss':.001,'magnet_loss':135}
u_z=0

prob=SPM_RotorThermalProblem(mat_dict,r_sh,d_ri,r_ro,d_sl,r_si,l_st,l_hub,T_ref,u_z,losses,omega)
ana=SPM_RotorThermalAnalyzer()

T=ana.analyze(prob)