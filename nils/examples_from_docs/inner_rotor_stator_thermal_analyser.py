#%% 5.2 Inputs from User

import sys
import numpy as np
import eMach.mach_eval.analyzers.mechanical.thermal_stator as sta
from matplotlib import pyplot as plt

Q= 6 #Number of Slots
g_sy = 10000 #Volumetric losses in yoke [W/m^3]
g_th = 100000 #Volumetric losses in tooth [W/m^3]
w_tooth = 0.02 #Tooth width [m]
l_st = 0.05 #Stack length
alpha_q = np.pi/Q #slot span [rad]
r_si =0.03 #Inner stator radius [m]
r_so = 0.1 #Outer stator radius [m]
r_sy = .08 #Inner stator yoke radius [m]
# r_sy = .04 #Inner stator yoke radius [m]
k_ins = 1 #Insulation thermal conductivity [W/m-K]
w_ins =.001 #Insulation Thickness [m]
k_fe = 38 #Stator iron thermal conductivity [W/m-k]
h = 100 #Exterior convection coefficient [W/m^2-k]
alpha_slot = .7 *alpha_q # back of slot span [rad]
Q_coil = 40 # Coil losses [W]
h_slot =0 #Inslot convection coefficient [W/m^2-K]
T_ref = 20 #Reference Temperature

problem = sta.StatorThermalProblem(
            g_sy,
            g_th,
            w_tooth,
            l_st,
            alpha_q,
            r_si,
            r_so,
            r_sy,
            k_ins,
            w_ins,
            k_fe,
            h,
            alpha_slot,
            Q_coil,
            h_slot,
            T_ref
        )
ana = sta.StatorThermalAnalyzer()

#%% 5.3 Outputs to User

results = ana.analyze(problem)
print(results)

sys.exit()

{'Coil temperature': 216.31038291260649, 'Stator yoke temperature': 204.06667848224436}

l_tooth_vect=np.linspace(.01,.8,100)
T_coil_vect=np.zeros_like(l_tooth_vect)
for ind,l_tooth in enumerate(l_tooth_vect):
    r_sy= l_tooth+r_si
    r_so= r_sy+.2
    problem = sta.StatorThermalProblem(
            g_sy,
            g_th,
            w_tooth,
            l_st,
            alpha_q,
            r_si,
            r_so,
            r_sy,
            k_ins,
            w_ins,
            k_fe,
            h,
            alpha_slot,
            Q_coil,
            h_slot,
            T_ref
        )
    ana = sta.StatorThermalAnalyzer()
    results = ana.analyze(problem)
    T_coil_vect[ind]=results['Coil temperature']

fig,ax=plt.subplots(1,1)
ax.plot(l_tooth_vect,T_coil_vect)
ax.set_xlabel('Stator tooth length [m]')
ax.set_ylabel('Coil temperature [C]')