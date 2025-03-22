#%% 7.2. Input from User

import numpy as np
from eMach.mach_eval.analyzers.electromagnetic.winding_factors import (
    WindingFactorsProblem,
    WindingFactorsAnalyzer,
    )

n = np.array([1,2,3,4,5])
# winding_layout = np.array([[0,0,0,0,-1,-1,0,0,0,0,1,1],[0,0,-1,-1,0,0,0,0,1,1,0,0]])
winding_layout = np.array([[0,0,0,0,-1,-1,0,0,0,0,1,1],[0,0,-1,-1,0,0,0,0,1,1,0,0],[-1,-1,0,0,0,0,1,1,0,0,0,0]])

# winding_layout = np.array([[1, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, -1, -1, 0, 0, 0, 0]])
# winding_layout = np.array([[1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0]])

winding_layout = np.array([[1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
alpha_1 = 0

#%% 7.3. Output to User

n = 3
p = 4
print(abs(np.exp(-1j*0) - np.exp(-1j * 2 * np.pi * p/9) + np.exp(-1j * 4 * np.pi * p/9)) / n)
print((2 * np.cos(2 * np.pi * p/9) - 1) / 3)

alpha_1 = 8*np.pi/9
kw_prob = WindingFactorsProblem(n,winding_layout,alpha_1)
kw_ana = WindingFactorsAnalyzer()

k_w = kw_ana.analyze(kw_prob)
print(k_w)
print(np.abs(k_w))

#%%

import sys

N_phases = 3
N_coils = 9
N_pole_pairs = 4
# N_pole_pairs = 5
mech_to_electr_angle = N_pole_pairs
multiplier = 2*np.pi/N_coils
signs = np.array([1, -1, 1])
mech_angles_phase_coils = np.array([0 * multiplier, 1 * multiplier, 2 * multiplier])
electr_angles_phase_coils = mech_angles_phase_coils * mech_to_electr_angle
# print(electr_angles_phase_coils)
# print(-2 * np.pi * 4/9, 4 * np.pi * 4/9)

k_w = np.sum(signs * np.exp(-1j * electr_angles_phase_coils))/N_phases
k_w_abs = np.abs(k_w)
print('k_w_abs =', k_w_abs)

# sys.exit()

mech_angles_phase_coils_2 = mech_angles_phase_coils + np.pi
electr_angles_phase_coils_2 = mech_angles_phase_coils_2 * mech_to_electr_angle
k_w_2 = np.sum(signs * np.exp(-1j * electr_angles_phase_coils_2))/N_phases
k_2_final = k_w + k_w_2
k_2_final_abs = np.abs(k_2_final)
print('k_2_final_abs =', k_2_final_abs)

#%%

# Calculate winding factor
N_phases = 3  # number of phases (Theepan's email from 17.01.2025)
n = np.array([1])  # harmonics to be calculated

N_layers = N_phases
N_poles * 6 / N_phases
N_slots_per_pole = 6

winding_layout = np.zeros([N_phases, N_poles * N_slots_per_pole/N_phases])
for i in range(len(N_layers)):
    winding_layout[i][i+0:i+2] = [-1, -1]
    winding_layout[i][i+6:i+8] = [1, 1]

winding_layout = np.array([
    [0,0,0,0,-1,-1,0,0,0,0,1,1],
    [0,0,-1,-1,0,0,0,0,1,1,0,0],
    [-1,-1,0,0,0,0,1,1,0,0,0,0]
])

#%%

N_poles = 8


N_phases = 3  # number of phases (Theepan's email from 17.01.2025)
n = np.array([1])  # harmonics to be calculated

N_layers = N_phases
# N_poles * 6 / N_phases
N_slots_per_pole = 6

# winding_layout = np.zeros([N_phases, 12])
winding_layout = np.zeros([N_phases, int(N_poles * N_slots_per_pole/N_phases)])
for i in range(N_layers):
    winding_layout[i][2*i+0:2*i+2] = [-1, -1]
    winding_layout[i][2*i+6:2*i+8] = [1, 1]
print(winding_layout)

#%%

import numpy as np

N_poles = 8
N_phases = 3  # number of phases
N_layers = N_phases
N_slots_per_pole = 6

# Define winding layout
N_slots_per_phase = int(N_poles * N_slots_per_pole / N_phases)
winding_layout = np.zeros([N_phases, N_slots_per_phase])

for i in range(N_layers):
    # Wrap around indices using modulo operation
    winding_layout[i][(2 * i + 0) % N_slots_per_phase] = -1
    winding_layout[i][(2 * i + 1) % N_slots_per_phase] = -1
    winding_layout[i][(2 * i + 6) % N_slots_per_phase] = 1
    winding_layout[i][(2 * i + 7) % N_slots_per_phase] = 1

print(winding_layout)
