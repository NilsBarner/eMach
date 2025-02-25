#%% 1. Import modules and define geometry.

import sys
sys.path.append("C:/Users/nmb48/Documents/GitHub/eMach")  # Adjust path as needed
import eMach
import numpy as np
import scipy.optimize as op
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from eMach.mach_eval.analyzers.mechanical.thermal_network import *
#################
#Define Materials
#################

k_1=10 #Base Material Thermal Conductivity W/m-K
k_2=100 #Conductive Material Thermal Conductivity W/m-K
k_3=.01 #Insulating Material Thermal Conductivity W/m-K

mat1=Material(k_1)
mat2=Material(k_2)
mat3=Material(k_3)
##################
#Define Convection
##################
h=10 #Convection coefficient W/m^2-K
#################
#Define Geometry
#################

w=0.1 #Width m
L=.75 #Length m
d=.1 #Depth m

L_1=.5*L #Length of base
L_2=L*3/4 #Length to mid section 2 and 4
A_1=w*d #Cross sectional area of base
A_2=w*d/2 #Cross sectional area of section 2 and 3
A_3=(L-L_1)*d #Cross sectional Area between section 2 and 3

#%% Create Resistance objects for this example.

N_nodes=6 #Number of Nodes

###################
#Define Resistances
###################
Resistances = []
##############
# Path 0
##############
Descr = "R_1,2"
Resistances.append(plane_wall(mat1, 1, 2, L_1, A_1))
Resistances[0].Descr = Descr
##############
# Path 1
##############
Descr = "R_2,3"
Resistances.append(plane_wall(mat2, 2, 3, L_2-L_1, A_2))
Resistances[1].Descr = Descr

##############
# Path 2
##############
Descr = "R_2,4"
Resistances.append(plane_wall(mat3, 2, 4, L_2-L_1, A_2))
Resistances[2].Descr = Descr

##############
# Path 3
##############
Descr = "R_3,5"
Resistances.append(plane_wall(mat2, 3, 5, w/4, A_3))
Resistances[3].Descr = Descr

##############
# Path 4
##############
Descr = "R_4,5"
Resistances.append(plane_wall(mat3, 4, 5, w/4, A_3))
Resistances[4].Descr = Descr

##############
# Path 5
##############
Descr = "R_3,0"
Resistances.append(conv(None, 3, 0, h, A_2))
Resistances[5].Descr = Descr

##############
# Path 6
##############
Descr = "R_4,0"
Resistances.append(conv(None, 4, 0, h, A_2))
Resistances[6].Descr = Descr

#%% 3. Specify the heat sources at nodes 1 and 5.

####################
#Define Heat Sources
####################
Q_dot=[0,]*N_nodes #create a list of 0's of length N_nodes
Q_dot[1]=10 #W
Q_dot[5]=10 #W

#%% 4. Specify the temperature at the reference node. For this example, only one reference temperatures is used (at node 0).

######################
#Define Reference Temps
######################
ref_node=0
ref_temp=25 #C
T_ref=[[ref_node,ref_temp],]

#%% 5. Create the problem and analyzer.

############################
#Create Problem and Analzyer
############################
prob=ThermalNetworkProblem(Resistances,Q_dot,T_ref,N_nodes)
ana=ThermalNetworkAnalyzer()

#%% Output to the user

############################
#Analyze the Problem
############################
T=ana.analyze(prob)

############################
#Make Plot
############################
x=[L*1.2,0,L_1,L_2,L_2,L_2]
y=[0,0,0,w/4,-w/4,0]
fig,ax=plt.subplots(1,1)
c1=ax.scatter(x,y,c=T,s=200)
h=fig.colorbar(c1,label='Temperature')
# Create a Rectangle patch
rect = Rectangle((0,-w/2),L,w,linewidth=1,edgecolor='k',facecolor='none')
# Add the patch to the Axes
ax.add_patch(rect)
# Create a Rectangle patch
rect = Rectangle((L_1,0),L-L_1,w/2,linewidth=1,edgecolor='k',facecolor='none')
# Add the patch to the Axes
ax.add_patch(rect)
# Create a Rectangle patch
rect = Rectangle((L_1,-w/2),L-L_1,w/2,linewidth=1,edgecolor='k',facecolor='none')
# Add the patch to the Axes
ax.add_patch(rect)
ax.plot([x[1],x[2]],[y[1],y[2]],'r--')
ax.plot([x[2],x[3]],[y[2],y[3]],'r--')
ax.plot([x[2],x[3]],[y[2],y[4]],'r--')
ax.plot([x[3],x[5]],[y[3],y[5]],'r--')
ax.plot([x[4],x[5]],[y[4],y[5]],'r--')
ax.plot([x[3],x[0]],[y[3],y[0]],'r--')
ax.plot([x[4],x[0]],[y[4],y[0]],'r--')
ax.set_yticks([])
ax.set_xticks([])

#%% Check node and system enthalpy conservation equations

def check_node_enthalpy_conservation(problem: ThermalNetworkProblem): # Q_dot: list, Resistances: list):

    Q_tot_in = sum(problem.Q_dot)
    Q_tot_out = 0

    for node_idx in range(N_nodes):
        Q_gen = problem.Q_dot[node_idx]
    
        Q_flows = []
        for r in problem.res:
            
            if node_idx == r.Node1:
                # print(r.Descr)  # check whether correct resistances have been assigned
                Q_flow = (T[r.Node2] - T[r.Node1]) / r.resistance_value
                Q_flows.append(Q_flow)
            elif node_idx == r.Node2:
                # print(r.Descr)  # check whether correct resistances have been assigned
                Q_flow = (T[r.Node1] - T[r.Node2]) / r.resistance_value
                Q_flows.append(Q_flow)
            else:
                continue
        
        if node_idx not in [row[0] for row in problem.T_ref]:
        
            # Check node enthalpy conservation
            atol = sum(Q_flows) + Q_gen
            if np.isclose(atol, 0):
                print(f"✅ Node '{node_idx}' is balanced within absolute tolerance {atol} W.")
            else:
                print(f"⚠️ Node '{node_idx}' is unbalanced with absolute tolerance {atol} W.")
                
        else:
            Q_out = sum(Q_flows) + Q_gen
            print(f"{Q_out} W of heat is leaving the network through node {node_idx}.")
            Q_tot_out += Q_out
                
    # Check system enthalpy conservation
    atol = -Q_tot_in + Q_tot_out
    print()
    if np.isclose(atol, 0):
        print(f"✅ System is balanced within absolute tolerance {atol} W.")
    else:
        print(f"⚠️ System is unbalanced with absolute tolerance {atol} W.")
        
check_node_enthalpy_conservation(prob)















