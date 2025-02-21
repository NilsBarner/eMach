import numpy as np
import scipy.optimize as op
import os
import sys

# add the directory  this file's directory to path for module import
sys.path.append(os.path.dirname(__file__))

import thermal_network as tb


class SPM_RotorThermalProblem:

    def __init__(self,
        rte,
        n_ps,
        r_wire,
        w_tooth,
        h_tooth,
        t_sb,
        lslot,
        rbsrt,
        rrim,
        g_rt,
        t_m,
        t_rb,
        ttrt,
        trim,
        Ku,
        Adi_per_str,
        A_str,
        Vol_nd_insul,
        Vol_nd_potting,
        pwp,
        P,
        pcp,
        pdp,
        l_rte,
        tmrt,
        omm,
        w_slot,
        h_slot,
        
        T_ref,
        losses
    ):
        self.rte = rte
        self.n_ps = n_ps
        self.r_wire = r_wire
        self.w_tooth = w_tooth
        self.h_tooth = h_tooth
        self.t_sb = t_sb
        self.lslot = lslot
        self.rbsrt = rbsrt
        self.rrim = rrim
        self.g_rt = g_rt
        self.t_m = t_m
        self.t_rb = t_rb
        self.ttrt = ttrt
        self.trim = trim
        self.Ku = Ku
        self.Adi_per_str = Adi_per_str
        self.A_str = A_str
        self.Vol_nd_insul = Vol_nd_insul
        self.Vol_nd_potting = Vol_nd_potting
        self.pwp = pwp
        self.P = P
        self.pcp = pcp
        self.pdp = pdp
        self.l_rte = l_rte
        self.tmrt = tmrt
        self.omm = omm
        self.w_slot = Aslot / (ttrt * rte)
        self.h_slot = t_ws
        
        self.T_ref = T_ref
        self.losses = losses
        self.omega = omm


class SPM_RotorThermalAnalyzer:

    def __init__(self):
        self.base_ana = tb.ThermalNetworkAnalyzer()

    def analyze(self, problem: SPM_RotorThermalProblem):

        N_nodes = 16
        Res = self.create_resistance_network(problem)

        ################################################
        #           Load Losses into loss Vector
        ################################################
        Q_dot = np.zeros([N_nodes, 1])
        Q_dot[0] = 0
        Q_dot[1] = 0
        Q_dot[2] = problem.losses["core"]
        Q_dot[3] = problem.losses["copper"]
        Q_dot[4] = problem.losses["core"]
        Q_dot[5] = problem.losses["copper"]
        Q_dot[6] = problem.losses["core"]
        Q_dot[7] = problem.losses["copper"]
        Q_dot[8] = problem.losses["core"]
        Q_dot[9] = problem.losses["core"]
        Q_dot[10] = 0
        Q_dot[11] = problem.losses["magnet"]
        Q_dot[12] = problem.losses["magnet"]
        Q_dot[13] = 0
        Q_dot[14] = 0
        Q_dot[15] = 0

        ################################################
        #    Create Reference Temperature Vector
        ################################################

        T_ref = [
            [0, problem.T_ref],
        ]

        base_prob = tb.ThermalNetworkProblem(Res, Q_dot, T_ref, N_nodes)
        T = self.base_ana.analyze(base_prob)
        return T

    def create_resistance_network(self, problem):
        ################################################
        #             Load Operating Point
        ################################################
        Omega = problem.omega

        ################################################
        #           Import Geometric Values
        ################################################
        # Widths
        w_slot = self.w_slot  # slot width
        w_bundle = 2 * self.r_wire  # (ASSUME CIRCULAR CROSS-SECTION)
        w_tooth = self.w_tooth
        w_so = self.w_tooth  # tooth opening width (NEGLECT CURVATURE FOR NOW)
        
        # Heights
        h_slot = self.h_slot  # slot height (ASSUME FLUSH FOR NOW)
        h_bundle = 2 * self.r_wire  # (ASSUME CIRCULAR CROSS-SECTION)
        h_tt = self.h_tooth  # tooth tip height (ASSUME FLUSH FOR NOW)
        
        # Lengths
        # l_endturn = None  # (NEGLECT FOR NOW)
        l = self.lslot  # (IGNORE END WINDINGS FOR NOW)
        nz = 1  # ASSUME SINGLE AXIAL DISCRETISATION FOR NOW
        l_div_nz = l/nz  # axial length of node
        
        # Thicknesses
        t_sbi = self.t_sb  # stator back iron mean thickness
        t_separator = 7 * 25.4e-6  # Figure 5-11 in Aidan's thesis
        trim = self.trim
        t_nomex = 7 * 25.4e-6  # Figure 5-11 in Aidan's thesis
        t_magnet = self.t_m
        g = self.g_rt * self.rte  # air gap thickness (page 27)
        
        # Radii
        r_t = self.rte  # electromagnetic tip radius
        r_sbi = self.rbsrt * self.rte  # stator back iron mean radius
        r_ag = self.rte  # stator airgap surface radius
        r_mm_node = r_ag + self.g_rt * self.rte + t_magnet/2  # mean radius of magnet node
        r_rimm_node = self.rrim
        r_ro = self.rte + self.g_rt * self.rte + self.t_m + self.t_rb  # rotor outer surface radius (see 'Outer Radius')
        r_hsis = self.rte - self.ttrt * self.rte - self.t_sb  # outer surface of lattice structure (and inner surface of outer heat exchanger shell) (see 'Innermost Radius')

        # Nondimensionals
        K_cb = self.Ku
        K_ib = self.Adi_per_str/self.A_str
        n_slots = self.n_ps
        # Slot fractions of copper, insulation, potting, slot separator, and nomex slot lining
        v_cu = self.Ku
        v_insul = self.Ku * self.Vol_nd_insul
        v_potting = self.Ku * self.Vol_nd_potting
        v_separator = 1e-30  # (NEGLECT FOR NOW)
        v_nomex = 1e-30  # (NEGLECT FOR NOW)

        ################################################
        #           Create Resistance Objects
        ################################################
        Resistances = []
        ##############
        # Path 12
        ##############
        Descr = "Shaft to rotor core interface (Approximated as Plane Wall)"
        d = 1e-3  # mean strut thickness (GUESS FOR NOW)
        a = 5e-3  # lattice cell vertex length  (GUESS FOR NOW)
        t_lattice = 5e-3  # radial thickness of lattice structure (GUESS FOR NOW)
        rho_lattice = 3 * np.pi/4 * (d/a)**2 + (1 - 3 * np.pi/4) * (d/a)**3  # (C.4)
        A_flow = np.pi * (r_hsis**2 - (r_hsis - t_lattice)**2) * (1 - gams_math.sqrt(rho_lattice))  # (C.5)
        rho_air_in = 1.225
        u_lattice = mdot_hex / (rho_air_in * A_flow)  # page 228
        nu = 1.48e-5
        Re_lattice = u_lattice * d/nu  # (C.6)
        Pr = 0.71
        Nu = Pr**0.36 * 0.033 * Re_lattice**0.8  # (C.7) IF-ELSE LOGIC TO BE INTEGRATED INTO GAMSPY!
        k_hex = 237
        h = k_air * Nu/d  # page 228
        Bi = h * d/k_hex  # page 228
        h_lattice = h * (2 * rho_lattice * Bi**(-0.5)) * gams_math.tanh(2 * Bi**0.5 * t_lattice/d)  # (C.8)
        V_lattice = np.pi * (r_hsis**2 - (r_hsis - t_lattice)**2) * l  # (C.12)
        alpha_A = 2 * gams_math.sqrt(3 * np.pi) / a * gams_math.sqrt(rho_lattice)  # (C.13)
        A_surf = alpha_A * V_lattice  # (C.11)
        t_shell = 1e-3  # page 225
        kappa = 8  # page 229
        f_lattice = kappa * (0.44 + 0.008 * a/d) / (a/d - 1)**(0.43 + 1.13 * d/a) * Re_lattice**(-0.15)  # (C.15)
        delta_p_lattice = 4 * l/a * f_lattice * rho_air_in * u_lattice**2  # (C.16)
        Resistances.append(tb.hex(
            1, 2,
            h_lattice, k_hex, mdot_hex, cp_air,
            r_hsis, n_slots, l_div_nz, t_shell, l
        ))
        Resistances[0].Descr = Descr
        ##############
        # Path 23
        ##############
        Descr = "Shaft/RC interface to rotor core center"
        Resistances.append(tb.hexrim_sbi(
            2, 3,
            k_core, k_hex, cp_air,
            t_shell, t_sbi, r_hsis, r_sbi, n_slots, l_div_nz, l
        ))
        Resistances[1].Descr = Descr
        ##############
        # Path 34
        ##############
        Descr = "Rotor core center to PM/RC interface"
        k_winding = k_cu  # page 83
        Resistances.append(tb.winding_sbi(
            3, 4,
            _bundle, k_winding_tv, k_core,
            w_bundle, l_div_nz, t_nomex, w_slot, t_sbi, l_div_nz
        ))
        Resistances[2].Descr = Descr
        ##############
        # Path 35
        ##############
        Descr = "PM/RC interface to PM center"
        Resistances.append(tb.tooth_sbi(
            3, 5,
            h_rt, k_core,
            r_sbi, n_slots, w_slot, l_div_nz, t_sbi
        ))
        Resistances[3].Descr = Descr
        ##############
        # Path 45
        ##############
        Descr = "PM center to PM/sleeve Interface"
        Resistances.append(tb.winding_tooth(
            4, 5,
            k_winding, k_nomex, k_core,
            w_bundle, h_bundle, l_div_nz, t_nomex, h_rt, w_tooth, h_rt
        ))
        Resistances[4].Descr = Descr
        ##############
        # Path 46
        ##############
        Descr = "PM/Sleeve Interface to sleeve center"
        K_pb = 1 - K_cb - K_ib  # page 214
        k_amalgam = k_potting * K_pb / (K_pb + K_ib) + k_insul * K_ib / (K_pb + K_ib)  # (B.18)
        k_winding_tv = k_amalgam * ((1 + K_cb) * k_cu + (1 - K_cb) * k_amalgam) / ((1 - K_cb) * k_cu + (1 + K_cb) * k_amalgam)  # (B.18)
        Resistances.append(tb.winding_rad(
            4, 6,
            k_amalgam, k_winding_tv, k_separator,
            K_pb, h_bundle, w_bundle, l_div_nz, t_separator
        ))
        Resistances[5].Descr = Descr
        ##############
        # Path 57
        ##############
        Descr = "Sleeve center to outer rotor edge"
        k_core = 73  # iron at room temperature
        Resistances.append(tb.tooth_tooth(
            5, 7,
            k_core,
            h_rt, k_core, w_tooth, l_div_nz
        ))
        Resistances[6].Descr = Descr
        ##############
        # Path 67
        ##############
        Descr = "Outer rotor edge to air"
        Resistances.append(tb.R_winding_tooth(
            6, 7,
            k_winding, k_nomex, k_core,
            w_bundle, h_bundle, l_div_nz, t_nomex, h_rt, w_tooth, h_rt
        )
        # print(Resistances[7].h)
        Resistances[7].Descr = Descr
        ##############
        # Path 68
        ##############
        Descr = "Rotor core center to Hub/RotorCore Interface"
        Resistances.append(tb.winding_rad(
            6, 8,
            k_amalgam, k_winding_tv, k_separator,
            K_pb, h_bundle, w_bundle, l_div_nz, t_separator
        ))
        Resistances[8].Descr = Descr
        ##############
        # Path 79
        ##############
        Descr = "PM center to Hub/PM Interface"
        Resistances.append(tb.tooth_tooth(
            7, 9,
            k_core,
            h_rt, k_core, w_tooth, l_div_nz
        ))
        Resistances[9].Descr = Descr
        ##############
        # Path 89
        ##############
        Descr = "Sleeve center to Hub/Sleeve Interface"
        Resistances.append(tb.R_winding_tooth(
            8, 9,
            k_winding, k_nomex, k_core,
            w_bundle, h_bundle, l_div_nz, t_nomex, h_rt, w_tooth, h_rt
        )
        Resistances[10].Descr = Descr
        ##############
        # Path 910
        ##############
        Descr = "Shaft Center to shaft Inline with Hub center"
        Resistances.append(tb.tooth_tt(
            9, 10,
            k_core,
            h_rt, r_ag, h_tt, h_rt, n_slots, w_slot, l_div_nz, w_tt
        ))
        Resistances[11].Descr = Descr
        ##############
        # Path 811
        ##############
        Descr = "Hub/Rotor Core interface to Center of Hub inline with Rotor Core"
        rho_a = 1.225
        Ma_air_in = 0.3  # page 84
        gamma_air = 1.3
        R_air = 287
        u_air_in = Ma_air_in * gams_math.sqrt(gamma_air * R_air * T_air_in)
        A_flow_ag = np.pi * ((r_t + g)**2 - r_t**2)  # (B.36)
        mdot_ag = rho_a * A_flow_ag * u_air_in  # (B.36)
        u_ag = 0.5 * Wdot_windage_inner / (Omega * mdot_ag * (r_t + g/2))  # (B.37)
        dH_ag = 4 * A_flow_ag / (2 * np.pi * r_t)  # (B.38)
        Re_ag_z = u_air_in * dH_ag/nu_a  # (B.38)
        Nu_ag = 0.022 * (1 + (dH_ag * u_ag / (np.pi * r_t * u_air_in)**2))**0.8714 * Re_ag_z**0.8 * Pr_a**0.5  # (B.38) NOTE: there is a bracket missing - I suspect this is the correct form by comparison with other Nusselt correlations but might have to confirm by looking at the OOM of the result
        h_ag = Nu_ag * k_air/dH_ag  # (B.38)
        # k_hs = None  # NEGLECT IN ABSENCE OF INFORMATION ABOUT k_hs
        cp_air = 1005
        Resistances.append(tb.winding_ag(
            8, 11,
            k_nomex, k_winding_tv, h_ag, mdot_ag, cp_air,
            h_rt, w_slot, l_div_nz, t_nomex, w_so, n_slots
        ))
        Resistances[12].Descr = Descr
        ##############
        # Path 1011
        ##############
        Descr = "Hub/PM interface to Center of Hub inline with PM"
        Resistances.append(tb.tt_ag(
            10, 15,
            k_core, h_ag, mdot_ag, cp_air,
            h_tt, r_ag, n_slots, w_so, l_div_nz, r_t
        ))
        Resistances[13].Descr = Descr
        ##############
        # Path 1112
        ##############
        Descr = "Hub/Sleeve interface to Center of Hub inline with Sleeve"
        Resistances.append(tb.magnet_ag(
            11, 16,
            k_magnet, h_ag, mdot_ag, cp_air,
            h_rm, r_t, g, n_slots, l_div_nz, r_t
        ))
        Resistances[14].Descr = Descr
        ##############
        # Path 1213
        ##############
        Descr = "Shaft inline with Hub center to Hub/Shaft interface"
        k_magnet = 13
        n_rm = 2  # number of radial magnet divisions (see Figure B-1)
        t_magnet_node = t_magnet/n_rm  # thickness of magnet node
        Resistances.append(tb.magnet_rad(
            12, 13,
            k_magnet,
            r_mm_node, t_magnet_node, n_slots, l_div_nz
        ))
        Resistances[15].Descr = Descr
        ##############
        # Path 1314
        ##############
        Descr = "Hub/Shaft interface to Center of Hub inline with Rotor Core"
        t_epoxy = 0.5e-3  # page 217
        k_epoxy = 0.14
        Resistances.append(tb.magnet_rim(
            13, 14,
            k_magnet, k_epoxy, k_rim,
            r_mm_node, t_magnet_node, n_slots, l_div_nz, r_t, g, r_rimm_node, t_rim_node
        ))
        Resistances[16].Descr = Descr
        ##############
        # Path 1415
        ##############
        Descr = "Hub inline with Rotor Core to Hub inline with PM"
        n_rrim = 2  # number of radial rim divisions (see Figure B-1) - I ADDED THIS VARIABLE
        t_rim_node = trim/n_rrim
        k_rim = 7  # Titanium Ti-6Al-4V Alloy (page 83)
        Resistances.append(tb.rim_rad(
            14, 15,
            k_rim,
            r_rimm_node, t_rim_node, n_slots, l_div_nz
        ))
        Resistances[17].Descr = Descr
        ##############
        # Path 1516
        ##############
        Descr = "Hub inline with PM to Hub inline with Sleeve "
        nu_a = 1.48e-5  # kinematic viscosity of air
        k_air = 0.03
        Pr_a = 0.71
        Re_ro = Omega * r_ro**2/nu_a  # (B.35)
        Nu_ro = 0.0296 * Re_ro**0.8 * Pr_a  # (B.35)
        h_ro = Nu_ro * k_air / (2 * r_ro)  # (B.35)
        Resistances.append(tb.rim_amb(
            15, 16
            k_rim, h_ro,
            r_rimm_node, r_ro, n_slots, l_div_nz
        ))
        Resistances[18].Descr = Descr
        
        return Resistances


if __name__ == "__main__":
    rte = 0.128
    n_ps = 60
    r_wire = 5e-3
    w_tooth = 6e-3
    h_tooth = 1e-3
    t_sb = 0.128 - 0.105 - 16e-3
    lslot = 0.198
    rbsrt = (0.128 - 16e-3 - t_sb/2) / rte
    rrim = 0.15 - 8e-3/2
    g_rt = 3e-3/rte
    t_m = 10e-3
    t_rb = 8e-3
    ttrt = 16e-3/rte
    trim = t_rb
    Ku = 0.5
    Adi_per_str = np.pi*(0.321/2/1000)**2 * 0.8  # guess
    A_str = np.pi*(0.321/2/1000)**2
    Vol_nd_insul = 0.1  # guess
    Vol_nd_potting = 0.1  # guess
    pwp = 6342/1e6  # Table 5.5
    P = 1e6
    pcp = 8820/1e6  # Table 5.5
    pdp = 1955 + 2382  # Table 5.5
    l_rte = l_slot/rte
    tmrt = t_m/rte
    omm = 12500 * np.pi/30
    w_slot = 6e-3
    h_slot = 16e-3
    
    T_ref = 25 + 273.15
    losses = {
        "core": pcp * P,
        "copper" : pdp * P,
        "magnet" : 0,
        "windage": pwp * P
    }
    
    prob=SPM_RotorThermalProblem(
        rte,
        n_ps,
        r_wire,
        w_tooth,
        h_tooth,
        t_sb,
        lslot,
        rbsrt,
        rrim,
        g_rt,
        t_m,
        t_rb,
        ttrt,
        trim,
        Ku,
        Adi_per_str,
        A_str,
        Vol_nd_insul,
        Vol_nd_potting,
        pwp,
        P,
        pcp,
        pdp,
        l_rte,
        tmrt,
        omm,
        w_slot,
        h_slot,
        
        T_ref,
        losses
    )
    
    ana=SPM_RotorThermalAnalyzer()
