class Resistance:
    """Base class for thermal resisance

    Attributes:
        Material: Material object holding material properties
        Node1: First node connected to resistance
        Node2: Second node connected to resistance
        resistance_value: Thermal resistance [K/W]
    """

    def __init__(self, Material: Material, Node1: int, Node2: int):
        self.Material = Material
        self.Node1 = Node1
        self.Node2 = Node2

    @property
    def resistance_value(self):
        return None
        

class hex(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            h_lattice, k_hex, mdot_hex, cp_air,
            r_hsis, n_slots, l_div_nz, t_shell, l
        ):
        super().__init__(None, Node1, Node2)
        self.h_lattice = h_lattice
        self.k_hex = k_hex
        self.mdot_hex = mdot_hex
        self.cp_air = cp_air
        self.r_hsis = r_hsis
        self.n_slots = n_slots
        self.l_div_nz = l_div_nz
        self.t_shell = t_shell
        self.l = l

    @property
    def resistance_value(self):
        R_lattice = 1 / (self.h_lattice * 2 * np.pi * self.r_hsis/self.n_slots * self.l_div_nz)  # (C.9)
        return R_lattice + self.t_shell / (self.k_hex * 2 * np.pi * self.r_hsis/self.n_slots * self.l) + 1 / (2 * self.mdot_hex * self.cp_air)  # (C.14)
        
        
class hexrim_sbi(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_core, k_hex, cp_air,
            t_shell, t_sbi, r_hsis, r_sbi, n_slots, l_div_nz, l
        ):
        super().__init__(None, Node1, Node2)
        self.k_core = k_core
        self.k_hex = k_hex
        self.cp_air = cp_air
        self.t_shell = t_shell
        self.t_sbi = t_sbi
        self.r_hsis = r_hsis
        self.r_sbi = r_sbi
        self.n_slots = n_slots
        self.l_div_nz = l_div_nz
        self.l = l
        
    @property
    def resistance_value(self):
        return t_shell / (k_hex * 2 * np.pi * r_hsis/n_slots * l) + t_sbi/2 / (k_core * 2 * np.pi * r_sbi/n_slots * l_div_nz)
        
        
class winding_sbi(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            h_bundle, k_winding_tv, k_core,
            w_bundle, l_div_nz, t_nomex, w_slot, t_sbi, l_div_nz
        ):
        super().__init__(None, Node1, Node2)
        self.h_bundle = h_bundle
        self.k_winding_tv = k_winding_tv
        self.k_core = k_core
        self.w_bundle = w_bundle
        self.l_div_nz = l_div_nz
        self.t_nomex = t_nomex
        self.w_slot = w_slot
        self.t_sbi = t_sbi
        
    @property
    def resistance_value(self):
        return h_bundle/2 / (k_winding_tv * w_bundle * l_div_nz) + \
            t_nomex / (k_nomex * w_slot * l_div_nz) + \
            t_sbi/2 / (k_core * w_slot * l_div_nz)  # (B.30)
            
            
class tooth_sbi(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            h_rt, k_core,
            r_sbi, n_slots, w_slot, l_div_nz, t_sbi
        ):
        super().__init__(None, Node1, Node2)
        self.h_rt = h_rt
        self.k_core = k_core
        self.n_slots = n_slots
        self.w_slot = w_slot
        self.l_div_nz = l_div_nz
        self.t_sbi = t_sbi
        
    @property
    def resistance_value(self):
        return 0.5 * (h_rt / ((2 * np.pi * (r_sbi + h_rt/2) / n_slots - w_slot) * k_core * l_div_nz) + t_sbi / (k_core * (2 * np.pi * r_sbi/n_slots - w_slot) * l_div_nz))  # (B.32)
        
        
class winding_tooth(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_winding, k_nomex, k_core,
            w_bundle, h_bundle, l_div_nz, t_nomex, h_rt, w_tooth, h_rt
        ):
        super().__init__(None, Node1, Node2)
        self.k_winding = k_winding
        self.k_nomex = k_nomex
        self.k_core = k_core
        self.w_bundle = w_bundle
        self.h_bundle = h_bundle
        self.l_div_nz = l_div_nz
        self.t_nomex = t_nomex
        self.h_rt = h_rt
        self.w_tooth = w_tooth
        self.h_rt = h_rt
        
    @property
    def resistance_value(self):
        return 0.5 * (self.w_bundle/2 / (self.k_winding * self.h_bundle * self.l_div_nz) + self.t_nomex / (self.k_nomex * self.h_rt * self.l_div_nz) + 0.5 * self.w_tooth / (self.k_core * self.h_rt * self.l_div_nz))  # (B.28)
        
        
class winding_rad(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_amalgam, k_winding_tv, k_separator,
            K_pb, h_bundle, w_bundle, l_div_nz, t_separator
        ):
        super().__init__(None, Node1, Node2)
        self.k_amalgam = k_amalgam
        self.k_winding_tv = k_winding_tv
        self.k_separator = k_separator
        self.K_pb = K_pb
        self.h_bundle = h_bundle
        self.w_bundle = w_bundle
        self.l_div_nz = l_div_nz
        self.t_separator = t_separator
    
    @property
    def resistance_value(self):
        return self.h_bundle / (self.k_winding_tv * self.w_bundle * self.l_div_nz) + self.t_separator / (self.k_separator * self.w_bundle * self.l_div_nz)  # (B.19)
        
        
class tooth_tooth(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_core,
            h_rt, k_core, w_tooth, l_div_nz
        ):
        super().__init__(None, Node1, Node2)
        self.k_core = k_core
        self.h_rt = h_rt
        self.k_core = k_core
        self.w_tooth = w_tooth
        self.l_div_nz = l_div_nz
    
    @property
    def resistance_value(self):
        return self.h_rt / (self.k_core * self.w_tooth * self.l_div_nz)  # (B.21)
        
        
class tooth_tt(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_core,
            h_rt, r_ag, h_tt, h_rt, n_slots, w_slot, l_div_nz, w_tt
        ):
        super().__init__(None, Node1, Node2)
        self.k_core = k_core
        self.h_rt = h_rt
        self.r_ag = r_ag
        self.h_tt = h_tt
        self.h_rt = h_rt
        self.n_slots = n_slots
        self.w_slot = w_slot
        self.l_div_nz = l_div_nz
        self.w_tt = w_tt
    
    @property
    def resistance_value(self):
        return 0.5 * (self.h_rt / ((2 * np.pi * (self.r_ag - self.h_tt - self.h_rt/2) / self.n_slots - self.w_slot) * self.k_core * self.l_div_nz) + self.h_tt / (self.k_core * self.w_tt * self.l_div_nz))  # (B.31) NOTE THEY TYPO - compare with (B.32)
        

class winding_ag(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_nomex, k_winding_tv, h_ag, mdot_ag, cp_air,
            h_rt, w_slot, l_div_nz, t_nomex, w_so, n_slots
        ):
        super().__init__(None, Node1, Node2)
        self.k_nomex = k_nomex
        self.k_winding_tv = k_winding_tv
        self.h_ag = h_ag
        self.mdot_ag = mdot_ag
        self.cp_air = cp_air
        self.h_rt = h_rt
        self.w_slot = w_slot
        self.l_div_nz = l_div_nz
        self.t_nomex = t_nomex
        self.w_so = w_so
        self.n_slots = n_slots
    
    @property
    def resistance_value(self):
        return h_rt/2 / (k_winding_tv * w_slot * l_div_nz) + \
            self.t_nomex / (self.k_nomex * self.w_so * self.l_div_nz) + \
            1 / (self.h_ag * self.w_so * self.l_div_nz) + \
            1 / (2 * self.mdot_ag/self.n_slots * self.cp_air)  # (B.39)
            # + h_tt / (k_hs * w_so * l_div_nz) + \  # NEGLECT IN ABSENCE OF INFORMATION ABOUT k_hs
            

class tt_ag(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_core, h_ag, mdot_ag, cp_air,
            h_tt, r_ag, n_slots, w_so, l_div_nz, r_t
        ):
        super().__init__(None, Node1, Node2)
        self.k_core = k_core
        self.h_ag = h_ag
        self.mdot_ag = mdot_ag
        self.cp_air = cp_air
        self.h_tt = h_tt
        self.r_ag = r_ag
        self.n_slots = n_slots
        self.w_so = w_so
        self.l_div_nz = l_div_nz
        self.r_t = r_t
    
    @property
    def resistance_value(self):
        return self.h_tt/2 / (self.k_core * (2 * np.pi * self.r_ag/self.n_slots - self.w_so) * self.l_div_nz) + \
            1 / (self.h_ag * (2 * np.pi * self.r_t/self.n_slots - self.w_so) * self.l_div_nz) + \
            1 / (2 * self.mdot_ag/self.n_slots * self.cp_air)  # (B.40)
            
            
class magnet_ag(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_magnet, h_ag, mdot_ag, cp_air,
            h_rm, r_t, g, n_slots, l_div_nz, r_t
        ):
        super().__init__(None, Node1, Node2)
        self.k_magnet = k_magnet
        self.h_ag = h_ag
        self.mdot_ag = mdot_ag
        self.cp_air = cp_air
        self.h_rm = h_rm
        self.r_t = r_t
        self.g = g
        self.n_slots = n_slots
        self.l_div_nz = l_div_nz
        self.r_t = r_t
    
    @property
    def resistance_value(self):
        return self.h_rm/2 / (self.k_magnet * 2 * np.pi * (self.r_t + self.g) / self.n_slots * self.l_div_nz) + \
            1 / (self.h_ag * 2 * np.pi * (self.r_t + self.g) / self.n_slots * self.l_div_nz) + \
            1 / (2 * self.mdot_ag/self.n_slots * self.cp_air)  # (B.41)
            
            
class magnet_rad(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_magnet,
            r_mm_node, t_magnet_node, n_slots, l_div_nz
        ):
        super().__init__(None, Node1, Node2)
        self.k_magnet = k_magnet
        self.r_mm_node = r_mm_node
        self.t_magnet_node = t_magnet_node
        self.n_slots = n_slots
        self.l_div_nz = l_div_nz
    
    @property
    def resistance_value(self):
        return np.log((self.r_mm_node + self.t_magnet_node/2) / (self.r_mm_node - self.t_magnet_node/2)) / (2 * np.pi/self.n_slots * self.k_magnet * self.l_div_nz)  # (B.25)
        
        
class magnet_rim(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_magnet, k_epoxy, k_rim,
            r_mm_node, t_magnet_node, n_slots, l_div_nz, r_t, g, r_rimm_node, t_rim_node
        ):
        super().__init__(None, Node1, Node2)
        self.k_magnet = k_magnet
        self.k_epoxy = k_epoxy
        self.k_rim = k_rim
        self.r_mm_node = r_mm_node
        self.t_magnet_node = t_magnet_node
        self.n_slots = n_slots
        self.l_div_nz = l_div_nz
        self.r_t = r_t
        self.g = g
        self.r_rimm_node = r_rimm_node
        self.t_rim_node = t_rim_node
    
    @property
    def resistance_value(self):
        return gams_math.log((self.r_mm_node + self.t_magnet_node/2) / self.r_mm_node) / (self.k_magnet * 2 * np.pi/self.n_slots * self.l_div_nz) + \
            self.t_epoxy / (self.k_epoxy * 2 * np.pi * (self.r_t + self.g + self.t_magnet) / self.n_slots * self.l_div_nz) + \
            gams_math.log(self.r_rimm_node / (self.r_rimm_node - self.t_rim_node/2)) / (self.k_rim * 2 * np.pi/self.n_slots * self.l_div_nz)  # (B.33)
            
            
class rim_rad(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_rim,
            r_rimm_node, t_rim_node, n_slots, l_div_nz
        ):
        super().__init__(None, Node1, Node2)
        self.k_rim = k_rim
        self.r_rimm_node = r_rimm_node
        self.t_rim_node = t_rim_node
        self.n_slots = n_slots
        self.l_div_nz = l_div_nz
        
    @property
    def resistance_value(self):
        return gams_math.log((self.r_rimm_node + self.t_rim_node/2) / (self.r_rimm_node - self.t_rim_node/2)) / (2 * np.pi/self.n_slots * self.k_rim * self.l_div_nz)  # (B.27)
        
        
class rim_amb(Resistance):
    
    def __init__(self,
            Node1: int, Node2: int,
            k_rim, h_ro,
            r_rimm_node, r_ro, n_slots, l_div_nz
        ):
        super().__init__(None, Node1, Node2)
        self.k_rim = k_rim
        self.h_ro = h_ro
        self.r_rimm_node = r_rimm_node
        self.r_ro = r_ro
        self.n_slots = n_slots
        self.l_div_nz = l_div_nz
        
    @property
    def resistance_value(self):
        return gams_math.log(self.r_ro/self.r_rimm_node) / (self.k_rim * 2 * np.pi/self.n_slots * self.l_div_nz) + 1 / (self.h_ro * 2 * np.pi * self.r_ro/self.n_slots * self.l_div_nz)  # (B.34)



    