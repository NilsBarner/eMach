import sys
import numpy as np

def generate_winding_matrix(poles, slots_per_pole, pattern_type="adjacent"):
    """
    Generate a winding matrix for an electrical machine based on the given pattern type.
    
    Parameters:
    - poles (int): Number of poles
    - slots_per_pole (int): Slots per pole
    - pattern_type (str): "adjacent" for equal-phase windings immediately adjacent,
                          "spaced" for equal-phase windings evenly spaced
    
    Returns:
    - numpy array: A winding matrix where each row represents a phase (A, B, C)
    """
    slots = poles * slots_per_pole
    phases = 3  # Three-phase system
    winding_matrix = np.zeros((phases, slots), dtype=int)
    
    # Phase sequence: A, B, C
    phase_order = [0, 1, 2]  # A, B, C
    
    if pattern_type == "adjacent":
        # Ensure "red-yellow-blue" (A-B-C) pattern repeats within each pole
        for slot in range(slots):
            phase = phase_order[slot % phases]  # Cycle through A-B-C
            sign = (-1) ** (slot % 2)  # Alternate +1 and -1 for each adjacent slot
            winding_matrix[phase, slot] = sign
    
    elif pattern_type == "spaced":
        # Distribute phases evenly around the stator
        for slot in range(slots):
            phase = phase_order[(slot % (phases * slots_per_pole // phases)) // (slots_per_pole // phases)]
            sign = (-1) ** ((slot // (slots // poles)) % 2)  # Alternate every pole pair
            winding_matrix[phase, slot] = sign
    else:
        raise ValueError("Invalid pattern_type. Choose 'adjacent' or 'spaced'.")
    
    return winding_matrix

# Example usage
# poles = 2
# slots_per_pole = 6

poles = 8
slots_per_pole = 6

winding_adjacent = generate_winding_matrix(poles, slots_per_pole, pattern_type="adjacent")
winding_spaced = generate_winding_matrix(poles, slots_per_pole, pattern_type="spaced")

print("Winding pattern (adjacent phases):")
print(winding_adjacent)
print("\nWinding pattern (spaced phases):")
print(winding_spaced)

sys.exit()

#%%

from eMach.mach_eval.analyzers.electromagnetic.winding_factors import (
    WindingFactorsProblem,
    WindingFactorsAnalyzer,
    )

n = np.array([1,2,3,4,5])

alpha_1 = np.pi/12
kw_prob = WindingFactorsProblem(n,winding_adjacent,alpha_1)
# kw_prob = WindingFactorsProblem(n,winding_adjacent.reshape(1, -1),alpha_1)
kw_ana = WindingFactorsAnalyzer()
print('kw_prob =', abs(kw_ana.analyze(kw_prob)))
kw_prob = WindingFactorsProblem(n,winding_spaced,alpha_1)
# kw_prob = WindingFactorsProblem(n,winding_spaced.reshape(1, -1),alpha_1)
kw_ana = WindingFactorsAnalyzer()
print('kw_prob =', abs(kw_ana.analyze(kw_prob)))

