import math

def ghk(R, T, F, P_K1, P_K2, P_K3, K_in, K_out, P_Na, Na_in, Na_out, P_Cl, Cl_in, Cl_out, P_Ca, Ca_in, Ca_out):
    """
    Computes membrane potential (Vm) in millivolts using the Goldman-Hodgkin-Katz equation,
    incorporating three potassium permeability values.
    
    Parameters:
        R     : float - Universal gas constant (J/(mol·K))
        T     : float - Temperature in Kelvin (K)
        F     : float - Faraday's constant (C/mol)
        P_K1  : float - Permeability of potassium (1st value)
        P_K2  : float - Permeability of potassium (2nd value)
        P_K3  : float - Permeability of potassium (3rd value)
        K_in  : float - Intracellular potassium concentration (mM)
        K_out : float - Extracellular potassium concentration (mM)
        P_Na  : float - Permeability of sodium
        Na_in : float - Intracellular sodium concentration (mM)
        Na_out: float - Extracellular sodium concentration (mM)
        P_Cl  : float - Permeability of chloride
        Cl_in : float - Intracellular chloride concentration (mM)
        Cl_out: float - Extracellular chloride concentration (mM)
        P_Ca  : float - Permeability of calcium
        Ca_in : float - Intracellular calcium concentration (mM)
        Ca_out: float - Extracellular calcium concentration (mM)
    
    Returns:
        Vm (float): Membrane potential in millivolts (mV)
    """
    
    # Calculate RT/F in millivolts
    V_T = (R * T) / F * 1000  # Convert to mV
    
    # Effective Potassium Permeability (weighted average or sum)
    P_K_eff = (P_K1 + P_K2 + P_K3) / 3
    
    numerator = (P_K_eff * K_out + P_Na * Na_out + P_Cl * Cl_in + P_Ca * Math.sqrt(Ca_out))
    denominator = (P_K_eff * K_in + P_Na * Na_in + P_Cl * Cl_out + P_Ca * Math.sqrt(Ca_in))
    
    if denominator == 0:
        raise ValueError("Denominator is zero; check input values.")
    
    Vm = V_T * math.log(numerator / denominator)
    return Vm

# Parameter Inputs
R = 8.314 # Constant Value
T = 294.15
F = 96485 # Constant Value

P_K1, P_K2, P_K3, K_in, K_out = 0.008, 0.01, 0.012, 140, 5
P_Na, Na_in, Na_out = 2, 15, 50
P_Cl, Cl_in, Cl_out = 0.45, 10, 120
P_Ca, Ca_in, Ca_out = 0.0001, 0.00001, 13.5

# Compute Membrane Potential
Vm = ghk(R, T, F, P_K1, P_K2, P_K3, K_in, K_out, P_Na, Na_in, Na_out, P_Cl, Cl_in, Cl_out, P_Ca, Ca_in, Ca_out)

print(f"\nMembrane potential: {Vm:.2f} mV")
