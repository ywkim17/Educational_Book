import math

def ghk(R, T, F, P_K, K_in, K_out, P_Na, Na_in, Na_out, P_Cl, Cl_in, Cl_out):
    """
    Computes membrane potential (Vm) in millivolts using the Goldman-Hodgkin-Katz equation.
    """
    V_T = (R * T) / F * 1000  # Convert to mV

    numerator = (P_K * K_out + P_Na * Na_out + P_Cl * Cl_in)
    denominator = (P_K * K_in + P_Na * Na_in + P_Cl * Cl_out)

    if denominator == 0:
        raise ValueError("Denominator is zero; check input values.")
    
    Vm = V_T * math.log(numerator / denominator)
    return Vm

# Constants
R = 8.314  # J/(mol·K)
T = 294.15  # Kelvin
F = 96485  # C/mol

# Ion permeabilities and concentrations
P_K = 2
K_in = 140
K_out_values = [5, 10, 15, 25, 35, 100, 150, 180, 200, 220]

P_Na = 0.1
Na_in = 15
Na_out = 145

P_Cl = 0.45
Cl_in = 10
Cl_out = 110

# Calculate and print Vm for each K_out value
print("K_out (mM) | Membrane Potential (mV)")
print("-" * 35)
for K_out in K_out_values:
    try:
        Vm = ghk(R, T, F, P_K, K_in, K_out, P_Na, Na_in, Na_out, P_Cl, Cl_in, Cl_out)
        print(f"{Vm:21.2f}")
    except ValueError as e:
        print(f"{K_out:10} | Error: {e}")