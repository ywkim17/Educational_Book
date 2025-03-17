import math

def nernst_equation(R, T, F, z, ion_in, ion_out):
    """
    Computes the Nernst equilibrium potential (E_ion) in millivolts.

    Parameters:
        R       : float - Universal gas constant (J/(mol·K))
        T       : float - Temperature in Kelvin (K)
        F       : float - Faraday's constant (C/mol)
        z       : float - Charge of the ion (e.g., +1 for K⁺, Na⁺; -1 for Cl⁻)
        ion_in  : float - Intracellular ion concentration (mM)
        ion_out : float - Extracellular ion concentration (mM)
    
    Returns:
        E_ion (float): Equilibrium potential in millivolts (mV)
    """

    if ion_in <= 0 or ion_out <= 0:
        raise ValueError("Ion concentrations must be greater than zero.")

    # Calculate RT/zF in millivolts
    V_T = (R * T) / (z * F) * 1000  # Convert to mV

    E_ion = V_T * math.log(ion_out / ion_in)
    return E_ion

# User Inputs
R = 8.314
T = 310
F = 6485

# Example Ion Data
z = 1
ion_in = 200
ion_out = 5

# Compute Nernst Potential
E_ion = nernst_equation(R, T, F, z, ion_in, ion_out)

print(f"\nEquilibrium potential: {E_ion:.2f} mV")