import numpy as np

def estimate_thermal_NO(lb_HHV, O2_mole_frac=0.05, T_flame_K=2000):
    """
    Estimate thermal NO formation using an empirical Zeldovich-based model.

    Parameters:
        lb_HHV (float): Heat of combustion in BTU per lb of fuel.
        O2_mole_frac (float): Oxygen mole fraction in the flame region.
        T_flame_K (float): Flame temperature in Kelvin.

    Returns:
        float: Estimated thermal NO in lb per lb of coal.
    """
    B = 1.6e5
    C = 38000
    m = 0.5
    lb_NO_per_MMBTU = B * (O2_mole_frac ** m) * np.exp(-C / T_flame_K)
    lb_NO_per_lb_coal = lb_NO_per_MMBTU * (lb_HHV / 1_000_000)
    return lb_NO_per_lb_coal

def estimate_flame_temp_Cp_method(products_mol, HHV_BTU_per_lb, T_ref_K=298.15):
    """
    Estimate adiabatic flame temperature using constant Cp and energy balance.

    Parameters:
        products_mol (dict): Molar amounts of combustion products.
        HHV_BTU_per_lb (float): Higher heating value in BTU/lb.
        T_ref_K (float): Reference temperature in Kelvin (default: 298.15 K).

    Returns:
        float: Estimated flame temperature in Kelvin.
    """
    Cp = {
        'CO2': 55.3,
        'H2O': 51.0,
        'SO2': 48.0,
        'CO': 35.0,
        'O2': 38.0,
        'N2': 33.0,
        'NO': 36.0
    }
    HHV_J = HHV_BTU_per_lb * 1055.06
    total_Cp = sum(products_mol[sp] * Cp.get(sp, 0) for sp in products_mol)
    delta_T = HHV_J / total_Cp
    return T_ref_K + delta_T

def coal_combustion_mass_based_with_NOx(ultimate, HHV_btu_per_lb=12900, excess_air=0.2, CO2_frac=0.9, NOx_eff=0.35):
    """
    Perform mass-based combustion analysis of coal including NOx and flame temperature estimation.

    Parameters:
        ultimate (dict): Ultimate analysis of coal (% by mass: C, H, O, N, S, Moisture).
        HHV_btu_per_lb (float): Higher heating value of coal in BTU/lb.
        excess_air (float): Fractional excess air (e.g., 0.2 for 20%).
        CO2_frac (float): Fraction of carbon forming CO2.
        NOx_eff (float): Fuel-nitrogen to NO conversion efficiency (e.g., 0.35).

    Returns:
        dict: Mass of combustion products per lb of coal and estimated flame temperature.
    """
    MW = {
        'C': 12.01,
        'H': 1.008,
        'O': 16.00,
        'N': 14.01,
        'S': 32.07,
        'CO2': 44.01,
        'CO': 28.01,
        'H2O': 18.02,
        'SO2': 64.07,
        'O2': 32.00,
        'N2': 28.02,
        'NO': 30.01
    }

    # Convert ultimate analysis to mass fractions
    total = sum(ultimate.values())
    frac = {k: v / total for k, v in ultimate.items()}
    dry_frac = 1.0 - frac.get('Moisture', 0)

    # Moles of elements per lb of dry coal
    mol_C = frac['C'] / MW['C']
    mol_H = frac['H'] / MW['H']
    mol_O = frac['O'] / MW['O']
    mol_S = frac['S'] / MW['S']
    mol_N = frac['N'] / MW['N']

    # Carbon split into CO2 and CO
    CO_frac = 1 - CO2_frac
    mol_CO2 = mol_C * CO2_frac
    mol_CO = mol_C * CO_frac

    # Oxygen requirement and combustion air
    mol_O2_required = mol_CO2 + mol_CO * 0.5 + mol_H / 4 + mol_S - mol_O / 2
    mol_O2_actual = mol_O2_required * (1 + excess_air)
    mol_N2_air = mol_O2_actual * (79 / 21)
    mol_H2O_air = mol_O2_actual * (0.015 / 0.207)

    # NO formation from fuel nitrogen
    mol_NO = mol_N * NOx_eff
    mol_N2_total = mol_N2_air + (mol_N - mol_NO)
    mol_O2_excess = mol_O2_actual - mol_O2_required
    mol_H2O = mol_H / 2 + mol_H2O_air

    # Combustion products in mols
    products_mol = {
        'CO2': mol_CO2,
        'CO': mol_CO,
        'H2O': mol_H2O,
        'SO2': mol_S,
        'N2': mol_N2_total,
        'O2': mol_O2_excess,
        'NO': mol_NO
    }

    # Estimate flame temperature
    T_flame_K = estimate_flame_temp_Cp_method(products_mol, HHV_btu_per_lb * dry_frac)

    # Estimate thermal NO
    O2_mole_frac = mol_O2_excess / sum(products_mol.values())
    lb_NO_thermal = estimate_thermal_NO(HHV_btu_per_lb * dry_frac, O2_mole_frac, T_flame_K)

    # Moles to pounds converter
    def moles_to_lbs(mol, mw): return mol * mw / 453.592

    # Final results in lbs per lb of coal
    results = {
        'CO2': moles_to_lbs(mol_CO2, MW['CO2']),
        'CO': moles_to_lbs(mol_CO, MW['CO']),
        'H2O': moles_to_lbs(mol_H2O, MW['H2O']),
        'SO2': moles_to_lbs(mol_S, MW['SO2']),
        'NO (from fuel-N)': moles_to_lbs(mol_NO, MW['NO']),
        'NO (thermal)': lb_NO_thermal,
        'NO (total)': moles_to_lbs(mol_NO, MW['NO']) + lb_NO_thermal,
        'N2': moles_to_lbs(mol_N2_total, MW['N2']),
        'O2 (excess)': moles_to_lbs(mol_O2_excess, MW['O2']),
        'Heat Released (BTU)': HHV_btu_per_lb * dry_frac * (CO2_frac + 0.3 * CO_frac),
        'Estimated Flame Temp (K)': T_flame_K
    }

    # Total flue gas mass
    results['Total Flue Gas'] = sum(v for k, v in results.items() if k not in ['Heat Released (BTU)', 'Estimated Flame Temp (K)', 'Total Flue Gas'])
    return results

# Example usage:
if __name__ == "__main__":
    # Define ultimate analysis (mass percent)
    ultimate_analysis = {
        'C': 75.0,
        'H': 5.0,
        'O': 10.0,
        'N': 1.5,
        'S': 1.0,
        'Ash': 5.0,
        'Moisture': 2.5
    }

    # Run combustion model
    results = coal_combustion_mass_based_with_NOx(
        ultimate_analysis,
        HHV_btu_per_lb=12900,
        excess_air=0.2,
        CO2_frac=0.9,
        NOx_eff=0.35
    )

    # Print results
    print("\n=== Mass-Based Combustion Output with NOx (per lb of coal) ===\n")
    for k, v in results.items():
        print(f"{k}: {v:.4f} lb" if 'lb' in k or 'Gas' in k else f"{k}: {v:.2f}")
