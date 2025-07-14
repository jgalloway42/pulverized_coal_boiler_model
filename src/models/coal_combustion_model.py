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
    lb_NO_per_lb_coal = lb_NO_per_MMBTU * (lb_HHV / 1e6)
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
    #  flame temperature in the 1800–2200 K range.
    Cp = {
        'CO2': 75.0,  # patched from 55.3
        'H2O': 70.0,  # patched from 51.0
        'SO2': 65.0,  # patched from 48.0
        'CO': 60.0,   # patched from 35.0
        'O2': 65.0,   # patched from 38.0
        'N2': 60.0,   # patched from 33.0
        'NO': 60.0    # patched from 36.0
    }
    HHV_J = HHV_BTU_per_lb * 1055.06
    total_Cp = sum(products_mol[sp] * Cp.get(sp, 0) for sp in products_mol)
    delta_T = HHV_J / total_Cp
    return T_ref_K + delta_T

def coal_combustion_from_mass_flow(ultimate, coal_lb_per_hr, air_scfh, HHV_btu_per_lb=12900, CO2_frac=0.9, NOx_eff=0.35):
    """
    Perform combustion analysis using mass flow rate of coal and volumetric air input.

    Parameters:
        ultimate (dict): Ultimate analysis of coal (% by mass: C, H, O, N, S, Moisture).
        coal_lb_per_hr (float): Coal feed rate in lb/hr.
        air_scfh (float): Combustion air in standard cubic feet per hour.
        HHV_btu_per_lb (float): Higher heating value of coal in BTU/lb.
        CO2_frac (float): Fraction of carbon forming CO2.
        NOx_eff (float): Fuel-nitrogen to NO conversion efficiency.

    Returns:
        dict: Combustion products and emissions per hour.
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

    # Air composition (mole fraction)
    air_O2_frac = 0.21
    air_N2_frac = 0.79
    air_mol_per_scf = 1.0 / 379.0  # mol/scf at STP

    # Convert ultimate analysis to mass fractions
    total = sum(ultimate.values())
    frac = {k: v / total for k, v in ultimate.items()}
    dry_frac = 1.0 - frac.get('Moisture', 0)

    coal_dry_lb_per_hr = coal_lb_per_hr * dry_frac

    # Moles of elements per lb of dry coal
    mol_C = frac['C'] / MW['C']
    mol_H = frac['H'] / MW['H']
    mol_O = frac['O'] / MW['O']
    mol_S = frac['S'] / MW['S']
    mol_N = frac['N'] / MW['N']

    total_mol_C = mol_C * coal_dry_lb_per_hr
    total_mol_H = mol_H * coal_dry_lb_per_hr
    total_mol_O = mol_O * coal_dry_lb_per_hr
    total_mol_S = mol_S * coal_dry_lb_per_hr
    total_mol_N = mol_N * coal_dry_lb_per_hr

    CO_frac = 1 - CO2_frac
    mol_CO2 = total_mol_C * CO2_frac
    mol_CO = total_mol_C * CO_frac

    mol_O2_required = mol_CO2 + mol_CO * 0.5 + total_mol_H / 4 + total_mol_S - total_mol_O / 2

    mol_air = air_scfh * air_mol_per_scf
    mol_O2_actual = mol_air * air_O2_frac
    mol_N2_air = mol_air * air_N2_frac
    mol_H2O_air = mol_air * (0.015 / 0.207)

    mol_NO = total_mol_N * NOx_eff
    mol_N2_total = mol_N2_air + (total_mol_N - mol_NO)
    mol_O2_excess = mol_O2_actual - mol_O2_required
    mol_H2O = total_mol_H / 2 + mol_H2O_air

    products_mol = {
        'CO2': mol_CO2,
        'CO': mol_CO,
        'H2O': mol_H2O,
        'SO2': total_mol_S,
        'N2': mol_N2_total,
        'O2': mol_O2_excess,
        'NO': mol_NO
    }

    T_flame_K = estimate_flame_temp_Cp_method(products_mol, HHV_btu_per_lb * dry_frac)

    O2_mole_frac = mol_O2_excess / sum(products_mol.values())
    lb_NO_thermal = estimate_thermal_NO(HHV_btu_per_lb * dry_frac, O2_mole_frac, T_flame_K) * coal_lb_per_hr

    def moles_to_lbs(mol, mw): return mol * mw / 453.592

    results = {
        'CO2 (LB)': moles_to_lbs(mol_CO2, MW['CO2']),
        'CO (LB)': moles_to_lbs(mol_CO, MW['CO']),
        'H2O (LB)': moles_to_lbs(mol_H2O, MW['H2O']),
        'SO2 (LB)': moles_to_lbs(total_mol_S, MW['SO2']),
        'NO (LB from fuel-N)': moles_to_lbs(mol_NO, MW['NO']),
        'NO (LB thermal)': lb_NO_thermal,
        'NO (LB total)': moles_to_lbs(mol_NO, MW['NO']) + lb_NO_thermal,
        'N2 (LB)': moles_to_lbs(mol_N2_total, MW['N2']),
        'O2 (excess %)': moles_to_lbs(mol_O2_excess, MW['O2']),
        'Heat Released (BTU/hr)': HHV_btu_per_lb * dry_frac * coal_lb_per_hr * (CO2_frac + 0.3 * CO_frac),
        'Estimated Flame Temp (F)': (T_flame_K - 273.15) * 9/5 + 32 # convert to Farrenheit
    }

    results['Total Flue Gas (lb/hr)'] = sum(v for k, v in results.items() if k not in ['Heat Released (BTU/hr)', 'Estimated Flame Temp (K)', 'Total Flue Gas (lb/hr)'])
    return results

# Example usage:
if __name__ == "__main__":
    # Sample ultimate analysis (mass percent)
    ultimate = {
        'C': 72.0,
        'H': 5.0,
        'O': 10.0,
        'N': 1.5,
        'S': 1.0,
        'Ash': 8.0,
        'Moisture': 2.5
    }

    # Input conditions
    coal_rate_lb_hr = 10000    # lb/hr of coal
    air_flow_scfh = 1800000    # standard cubic feet per hour of air

    # Run combustion model
    results = coal_combustion_from_mass_flow(
        ultimate,
        coal_lb_per_hr=coal_rate_lb_hr,
        air_scfh=air_flow_scfh,
        HHV_btu_per_lb=12900,
        CO2_frac=0.9,
        NOx_eff=0.35
    )

    # Output results
    print("\n=== Combustion Results (per hour) ===\n")
    for k, v in results.items():
        print(f"{k}: {v:,.2f}")
