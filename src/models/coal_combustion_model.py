import numpy as np

def estimate_thermal_NO(lb_HHV, O2_mole_frac=0.05, T_flame_K=2000):
    B = 1.6e5
    C = 38000
    m = 0.5
    lb_NO_per_MMBTU = B * (O2_mole_frac ** m) * np.exp(-C / T_flame_K)
    lb_NO_per_lb_coal = lb_NO_per_MMBTU * (lb_HHV / 1_000_000)
    return lb_NO_per_lb_coal

def estimate_flame_temp_Cp_method(products_mol, HHV_BTU_per_lb, T_ref_K=298.15):
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

    total = sum(ultimate.values())
    frac = {k: v / total for k, v in ultimate.items()}
    dry_frac = 1.0 - frac.get('Moisture', 0)

    mol_C = frac['C'] / MW['C']
    mol_H = frac['H'] / MW['H']
    mol_O = frac['O'] / MW['O']
    mol_S = frac['S'] / MW['S']
    mol_N = frac['N'] / MW['N']

    CO_frac = 1 - CO2_frac
    mol_CO2 = mol_C * CO2_frac
    mol_CO = mol_C * CO_frac

    mol_O2_required = mol_CO2 * 1 + mol_CO * 0.5 + mol_H / 4 + mol_S - mol_O / 2
    mol_O2_actual = mol_O2_required * (1 + excess_air)
    mol_N2_air = mol_O2_actual * (79 / 21)
    mol_H2O_air = mol_O2_actual * (0.015 / 0.207)

    mol_NO = mol_N * NOx_eff
    mol_N2_total = mol_N2_air + (mol_N - mol_NO)
    mol_O2_excess = mol_O2_actual - mol_O2_required
    mol_H2O = mol_H / 2 + mol_H2O_air

    products_mol = {
        'CO2': mol_CO2,
        'CO': mol_CO,
        'H2O': mol_H2O,
        'SO2': mol_S,
        'N2': mol_N2_total,
        'O2': mol_O2_excess,
        'NO': mol_NO
    }

    T_flame_K = estimate_flame_temp_Cp_method(products_mol, HHV_btu_per_lb * dry_frac)

    O2_mole_frac = mol_O2_excess / sum(products_mol.values())
    lb_NO_thermal = estimate_thermal_NO(HHV_btu_per_lb * dry_frac, O2_mole_frac, T_flame_K)

    def moles_to_lbs(mol, mw): return mol * mw / 453.592

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

    results['Total Flue Gas'] = sum(v for k, v in results.items() if k not in ['Heat Released (BTU)', 'Estimated Flame Temp (K)', 'Total Flue Gas'])
    return results
