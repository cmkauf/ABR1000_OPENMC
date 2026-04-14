"""
ABR-1000 OpenMC Depletion Model — Equilibrium TRU-Recycled Metal Core
One full irradiation cycle (328.5 EFPD) as a function of burnup.

References:
  [ANL]  ANL-AFCI-202, Argonne National Laboratory, September 2007.
  [NEA]  NEA-NSC-R(2015)9, Nuclear Energy Agency, February 2016.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DEPLETION DESIGN BASIS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Thermal power:         1000 MWt  (constant throughout cycle)
  Cycle length:          328.5 EFPD  (12 months × 90% capacity factor)
  HM inventory (BOEC):  ~13.2 MT  ([ANL] Table II.1-3)
  Specific power:        ~73 kW/kg HM
  Expected k swing:      ~2.2% Δk  (BOEC→EOEC, [ANL] Table II.1-3)
  Expected discharge BU: ~93 MWd/kg  (4-batch fuel management)

  Timesteps: 12 geometric steps (0.36 → 137.7 d), finer at start to
  capture short-lived fission product buildup (Xe-135, Sm-149).

  Integrator: CECMIntegrator (predictor-corrector, 2nd-order accurate).
  Substeps:   2 OpenMC transport solves per timestep.

DEPLETABLE MATERIALS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Only the 10 active-core sub-zone materials (5 axial × 2 radial) are
  marked depletable=True. Structural, coolant, shield, and reflector
  materials are not depleted (no fission, negligible transmutation).

CHAIN FILE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Set CHAIN_FILE at the top of this script to your local path. Common
  locations:
    ENDF/B-VIII.0 SFR chain: chain_endfb80_sfr.xml
    ENDF/B-VIII.0 full:      chain_endfb80.xml
  The SFR chain is preferred for fast-spectrum depletion; it includes
  optimised branching ratios for heavy actinides at fast energies.

POST-PROCESSING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  After the run, results are read from openmc.deplete.Results and
  printed / saved to CSV:
    - k-eff vs. time and burnup
    - Actinide inventory (U, Pu, Am, Cm, TRU, HM) vs. burnup
    - TRU consumption rate and conversion ratio

  The geometry and neutronics model are identical to v5 (all three
  bias-reduction tweaks applied).
"""

import math
import numpy as np
import openmc
import openmc.deplete

# ======================================================================
# USER CONFIGURATION
# ======================================================================
# Path to your OpenMC depletion chain file.
# Use the SFR-optimised chain if available; fall back to the full chain.
CHAIN_FILE = 'chain_endfb80_sfr.xml'   # <-- set this to your path

# Thermal power [W]
POWER_W = 1000.0e6

# Cycle length [effective full-power days] = 12 months × 90% capacity factor
CYCLE_EFPD = 328.5

# Depletion particles per batch (fewer than k-eff run is acceptable
# since depletion only needs accurate one-group cross sections, not k)
DEPL_PARTICLES = 20000

# k-eff transport solve settings (used at each depletion timestep)
DEPL_BATCHES  = 60
DEPL_INACTIVE = 20

# ======================================================================
# TEMPERATURES  [NEA] Table 2.15
# ======================================================================
T_cool   = 432.5 + 273.15
T_struct = T_cool
T_fuel   = 534.0 + 273.15

# ======================================================================
# GEOMETRY CONSTANTS
# ======================================================================
pitch = 16.142
A_hex = (math.sqrt(3) / 2.0) * pitch**2

def r_eff(n):
    return math.sqrt(n * A_hex / math.pi)

R_inner  = r_eff(85)
R_outer  = r_eff(199)
R_refl   = r_eff(313)
R_shield = r_eff(379)

H_active  = 85.82
H_lower_r = 125.16
H_lower_s = 35.76
H_bond_na = 20.06
H_gas_pl  = 101.01
H_upper_s = 112.39

z_fuel_bot = -H_active / 2.0
z_fuel_top = +H_active / 2.0
z_lr_bot   = z_fuel_bot - H_lower_r
z_ls_bot   = z_lr_bot   - H_lower_s
z_bond_top = z_fuel_top + H_bond_na
z_pl_top   = z_bond_top + H_gas_pl
z_us_top   = z_pl_top   + H_upper_s

N_AX = 5
dz   = H_active / N_AX
ax_z = [z_fuel_bot + i * dz for i in range(N_AX + 1)]

# ======================================================================
# VOLUME FRACTIONS  [NEA] Table 2.20
# ======================================================================
VF_FUEL = 0.3900
VF_NA   = 0.3534
VF_HT9  = 0.2566

# ======================================================================
# PURE-MATERIAL NUMBER DENSITIES  [NEA] Table 2.21  (a/b-cm)
# ======================================================================
ND_NA      = 2.2272e-02
ND_HT9_Fe  = 6.9715e-02
ND_HT9_Ni  = 4.2984e-04
ND_HT9_Cr  = 1.0366e-02
ND_HT9_Mn  = 4.5921e-04
ND_HT9_Mo  = 4.9007e-04
ND_LS_Na   = 1.5591e-02
ND_LS_Fe   = 1.5878e-02
ND_LS_Ni   = 3.2604e-03
ND_LS_Cr   = 3.2355e-03
ND_LS_Mn   = 5.0846e-04
ND_LS_Mo   = 4.3524e-04

# ======================================================================
# FUEL-PIN NUMBER DENSITIES  [NEA] Tables 2.22 / 2.23
# ======================================================================
_mean = lambda lst: sum(lst) / len(lst)

_I_zones = {
    'U234':     [1.1369e-6, 1.0856e-6, 1.0727e-6, 1.1028e-6, 1.1759e-6],
    'U235':     [3.0421e-5, 2.9338e-5, 2.8961e-5, 3.0070e-5, 3.2571e-5],
    'U236':     [2.4896e-6, 2.5117e-6, 2.5536e-6, 2.3779e-6, 2.0226e-6],
    'U238':     [1.9613e-2, 1.9474e-2, 1.9433e-2, 1.9550e-2, 1.9801e-2],
    'Np237':    [4.6686e-5, 4.6962e-5, 4.6782e-5, 4.7603e-5, 4.8895e-5],
    'Pu238':    [1.1695e-4, 1.1284e-4, 1.1196e-4, 1.1370e-4, 1.1829e-4],
    'Pu239':    [2.2076e-3, 2.1814e-3, 2.1754e-3, 2.1813e-3, 2.2011e-3],
    'Pu240':    [1.3244e-3, 1.2955e-3, 1.2902e-3, 1.2986e-3, 1.3248e-3],
    'Pu241':    [1.9375e-4, 1.8610e-4, 1.8518e-4, 1.8537e-4, 1.8845e-4],
    'Pu242':    [2.9277e-4, 2.8911e-4, 2.8818e-4, 2.9038e-4, 2.9569e-4],
    'Am241':    [1.0791e-4, 1.0465e-4, 1.0353e-4, 1.0686e-4, 1.1421e-4],
    'Am242_m1': [9.2989e-6, 9.0848e-6, 9.0224e-6, 9.1756e-6, 9.4890e-6],
    'Am243':    [1.0017e-4, 9.8324e-5, 9.7993e-5, 9.8630e-5, 1.0032e-4],
    'Cm242':    [5.6250e-6, 5.8208e-6, 5.9476e-6, 5.4901e-6, 4.5416e-6],
    'Cm243':    [5.4321e-7, 5.0246e-7, 5.0136e-7, 4.8876e-7, 4.8480e-7],
    'Cm244':    [6.7240e-5, 6.5722e-5, 6.5622e-5, 6.5349e-5, 6.5394e-5],
    'Cm245':    [1.7397e-5, 1.6743e-5, 1.6663e-5, 1.6696e-5, 1.7026e-5],
    'Cm246':    [9.2285e-6, 9.1426e-6, 9.1307e-6, 9.1364e-6, 9.1805e-6],
    'Zr':       [7.2802e-3, 7.2802e-3, 7.2802e-3, 7.2802e-3, 7.2802e-3],
    'Mo_fp':    [9.2873e-4, 1.1464e-3, 1.2031e-3, 1.0625e-3, 7.4065e-4],
}

_O_zones = {
    'U234':     [1.6317e-6, 1.5766e-6, 1.5638e-6, 1.5894e-6, 1.6552e-6],
    'U235':     [3.0822e-5, 2.9870e-5, 2.9561e-5, 3.0391e-5, 3.2250e-5],
    'U236':     [1.7881e-6, 1.8534e-6, 1.8941e-6, 1.7528e-6, 1.4710e-6],
    'U238':     [1.8244e-2, 1.8144e-2, 1.8115e-2, 1.8191e-2, 1.8359e-2],
    'Np237':    [9.8244e-5, 9.7300e-5, 9.6775e-5, 9.8481e-5, 1.0175e-4],
    'Pu238':    [1.6436e-4, 1.6026e-4, 1.5949e-4, 1.6063e-4, 1.6416e-4],
    'Pu239':    [2.8147e-3, 2.7664e-3, 2.7538e-3, 2.7786e-3, 2.8416e-3],
    'Pu240':    [1.7467e-3, 1.7191e-3, 1.7135e-3, 1.7231e-3, 1.7508e-3],
    'Pu241':    [2.8976e-4, 2.8138e-4, 2.8012e-4, 2.8135e-4, 2.8697e-4],
    'Pu242':    [4.0754e-4, 4.0412e-4, 4.0321e-4, 4.0530e-4, 4.1028e-4],
    'Am241':    [1.8607e-4, 1.8127e-4, 1.7970e-4, 1.8397e-4, 1.9339e-4],
    'Am242_m1': [1.2185e-5, 1.2045e-5, 1.2021e-5, 1.2039e-5, 1.2064e-5],
    'Am243':    [1.3234e-4, 1.3019e-4, 1.2985e-4, 1.3036e-4, 1.3206e-4],
    'Cm242':    [6.4688e-6, 6.8630e-6, 7.0553e-6, 6.4446e-6, 5.1976e-6],
    'Cm243':    [6.3471e-7, 6.0893e-7, 6.0901e-7, 5.9753e-7, 5.9372e-7],
    'Cm244':    [8.0107e-5, 7.8889e-5, 7.8847e-5, 7.8479e-5, 7.8359e-5],
    'Cm245':    [2.0200e-5, 1.9678e-5, 1.9613e-5, 1.9635e-5, 1.9913e-5],
    'Cm246':    [1.0443e-5, 1.0371e-5, 1.0361e-5, 1.0367e-5, 1.0410e-5],
    'Zr':       [7.2802e-3, 7.2802e-3, 7.2802e-3, 7.2802e-3, 7.2802e-3],
    'Mo_fp':    [8.1524e-4, 1.0174e-3, 1.0697e-3, 9.4870e-4, 6.6172e-4],
}

# ======================================================================
# MATERIAL BUILDERS
# ======================================================================
_ACTINIDES = [
    'U234','U235','U236','U238',
    'Np237',
    'Pu238','Pu239','Pu240','Pu241','Pu242',
    'Am241','Am242_m1','Am243',
    'Cm242','Cm243','Cm244','Cm245','Cm246',
]

def make_active_zone(zone_data, iz, label, T, depletable=False):
    m = openmc.Material(name=f'active_{label}_z{iz}')
    for nuc in _ACTINIDES:
        m.add_nuclide(nuc, zone_data[nuc][iz] * VF_FUEL)
    m.add_element('Zr', zone_data['Zr'][iz] * VF_FUEL)
    mo_total = zone_data['Mo_fp'][iz] * VF_FUEL + ND_HT9_Mo * VF_HT9
    m.add_element('Mo', mo_total)
    m.add_nuclide('Na23', ND_NA * VF_NA)
    m.add_element('Fe',   ND_HT9_Fe * VF_HT9)
    m.add_element('Ni',   ND_HT9_Ni * VF_HT9)
    m.add_element('Cr',   ND_HT9_Cr * VF_HT9)
    m.add_nuclide('Mn55', ND_HT9_Mn * VF_HT9)
    total = (
        sum(zone_data[nuc][iz] * VF_FUEL for nuc in _ACTINIDES)
        + zone_data['Zr'][iz] * VF_FUEL
        + mo_total
        + ND_NA * VF_NA
        + (ND_HT9_Fe + ND_HT9_Ni + ND_HT9_Cr + ND_HT9_Mn) * VF_HT9
    )
    m.set_density('atom/b-cm', total)
    m.temperature = T
    m.depletable = depletable
    return m


def make_na_ht9(vf_na, vf_ht9, name, T):
    assert abs(vf_na + vf_ht9 - 1.0) < 1e-5, f"{name}: {vf_na+vf_ht9}"
    m = openmc.Material(name=name)
    m.add_nuclide('Na23',  ND_NA     * vf_na)
    m.add_element('Fe',    ND_HT9_Fe * vf_ht9)
    m.add_element('Ni',    ND_HT9_Ni * vf_ht9)
    m.add_element('Cr',    ND_HT9_Cr * vf_ht9)
    m.add_nuclide('Mn55',  ND_HT9_Mn * vf_ht9)
    m.add_element('Mo',    ND_HT9_Mo * vf_ht9)
    total = (ND_NA * vf_na
             + (ND_HT9_Fe + ND_HT9_Ni + ND_HT9_Cr + ND_HT9_Mn + ND_HT9_Mo) * vf_ht9)
    m.set_density('atom/b-cm', total)
    m.temperature = T
    return m

# ======================================================================
# NON-FUEL MATERIALS  (depletable=False, default)
# ======================================================================
lower_refl_mat   = make_na_ht9(0.3534, 0.6466, 'lower_reflector',  T_struct)
upper_struct_mat = make_na_ht9(0.3534, 0.6466, 'upper_structure',   T_struct)
bond_na_mat      = make_na_ht9(0.7434, 0.2566, 'displaced_bond_na', T_struct)
plenum_mat       = make_na_ht9(0.5793, 0.4207, 'gas_plenum',        T_struct)

lower_struct_mat = openmc.Material(name='lower_structure')
lower_struct_mat.add_nuclide('Na23', ND_LS_Na)
lower_struct_mat.add_element('Fe',   ND_LS_Fe)
lower_struct_mat.add_element('Ni',   ND_LS_Ni)
lower_struct_mat.add_element('Cr',   ND_LS_Cr)
lower_struct_mat.add_nuclide('Mn55', ND_LS_Mn)
lower_struct_mat.add_element('Mo',   ND_LS_Mo)
lower_struct_mat.set_density('atom/b-cm',
    ND_LS_Na + ND_LS_Fe + ND_LS_Ni + ND_LS_Cr + ND_LS_Mn + ND_LS_Mo)
lower_struct_mat.temperature = T_struct

radial_refl_mat = make_na_ht9(0.1550, 0.8450, 'radial_reflector', T_struct)

_VF_sh_Na  = 0.1710; _VF_sh_HT9 = 0.2968; _VF_sh_B4C = 1.0 - _VF_sh_Na - _VF_sh_HT9
radial_shield = openmc.Material(name='radial_shield')
radial_shield.add_nuclide('Na23',  ND_NA      * _VF_sh_Na)
radial_shield.add_element('Fe',    ND_HT9_Fe  * _VF_sh_HT9)
radial_shield.add_element('Ni',    ND_HT9_Ni  * _VF_sh_HT9)
radial_shield.add_element('Cr',    ND_HT9_Cr  * _VF_sh_HT9)
radial_shield.add_nuclide('Mn55',  ND_HT9_Mn  * _VF_sh_HT9)
radial_shield.add_element('Mo',    ND_HT9_Mo  * _VF_sh_HT9)
radial_shield.add_nuclide('B10',   1.5018e-02 * _VF_sh_B4C)
radial_shield.add_nuclide('B11',   6.3609e-02 * _VF_sh_B4C)
radial_shield.add_nuclide('C12',   1.9657e-02 * _VF_sh_B4C)
radial_shield.set_density('atom/b-cm',
    ND_NA * _VF_sh_Na
    + (ND_HT9_Fe+ND_HT9_Ni+ND_HT9_Cr+ND_HT9_Mn+ND_HT9_Mo) * _VF_sh_HT9
    + (1.5018e-2+6.3609e-2+1.9657e-2) * _VF_sh_B4C)
radial_shield.temperature = T_struct

_VF_ca_Na=0.2883; _VF_ca_HT9=0.2077; _VF_ca_B4C=1.0-_VF_ca_Na-_VF_ca_HT9
ctrl_absorber = openmc.Material(name='ctrl_absorber')
ctrl_absorber.add_nuclide('Na23',  ND_NA      * _VF_ca_Na)
ctrl_absorber.add_element('Fe',    ND_HT9_Fe  * _VF_ca_HT9)
ctrl_absorber.add_element('Ni',    ND_HT9_Ni  * _VF_ca_HT9)
ctrl_absorber.add_element('Cr',    ND_HT9_Cr  * _VF_ca_HT9)
ctrl_absorber.add_nuclide('Mn55',  ND_HT9_Mn  * _VF_ca_HT9)
ctrl_absorber.add_element('Mo',    ND_HT9_Mo  * _VF_ca_HT9)
ctrl_absorber.add_nuclide('B10',   5.3642e-02 * _VF_ca_B4C)
ctrl_absorber.add_nuclide('B11',   2.8884e-02 * _VF_ca_B4C)
ctrl_absorber.add_nuclide('C12',   2.0632e-02 * _VF_ca_B4C)
ctrl_absorber.set_density('atom/b-cm',
    ND_NA * _VF_ca_Na
    + (ND_HT9_Fe+ND_HT9_Ni+ND_HT9_Cr+ND_HT9_Mn+ND_HT9_Mo) * _VF_ca_HT9
    + (5.3642e-2+2.8884e-2+2.0632e-2) * _VF_ca_B4C)
ctrl_absorber.temperature = T_struct

ctrl_empty = make_na_ht9(0.9074, 0.0926, 'ctrl_empty_duct', T_struct)

# ======================================================================
# CONTROL FRACTIONS AND ROD INSERTION
# ======================================================================
f_ci = 7.0  / 85.0
f_co = 12.0 / 114.0
F_INSERT = 0.50  # BOEC partial rod insertion fraction

def make_active_blended(fuel_mat, f_ctrl, label, T):
    """Three-way blend: fuel + partial absorber + empty duct."""
    f_fuel  = 1.0 - f_ctrl
    f_abs   = f_ctrl * F_INSERT
    f_empty = f_ctrl * (1.0 - F_INSERT)
    m = openmc.Material.mix_materials(
        [fuel_mat,  ctrl_absorber, ctrl_empty],
        [f_fuel,    f_abs,         f_empty],
        'vo'
    )
    m.name  = label
    m.temperature = T
    # Preserve depletable flag from the fuel material
    m.depletable = fuel_mat.depletable
    return m

# ======================================================================
# ACTIVE-CORE MATERIALS — marked depletable=True
# 10 materials: 5 axial zones × 2 radial zones
# ======================================================================
inner_fuel_mats = [
    make_active_zone(_I_zones, iz, 'inner', T_fuel, depletable=True)
    for iz in range(N_AX)
]
outer_fuel_mats = [
    make_active_zone(_O_zones, iz, 'outer', T_fuel, depletable=True)
    for iz in range(N_AX)
]

inner_ac_mats = [
    make_active_blended(inner_fuel_mats[iz], f_ci, f'inner_ac_z{iz}', T_fuel)
    for iz in range(N_AX)
]
outer_ac_mats = [
    make_active_blended(outer_fuel_mats[iz], f_co, f'outer_ac_z{iz}', T_fuel)
    for iz in range(N_AX)
]

# ======================================================================
# NON-ACTIVE AXIAL ZONE STRUCTURAL BLENDS
# ======================================================================
def blend_struct(base_mat, f_ctrl, name):
    m = openmc.Material.mix_materials(
        [base_mat, ctrl_empty], [1.0-f_ctrl, f_ctrl], 'vo'
    )
    m.name = name
    m.temperature = T_struct
    return m

inner_ls = blend_struct(lower_struct_mat, f_ci, 'inner_ls')
outer_ls = blend_struct(lower_struct_mat, f_co, 'outer_ls')
inner_lr = blend_struct(lower_refl_mat,   f_ci, 'inner_lr')
outer_lr = blend_struct(lower_refl_mat,   f_co, 'outer_lr')
inner_bn = blend_struct(bond_na_mat,      f_ci, 'inner_bn')
outer_bn = blend_struct(bond_na_mat,      f_co, 'outer_bn')
inner_gp = blend_struct(plenum_mat,       f_ci, 'inner_gp')
outer_gp = blend_struct(plenum_mat,       f_co, 'outer_gp')
inner_us = blend_struct(upper_struct_mat, f_ci, 'inner_us')
outer_us = blend_struct(upper_struct_mat, f_co, 'outer_us')

# ======================================================================
# MATERIALS COLLECTION
# ======================================================================
materials = openmc.Materials(
    inner_ac_mats + outer_ac_mats
    + [inner_ls, outer_ls, inner_lr, outer_lr,
       inner_bn, outer_bn, inner_gp, outer_gp,
       inner_us, outer_us,
       radial_refl_mat, radial_shield]
)

# ======================================================================
# GEOMETRY
# ======================================================================
cyl_inner  = openmc.ZCylinder(r=R_inner)
cyl_outer  = openmc.ZCylinder(r=R_outer)
cyl_refl   = openmc.ZCylinder(r=R_refl)
cyl_shield = openmc.ZCylinder(r=R_shield, boundary_type='vacuum')

zp_ls_bot = openmc.ZPlane(z0=z_ls_bot,   boundary_type='vacuum')
zp_lr_bot = openmc.ZPlane(z0=z_lr_bot)
zp_fu_bot = openmc.ZPlane(z0=ax_z[0])
zp_fu_top = openmc.ZPlane(z0=ax_z[N_AX])
zp_bn_top = openmc.ZPlane(z0=z_bond_top)
zp_pl_top = openmc.ZPlane(z0=z_pl_top)
zp_us_top = openmc.ZPlane(z0=z_us_top,   boundary_type='vacuum')

zp_ax = [openmc.ZPlane(z0=ax_z[i]) for i in range(1, N_AX)]

def make_zone_cells(rad_reg, label, mat_ls, mat_lr, ac_mats,
                    mat_bn, mat_gp, mat_us):
    cells = []
    cells.append(openmc.Cell(name=f'{label}_ls', fill=mat_ls,
        region=rad_reg & +zp_ls_bot & -zp_lr_bot))
    cells.append(openmc.Cell(name=f'{label}_lr', fill=mat_lr,
        region=rad_reg & +zp_lr_bot & -zp_fu_bot))
    cells.append(openmc.Cell(name=f'{label}_ac0', fill=ac_mats[0],
        region=rad_reg & +zp_fu_bot & -zp_ax[0]))
    for i in range(1, N_AX - 1):
        cells.append(openmc.Cell(name=f'{label}_ac{i}', fill=ac_mats[i],
            region=rad_reg & +zp_ax[i-1] & -zp_ax[i]))
    cells.append(openmc.Cell(name=f'{label}_ac4', fill=ac_mats[4],
        region=rad_reg & +zp_ax[N_AX-2] & -zp_fu_top))
    cells.append(openmc.Cell(name=f'{label}_bn', fill=mat_bn,
        region=rad_reg & +zp_fu_top & -zp_bn_top))
    cells.append(openmc.Cell(name=f'{label}_gp', fill=mat_gp,
        region=rad_reg & +zp_bn_top & -zp_pl_top))
    cells.append(openmc.Cell(name=f'{label}_us', fill=mat_us,
        region=rad_reg & +zp_pl_top & -zp_us_top))
    return cells

inner_cells = make_zone_cells(-cyl_inner, 'inner',
    inner_ls, inner_lr, inner_ac_mats, inner_bn, inner_gp, inner_us)
outer_cells = make_zone_cells(+cyl_inner & -cyl_outer, 'outer',
    outer_ls, outer_lr, outer_ac_mats, outer_bn, outer_gp, outer_us)

cell_refl   = openmc.Cell(name='radial_refl',   fill=radial_refl_mat,
    region=+cyl_outer & -cyl_refl  & +zp_ls_bot & -zp_us_top)
cell_shield = openmc.Cell(name='radial_shield',  fill=radial_shield,
    region=+cyl_refl  & -cyl_shield & +zp_ls_bot & -zp_us_top)

universe = openmc.Universe(
    cells=inner_cells + outer_cells + [cell_refl, cell_shield])
geometry = openmc.Geometry(universe)

# ======================================================================
# SETTINGS  (used at every depletion transport step)
# ======================================================================
settings = openmc.Settings()
settings.run_mode  = 'eigenvalue'
settings.batches   = DEPL_BATCHES
settings.inactive  = DEPL_INACTIVE
settings.particles = DEPL_PARTICLES

settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-R_outer, -R_outer, z_fuel_bot],
        [ R_outer,  R_outer, z_fuel_top],
        only_fissionable=True
    )
)

entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left  = [-R_outer, -R_outer, z_fuel_bot]
entropy_mesh.upper_right = [ R_outer,  R_outer, z_fuel_top]
entropy_mesh.dimension   = [20, 20, 20]
settings.entropy_mesh = entropy_mesh
settings.temperature  = {'default': T_struct, 'method': 'interpolation'}

# ======================================================================
# TALLIES
# (Depletion uses its own internal reaction-rate tallies automatically;
#  these are additional output tallies for physics analysis.)
# ======================================================================
inner_ac_cells = inner_cells[2:7]
outer_ac_cells = outer_cells[2:7]

fuel_filt  = openmc.CellFilter(inner_ac_cells + outer_ac_cells)
inner_filt = openmc.CellFilter(inner_ac_cells)
outer_filt = openmc.CellFilter(outer_ac_cells)
leak_filt  = openmc.SurfaceFilter([cyl_shield, zp_ls_bot, zp_us_top])
neut_filt  = openmc.ParticleFilter(['neutron'])

grpstr33 = [
    1.00000E-05, 4.17458E-01, 5.31578E-01, 3.92786E+00, 8.31528E+00,
    1.37096E+01, 2.26033E+01, 3.72665E+01, 6.14421E+01, 1.01301E+02,
    1.67017E+02, 2.75364E+02, 4.53999E+02, 7.48518E+02, 1.23410E+03,
    2.03468E+03, 3.35462E+03, 5.53084E+03, 9.11881E+03, 1.50344E+04,
    2.47875E+04, 4.08677E+04, 6.73794E+04, 1.11090E+05, 1.83156E+05,
    3.01974E+05, 4.97871E+05, 8.20850E+05, 1.35335E+06, 2.23130E+06,
    3.67879E+06, 6.06531E+06, 1.00000E+07, 2.00000E+07
]
E_filt = openmc.EnergyFilter(grpstr33)

tallies = openmc.Tallies()

t = openmc.Tally(name='total_power')
t.filters = [fuel_filt]
t.scores  = ['kappa-fission', 'fission', 'nu-fission']
tallies.append(t)

t = openmc.Tally(name='inner_power')
t.filters = [inner_filt]; t.scores = ['kappa-fission', 'fission']
tallies.append(t)

t = openmc.Tally(name='outer_power')
t.filters = [outer_filt]; t.scores = ['kappa-fission', 'fission']
tallies.append(t)

t = openmc.Tally(name='system_leakage')
t.filters = [neut_filt, leak_filt]; t.scores = ['current']
tallies.append(t)

t = openmc.Tally(name='flux_33g')
t.filters = [fuel_filt, E_filt]; t.scores = ['flux', 'nu-fission']
tallies.append(t)

t = openmc.Tally(name='axial_power_inner')
t.filters = [openmc.CellFilter(inner_ac_cells)]; t.scores = ['kappa-fission']
tallies.append(t)

t = openmc.Tally(name='axial_power_outer')
t.filters = [openmc.CellFilter(outer_ac_cells)]; t.scores = ['kappa-fission']
tallies.append(t)

# ======================================================================
# ASSEMBLE MODEL
# ======================================================================
model = openmc.Model(geometry=geometry, materials=materials,
                     settings=settings, tallies=tallies)

# ======================================================================
# DEPLETION TIMESTEP SCHEDULE
# 12 geometric steps spanning 328.5 EFPD (0.36 d → 137.7 d).
# Geometric spacing captures both the short-lived FP transient (first
# ~10 days) and the slow actinide burnup trend (last ~200 days).
# ======================================================================
# Reference HM mass for burnup labelling [kg]
# Computed from benchmark number densities and cylindrical geometry
_Na_av = 6.022e23
_A_hm  = {  # approximate atomic mass per actinide
    'U234':234,'U235':235,'U236':236,'U238':238,'Np237':237,
    'Pu238':238,'Pu239':239,'Pu240':240,'Pu241':241,'Pu242':242,
    'Am241':241,'Am242_m1':242,'Am243':243,
    'Cm242':242,'Cm243':243,'Cm244':244,'Cm245':245,'Cm246':246,
}
def _hm_mass_kg(zone_data, V_zone):
    m = 0.0
    for nuc in _ACTINIDES:
        nd = _mean(zone_data[nuc]) if isinstance(zone_data[nuc], list) else zone_data[nuc]
        m += nd * VF_FUEL * 1e24 * _A_hm[nuc] / _Na_av * V_zone
    return m / 1e3  # g -> kg

V_inner = math.pi * R_inner**2 * H_active
V_outer = math.pi * (R_outer**2 - R_inner**2) * H_active
M_HM_kg = _hm_mass_kg(_I_zones, V_inner) + _hm_mass_kg(_O_zones, V_outer)
M_HM_MT = M_HM_kg / 1e6
print(f"Computed HM mass: {M_HM_MT:.2f} MT  (ANL reference: 13.2 MT)")

# Geometric timestep list [days]
_t_pts = np.geomspace(0.5, CYCLE_EFPD, 13)   # 13 points = 12 intervals
timesteps_d = list(np.diff(_t_pts))
# Adjust last step so the total is exactly CYCLE_EFPD
timesteps_d[-1] += CYCLE_EFPD - sum(timesteps_d)
timesteps_d = [float(dt) for dt in timesteps_d]

# Cumulative burnup at each step endpoint for reporting
cum_EFPD = np.cumsum([0.0] + timesteps_d)
cum_BU   = POWER_W / 1e6 * cum_EFPD / (M_HM_kg / 1e3)  # MWd/kg

print("\nDepletion timestep schedule:")
print(f"  {'Step':>4}  {'dt (d)':>8}  {'cum (d)':>8}  {'BU (MWd/kg)':>12}")
for i, dt in enumerate(timesteps_d):
    print(f"  {i+1:4d}  {dt:8.2f}  {cum_EFPD[i+1]:8.1f}  {cum_BU[i+1]:12.1f}")
print(f"\n  Total: {sum(timesteps_d):.2f} EFPD  "
      f"Final BU: {cum_BU[-1]:.1f} MWd/kg")
print(f"  Expected discharge (×4 batches): {cum_BU[-1]*4:.0f} MWd/kg  "
      f"(ANL: 93 MWd/kg)")

# ======================================================================
# RUN DEPLETION
# ======================================================================
print("\n" + "="*65)
print("ABR-1000  Depletion Run  (CECM integrator, 12 timesteps)")
print("="*65)

integrator = openmc.deplete.CECMIntegrator(
    model,
    CHAIN_FILE,
    timesteps=timesteps_d,
    timestep_units='d',
    power=POWER_W,
)
integrator.integrate()

# ======================================================================
# POST-PROCESSING
# ======================================================================
print("\n" + "="*65)
print("POST-PROCESSING RESULTS")
print("="*65)

results = openmc.deplete.Results('depletion_results.h5')
n_steps = len(results)

# Nuclides to track for inventory
_track_u   = ['U235', 'U238']
_track_pu  = ['Pu239', 'Pu240', 'Pu241', 'Pu242']
_track_am  = ['Am241', 'Am243']
_track_cm  = ['Cm244', 'Cm245']
_track_all = _ACTINIDES  # full actinide list

# Collect all depletable material IDs
depl_mat_ids = [m.id for m in (inner_ac_mats + outer_ac_mats)]

def get_total_mass_kg(results, step_idx, nuclide):
    """Sum nuclide mass [kg] across all depletable materials at a given step."""
    total = 0.0
    _, mats = results[step_idx]
    for mat in mats:
        if mat.id in depl_mat_ids:
            try:
                total += mat.get_mass(nuclide)  # grams
            except Exception:
                pass
    return total / 1e3  # g -> kg

# Build results table
rows = []
for si in range(n_steps):
    time_d = results.get_times('d')[si]
    keff, keff_unc = results.get_keff()[si]

    # Actinide masses [kg]
    m_U235  = get_total_mass_kg(results, si, 'U235')
    m_U238  = get_total_mass_kg(results, si, 'U238')
    m_Pu239 = get_total_mass_kg(results, si, 'Pu239')
    m_Pu240 = get_total_mass_kg(results, si, 'Pu240')
    m_Pu241 = get_total_mass_kg(results, si, 'Pu241')
    m_Pu242 = get_total_mass_kg(results, si, 'Pu242')
    m_Am241 = get_total_mass_kg(results, si, 'Am241')
    m_Am243 = get_total_mass_kg(results, si, 'Am243')
    m_Cm244 = get_total_mass_kg(results, si, 'Cm244')
    m_Cm245 = get_total_mass_kg(results, si, 'Cm245')
    m_Np237 = get_total_mass_kg(results, si, 'Np237')

    # TRU = Np + all Pu + Am + Cm
    m_TRU = (m_Np237
             + m_Pu239 + m_Pu240 + m_Pu241 + m_Pu242
             + m_Am241 + m_Am243
             + m_Cm244 + m_Cm245)
    # Total HM ≈ U + TRU (neglect minor isotopes)
    m_HM  = m_U235 + m_U238 + m_TRU

    # Burnup [MWd/kgHM] — referenced to initial HM mass
    burnup = POWER_W / 1e6 * time_d / (M_HM_kg / 1e3)

    rows.append({
        'step': si,
        'time_d': time_d,
        'burnup': burnup,
        'keff': keff,
        'keff_unc': keff_unc,
        'U235_kg': m_U235,
        'U238_kg': m_U238,
        'Pu239_kg': m_Pu239,
        'Pu240_kg': m_Pu240,
        'Pu241_kg': m_Pu241,
        'Pu242_kg': m_Pu242,
        'Am241_kg': m_Am241,
        'Am243_kg': m_Am243,
        'Cm244_kg': m_Cm244,
        'Cm245_kg': m_Cm245,
        'Np237_kg': m_Np237,
        'TRU_kg': m_TRU,
        'HM_kg': m_HM,
    })

# Print k-eff vs. burnup table
print("\n  k-eff vs. Burnup:")
print(f"  {'Step':>4}  {'Time(d)':>8}  {'BU(MWd/kg)':>11}  "
      f"{'k-eff':>7}  {'±1σ':>6}")
print("  " + "-"*50)
for r in rows:
    print(f"  {r['step']:4d}  {r['time_d']:8.1f}  {r['burnup']:11.1f}  "
          f"  {r['keff']:.5f}  {r['keff_unc']:.5f}")

# k-eff reactivity swing
k_BOEC = rows[0]['keff']
k_EOEC = rows[-1]['keff']
dk_pct = (k_BOEC - k_EOEC) / (k_BOEC * k_EOEC) * 100  # Δρ in %Δk/k
print(f"\n  Burnup reactivity loss: {dk_pct:.2f}% Δk  "
      f"(ANL reference: 2.2% Δk)")
print(f"  k_BOEC = {k_BOEC:.5f}")
print(f"  k_EOEC = {k_EOEC:.5f}")

# Actinide inventory change
print("\n  Actinide inventory change (BOEC → EOEC):")
r0 = rows[0]; rn = rows[-1]
for nuc in ['U235','U238','Pu239','Pu240','Pu241','Pu242','Am241','Am243',
            'Cm244','Cm245','TRU','HM']:
    k = f'{nuc}_kg'
    m0 = r0[k]; mf = rn[k]
    dm = mf - m0
    print(f"    {nuc:<8}: {m0:8.1f} kg  ->  {mf:8.1f} kg  "
          f"(Δ = {dm:+7.1f} kg)")

# TRU conversion ratio
# CR_TRU = (TRU produced) / (TRU destroyed)
tru_destroyed = r0['TRU_kg'] - rn['TRU_kg']
# Pu239 produced from U238 capture (estimate: Δ(all Pu+Am+Cm beyond initial))
# Simplified: CR = 1 - |ΔTRU| / (TRU_initial * burnup_fraction)
print(f"\n  TRU destroyed: {tru_destroyed:.1f} kg/cycle  "
      f"(ANL reference: ~82 kg/year = 82 kg/cycle)")
tru_rate_per_year = tru_destroyed * 365.25 / CYCLE_EFPD
print(f"  TRU consumption rate: {tru_rate_per_year:.1f} kg/year")

# Save to CSV
import csv
csv_file = 'ABR1000_depletion_results.csv'
with open(csv_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(f"\n  Results saved to: {csv_file}")

print("\n" + "="*65)
print("Depletion run complete.")
print("="*65)
