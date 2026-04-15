"""
ABR-1000 OpenMC Depletion Model — Equilibrium TRU-Recycled Metal Core
One full irradiation cycle as a function of burnup (328.5 EFPD).

References:
  [ANL]  ANL-AFCI-202, Argonne National Laboratory, September 2007.
  [NEA]  NEA-NSC-R(2015)9, Nuclear Energy Agency, February 2016.

"""

import math, csv
import numpy as np
import openmc
import openmc.deplete

# ======================================================================
# USER CONFIGURATION  
# ======================================================================
CHAIN_FILE     = 'chain_endfb80_sfr.xml'  # path to depletion chain file
POWER_W        = 1_000e6                  # reactor thermal power [W]
CYCLE_EFPD     = 328.5                    # cycle length [effective full-power days]
M_HM_KG        = 13_200.0                 # reference HM mass [kg] for burnup axis
THREADS        = 6                        # OpenMC thread count

# Transport settings per depletion step
DEPL_BATCHES   = 60
DEPL_INACTIVE  = 20
DEPL_PARTICLES = 20_000

# ======================================================================
# TEMPERATURES  
# ======================================================================
T_cool   = 432.5 + 273.15   # K  coolant + structural average
T_struct = T_cool
T_fuel   = 534.0 + 273.15   # K  average fuel temperature

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

H_lower_struct = 35.76    # cm
H_lower_refl   = 125.16   # cm
H_active       = 85.82    # cm 
H_bond_na      = 20.06    # cm
H_gas_plenum   = 101.01   # cm
H_upper_struct = 112.39   # cm

N_AX       = 5                            # axial fuel composition zones
dz         = H_active / N_AX              # 17.164 cm per zone
z_act_bot  = -H_active / 2.0              # -42.91 cm
z_act_top  = +H_active / 2.0              # +42.91 cm
z_ax       = [z_act_bot + i * dz for i in range(N_AX + 1)]
z_bond_top = z_act_top  + H_bond_na
z_plen_top = z_bond_top + H_gas_plenum
z_top      = z_plen_top + H_upper_struct
z_lrefl_bot= z_act_bot  - H_lower_refl
z_bot      = z_lrefl_bot - H_lower_struct

# ======================================================================
# VOLUME FRACTIONS  
# ======================================================================
VF_FUEL = 0.3900
VF_NA   = 0.3534
VF_HT9  = 0.2566

# ======================================================================
# PURE-MATERIAL NUMBER DENSITIES  
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
# FUEL-PIN NUMBER DENSITIES 
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
    'Zr':       [7.2802e-3] * N_AX,
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
    'Zr':       [7.2802e-3] * N_AX,
    'Mo_fp':    [8.1524e-4, 1.0174e-3, 1.0697e-3, 9.4870e-4, 6.6172e-4],
}
# Zone-averaged scalars for summary printout
_I = {k: (_mean(v) if isinstance(v, list) else v) for k, v in _I_zones.items()}
_O = {k: (_mean(v) if isinstance(v, list) else v) for k, v in _O_zones.items()}

# ======================================================================
# MATERIAL BUILDERS
# ======================================================================
_ACTINIDES = [
    'U234', 'U235', 'U236', 'U238',
    'Np237',
    'Pu238', 'Pu239', 'Pu240', 'Pu241', 'Pu242',
    'Am241', 'Am242_m1', 'Am243',
    'Cm242', 'Cm243', 'Cm244', 'Cm245', 'Cm246',
]

def make_active_zone(zone_data, iz, label, T, depletable=False):
    """Homogenised active-core material for axial zone iz (0..N_AX-1)."""
    m = openmc.Material(name=f'active_{label}_z{iz}')
    for nuc in _ACTINIDES:
        m.add_nuclide(nuc, zone_data[nuc][iz] * VF_FUEL)
    m.add_element('Zr', zone_data['Zr'][iz] * VF_FUEL)
    mo_total = zone_data['Mo_fp'][iz] * VF_FUEL + ND_HT9_Mo * VF_HT9
    m.add_element('Mo', mo_total)
    m.add_nuclide('Na23', ND_NA * VF_NA)
    m.add_element('Fe',   ND_HT9_Fe * VF_HT9)   # natural isotope expansion
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
    m.depletable   = depletable
    return m


def make_na_ht9(vf_na, vf_ht9, name, T):
    """Na + HT9 structural mixture (not depletable)."""
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
# NON-FUEL MATERIALS  (depletable=False)
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

_VF_sh_Na = 0.1710; _VF_sh_HT9 = 0.2968; _VF_sh_B4C = 1.0 - _VF_sh_Na - _VF_sh_HT9
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
    + (ND_HT9_Fe + ND_HT9_Ni + ND_HT9_Cr + ND_HT9_Mn + ND_HT9_Mo) * _VF_sh_HT9
    + (1.5018e-2 + 6.3609e-2 + 1.9657e-2) * _VF_sh_B4C)
radial_shield.temperature = T_struct

_VF_ca_Na = 0.2883; _VF_ca_HT9 = 0.2077; _VF_ca_B4C = 1.0 - _VF_ca_Na - _VF_ca_HT9
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
    + (ND_HT9_Fe + ND_HT9_Ni + ND_HT9_Cr + ND_HT9_Mn + ND_HT9_Mo) * _VF_ca_HT9
    + (5.3642e-2 + 2.8884e-2 + 2.0632e-2) * _VF_ca_B4C)
ctrl_absorber.temperature = T_struct

ctrl_empty = make_na_ht9(0.9074, 0.0926, 'ctrl_empty_duct', T_struct)

# ======================================================================
# CONTROL ASSEMBLY FRACTIONS
# ======================================================================
f_ci = 7.0  / 85.0     # inner zone: 7 control / 85 total assemblies
f_co = 12.0 / 114.0    # outer zone: 12 control / 114 total assemblies

def blend_two(m_base, f_base, m_ctrl, f_ctrl, name, T, depletable=False):
    """Volume-fraction blend of two materials."""
    assert abs(f_base + f_ctrl - 1.0) < 1e-9
    m = openmc.Material.mix_materials([m_base, m_ctrl], [f_base, f_ctrl], 'vo')
    m.name       = name
    m.temperature = T
    m.depletable  = depletable
    return m

# ======================================================================
# DEPLETABLE ACTIVE-CORE MATERIALS  (10 total: 5 axial x 2 radial)
# Rods fully withdrawn in active zone (BOEC convention for k_BOC).
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
    blend_two(inner_fuel_mats[iz], 1 - f_ci, ctrl_empty, f_ci,
              f'inner_ac_z{iz}', T_fuel, depletable=True)
    for iz in range(N_AX)
]
outer_ac_mats = [
    blend_two(outer_fuel_mats[iz], 1 - f_co, ctrl_empty, f_co,
              f'outer_ac_z{iz}', T_fuel, depletable=True)
    for iz in range(N_AX)
]

# ── Material volumes  [cm^3]  ─────────────────────────────────────────
#   V_inner_zone = pi * R_inner^2 * (H_active / N_AX)
#   V_outer_zone = pi * (R_outer^2 - R_inner^2) * (H_active / N_AX)
_V_inner_zone = math.pi * R_inner**2               * (H_active / N_AX)
_V_outer_zone = math.pi * (R_outer**2 - R_inner**2) * (H_active / N_AX)

for iz in range(N_AX):
    inner_ac_mats[iz].volume = _V_inner_zone   # cm^3
    outer_ac_mats[iz].volume = _V_outer_zone   # cm^3

# Structural zone blends (not depletable)
def zblend(base, f_ctrl, name):
    return blend_two(base, 1 - f_ctrl, ctrl_empty, f_ctrl, name, T_struct)

inner_ls = zblend(lower_struct_mat, f_ci, 'inner_ls')
outer_ls = zblend(lower_struct_mat, f_co, 'outer_ls')
inner_lr = zblend(lower_refl_mat,   f_ci, 'inner_lr')
outer_lr = zblend(lower_refl_mat,   f_co, 'outer_lr')
inner_bn = zblend(bond_na_mat,      f_ci, 'inner_bn')
outer_bn = zblend(bond_na_mat,      f_co, 'outer_bn')
inner_gp = zblend(plenum_mat,       f_ci, 'inner_gp')
outer_gp = zblend(plenum_mat,       f_co, 'outer_gp')
inner_us = zblend(upper_struct_mat, f_ci, 'inner_us')
outer_us = zblend(upper_struct_mat, f_co, 'outer_us')

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

z0     = openmc.ZPlane(z0=z_bot,       boundary_type='vacuum')
z1     = openmc.ZPlane(z0=z_lrefl_bot)
z2     = openmc.ZPlane(z0=z_act_bot)
z3     = openmc.ZPlane(z0=z_act_top)
z4     = openmc.ZPlane(z0=z_bond_top)
z5     = openmc.ZPlane(z0=z_plen_top)
z6     = openmc.ZPlane(z0=z_top,       boundary_type='vacuum')
zp_ax  = [openmc.ZPlane(z0=z_ax[i]) for i in range(1, N_AX)]  # 4 interior planes


def make_zone_cells(rad_reg, label, mat_ls, mat_lr, ac_mats,
                    mat_bn, mat_gp, mat_us):
    """Build 10 axial cells for one radial annulus."""
    cells = []
    cells.append(openmc.Cell(name=f'{label}_ls',  fill=mat_ls,
                             region=rad_reg & +z0 & -z1))
    cells.append(openmc.Cell(name=f'{label}_lr',  fill=mat_lr,
                             region=rad_reg & +z1 & -z2))
    # 5 active sub-zones
    cells.append(openmc.Cell(name=f'{label}_ac0', fill=ac_mats[0],
                             region=rad_reg & +z2 & -zp_ax[0]))
    for i in range(1, N_AX - 1):
        cells.append(openmc.Cell(name=f'{label}_ac{i}', fill=ac_mats[i],
                                 region=rad_reg & +zp_ax[i-1] & -zp_ax[i]))
    cells.append(openmc.Cell(name=f'{label}_ac4', fill=ac_mats[4],
                             region=rad_reg & +zp_ax[N_AX-2] & -z3))
    cells.append(openmc.Cell(name=f'{label}_bn',  fill=mat_bn,
                             region=rad_reg & +z3 & -z4))
    cells.append(openmc.Cell(name=f'{label}_gp',  fill=mat_gp,
                             region=rad_reg & +z4 & -z5))
    cells.append(openmc.Cell(name=f'{label}_us',  fill=mat_us,
                             region=rad_reg & +z5 & -z6))
    return cells   # [ls, lr, ac0, ac1, ac2, ac3, ac4, bn, gp, us]


inner_cells = make_zone_cells(
    -cyl_inner, 'inner',
    inner_ls, inner_lr, inner_ac_mats, inner_bn, inner_gp, inner_us)
outer_cells = make_zone_cells(
    +cyl_inner & -cyl_outer, 'outer',
    outer_ls, outer_lr, outer_ac_mats, outer_bn, outer_gp, outer_us)

cell_refl   = openmc.Cell(name='radial_refl',   fill=radial_refl_mat,
                          region=+cyl_outer & -cyl_refl  & +z0 & -z6)
cell_shield = openmc.Cell(name='radial_shield',  fill=radial_shield,
                          region=+cyl_refl  & -cyl_shield & +z0 & -z6)

universe = openmc.Universe(
    cells=inner_cells + outer_cells + [cell_refl, cell_shield])
geometry = openmc.Geometry(universe)

# ======================================================================
# SETTINGS  (eigenvalue, used at every depletion transport step)
# ======================================================================
settings = openmc.Settings()
settings.run_mode  = 'eigenvalue'
settings.batches   = DEPL_BATCHES
settings.inactive  = DEPL_INACTIVE
settings.particles = DEPL_PARTICLES

settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-R_outer, -R_outer, z_act_bot],
        [ R_outer,  R_outer, z_act_top],
    ),
    constraints={'fissionable': True},
)
entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left  = [-R_outer, -R_outer, z_act_bot]
entropy_mesh.upper_right = [ R_outer,  R_outer, z_act_top]
entropy_mesh.dimension   = [20, 20, 20]
settings.entropy_mesh = entropy_mesh
settings.temperature  = {'default': T_struct, 'method': 'interpolation'}

# ======================================================================
# TALLIES  (physics output recorded at every depletion step)
# ======================================================================
grpstr33 = [
    1.00000E-05, 4.17458E-01, 5.31578E-01, 3.92786E+00, 8.31528E+00,
    1.37096E+01, 2.26033E+01, 3.72665E+01, 6.14421E+01, 1.01301E+02,
    1.67017E+02, 2.75364E+02, 4.53999E+02, 7.48518E+02, 1.23410E+03,
    2.03468E+03, 3.35462E+03, 5.53084E+03, 9.11881E+03, 1.50344E+04,
    2.47875E+04, 4.08677E+04, 6.73794E+04, 1.11090E+05, 1.83156E+05,
    3.01974E+05, 4.97871E+05, 8.20850E+05, 1.35335E+06, 2.23130E+06,
    3.67879E+06, 6.06531E+06, 1.00000E+07, 2.00000E+07
]

inner_ac_cells = inner_cells[2:7]   
outer_ac_cells = outer_cells[2:7]

fuel_filt  = openmc.CellFilter(inner_ac_cells + outer_ac_cells)
inner_filt = openmc.CellFilter(inner_ac_cells)
outer_filt = openmc.CellFilter(outer_ac_cells)
all_filt   = openmc.CellFilter(
    inner_cells + outer_cells + [cell_refl, cell_shield])
neut_filt  = openmc.ParticleFilter(['neutron'])
E_filt     = openmc.EnergyFilter(grpstr33)
leak_filt  = openmc.SurfaceFilter([cyl_shield, z0, z6])

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

t = openmc.Tally(name='neutron_balance')
t.filters = [neut_filt, all_filt]
t.scores  = ['absorption', 'nu-fission', '(n,2n)', '(n,3n)']
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
model = openmc.Model(
    geometry=geometry,
    materials=materials,
    settings=settings,
    tallies=tallies,
)

# ======================================================================
# DEPLETION TIMESTEP SCHEDULE
# ======================================================================
_t_pts    = np.geomspace(0.5, CYCLE_EFPD, 13)   
timesteps = list(np.diff(_t_pts))
timesteps[-1] += CYCLE_EFPD - sum(timesteps)     # absorb floating-point error
timesteps = [float(dt) for dt in timesteps]

_cum_d = np.cumsum([0.0] + timesteps)
_cum_bu = POWER_W / 1e6 * _cum_d / (M_HM_KG / 1e3)

print("\n" + "=" * 68)
print("ABR-1000  Full Cycle Depletion Run")
print("=" * 68)
print(f"  Chain file      : {CHAIN_FILE}")
print(f"  Power           : {POWER_W/1e6:.0f} MWt")
print(f"  Cycle length    : {CYCLE_EFPD:.1f} EFPD  (12 months x 90% CF)")
print(f"  Reference HM    : {M_HM_KG:.0f} kg  ([ANL] Table II.1-3)")
print(f"  Transport/step  : {DEPL_BATCHES} batches x {DEPL_PARTICLES:,} particles "
      f"({DEPL_INACTIVE} inactive)")
print(f"  Depletable mats : {N_AX * 2}  (5 axial x 2 radial zones)")
print()
print(f"  Timestep schedule ({len(timesteps)} steps):")
print(f"  {'#':>3}  {'dt [d]':>8}  {'cum [d]':>8}  {'BU [MWd/kg]':>12}")
for i, dt in enumerate(timesteps):
    print(f"  {i+1:3d}  {dt:8.2f}  {_cum_d[i+1]:8.1f}  {_cum_bu[i+1]:12.1f}")
print()
print(f"  Expected k swing : ~2210 pcm ([NEA] Table 4.3 mean +/- 422 pcm)")
print(f"  Expected TRU loss: ~82 kg/cycle ([ANL] Table II.1-4)")
print("=" * 68 + "\n")

# ======================================================================
# RUN DEPLETION
# ======================================================================
operator = openmc.deplete.CoupledOperator(model, CHAIN_FILE)

integrator = openmc.deplete.CECMIntegrator(
    operator,
    timesteps,
    POWER_W,
    timestep_units='d',
)
integrator.integrate()

# ======================================================================
# POST-PROCESSING
# ======================================================================
print("\n" + "=" * 68)
print("POST-PROCESSING DEPLETION RESULTS")
print("=" * 68)

results   = openmc.deplete.Results('depletion_results.h5')
times_d   = results.get_times('d')
keffs, kuncs = zip(*results.get_keff())

# Material IDs of all depletable (active-core) materials
depl_ids = set(m.id for m in inner_ac_mats + outer_ac_mats)

def total_mass_kg(si, nuclide):
    """Sum nuclide mass [kg] across all depletable materials at step si."""
    _, mats = results[si]
    g = 0.0
    for mat in mats:
        if mat.id in depl_ids:
            try:
                g += mat.get_mass(nuclide)
            except Exception:
                pass
    return g / 1e3

# Nuclides tracked in the summary
TRACK    = ['U235','U238','Np237','Pu238','Pu239','Pu240','Pu241','Pu242',
            'Am241','Am243','Cm244','Cm245']
TRU_NUCS = ['Np237','Pu238','Pu239','Pu240','Pu241','Pu242',
            'Am241','Am242_m1','Am243','Cm242','Cm243','Cm244','Cm245','Cm246']

rows = []
for si in range(len(results)):
    t_d  = times_d[si]
    bu   = POWER_W / 1e6 * t_d / (M_HM_KG / 1e3)
    k    = keffs[si]
    ku   = kuncs[si]
    mass = {n: total_mass_kg(si, n) for n in TRACK}
    mass['Am242_m1'] = total_mass_kg(si, 'Am242_m1')
    mass['Cm242']    = total_mass_kg(si, 'Cm242')
    mass['Cm243']    = total_mass_kg(si, 'Cm243')
    mass['Cm246']    = total_mass_kg(si, 'Cm246')
    m_TRU = sum(total_mass_kg(si, n) for n in TRU_NUCS)
    m_HM  = mass['U235'] + mass['U238'] + m_TRU
    rows.append({
        'step': si, 'time_d': t_d, 'burnup_MWd_kg': bu,
        'keff': k, 'keff_1sig': ku,
        **{f'm_{n}_kg': mass.get(n, 0.0) for n in TRACK},
        'm_TRU_kg': m_TRU, 'm_HM_kg': m_HM,
    })

# k-eff vs burnup table
print(f"\n  {'Step':>4}  {'Time (d)':>9}  {'BU (MWd/kg)':>12}  "
      f"{'k-eff':>8}  {'+/-1sig':>8}")
print("  " + "-" * 55)
for r in rows:
    print(f"  {r['step']:4d}  {r['time_d']:9.1f}  {r['burnup_MWd_kg']:12.1f}"
          f"  {r['keff']:.5f}  {r['keff_1sig']:.5f}")

# Reactivity swing
k0   = rows[0]['keff'];  kf = rows[-1]['keff']
dk   = (k0 - kf) / (k0 * kf) * 1e5
print(f"\n  Burnup reactivity swing (Drho_cycle): {dk:.0f} pcm")
print(f"    k_BOC = {k0:.5f},  k_EOC = {kf:.5f}")
print(f"    Benchmark mean: 2210 +/- 422 pcm  ([NEA] Table 4.3)")

# Actinide inventory change
r0 = rows[0];  rn = rows[-1]
print(f"\n  Actinide inventory change  (BOC -> EOC):")
print(f"  {'Nuclide':<10}  {'BOC [kg]':>10}  {'EOC [kg]':>10}  {'Delta [kg]':>11}")
for nuc in TRACK:
    k  = f'm_{nuc}_kg'
    m0 = r0[k];  mf = rn[k]
    print(f"  {nuc:<10}  {m0:10.2f}  {mf:10.2f}  {mf-m0:+11.2f}")
print(f"  {'TRU total':<10}  {r0['m_TRU_kg']:10.2f}  {rn['m_TRU_kg']:10.2f}"
      f"  {rn['m_TRU_kg']-r0['m_TRU_kg']:+11.2f}")
print(f"  {'HM total':<10}  {r0['m_HM_kg']:10.2f}  {rn['m_HM_kg']:10.2f}"
      f"  {rn['m_HM_kg']-r0['m_HM_kg']:+11.2f}")

tru_consumed = r0['m_TRU_kg'] - rn['m_TRU_kg']
tru_rate     = tru_consumed * 365.25 / CYCLE_EFPD
print(f"\n  TRU destroyed this cycle : {tru_consumed:.1f} kg")
print(f"  TRU consumption rate     : {tru_rate:.1f} kg/year")
print(f"    Benchmark reference    : 81.6 kg/year  ([ANL] Table II.1-4)")

# Save CSV
csv_path = 'ABR1000_depletion_summary.csv'
with open(csv_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(f"\n  Results saved -> {csv_path}")
print("=" * 68 + "\n")
