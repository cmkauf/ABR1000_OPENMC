"""
ABR-1000 OpenMC Model; Equilibrium TRU-Recycled Metal Core
HETEROGENEOUS AXIAL LAYOUT  v4  —  direct benchmark number densities

References:
  [ANL]  ANL-AFCI-202, Argonne National Laboratory, September 2007.
  [NEA]  NEA-NSC-R(2015)9, Nuclear Energy Agency, February 2016.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CHANGE LOG
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  v1  k≈1.046  Entire ±150 cm axial extent filled with active fuel.
  v2  k≈0.950  Axial zones added but fabrication vol fracs used (fuel=29.2%).
  v3  k≈0.970  Operating vol fracs fixed (fuel=39%) but fuel isotopics still
               built from the fresh LWR-SF TRU feed vector.  The equilibrium
               cycle fuel is heavily depleted and has a substantially different
               isotopic mix (lower Pu-239 fraction, more Pu-240/Am/Cm) that
               is less reactive per unit mass than fresh TRU.  This suppressed
               k by ~650 pcm relative to the benchmark composition.
  v4  (this)   Active-core materials specified directly from [NEA] Tables 2.22
               and 2.23 zone-averaged fuel pin number densities, homogenised
               with Table 2.20 volume fractions and Table 2.21 structural
               densities.  All other zones (lower structure, lower reflector,
               displaced bond Na, gas plenum, upper structure) use the same
               architecture as v3.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ACTIVE-CORE MATERIAL CONSTRUCTION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  The homogenised number density N_hom for any nuclide is:

      N_hom = vf_fuel × N_pin  +  vf_Na × N_Na  +  vf_HT9 × N_HT9

  where:
    vf_fuel = 0.39,  vf_Na = 0.3534,  vf_HT9 = 0.2566  (Table 2.20)
    N_pin  = zone-averaged PIN number density from Tables 2.22/2.23
    N_Na   = 2.2272e-2 a/b-cm  (Table 2.21)
    N_HT9  = per-element densities from Table 2.21

  The axial variation across the 5 zones is collapsed to a single zone-
  average for this radially-homogenised cylindrical model.  The benchmark
  participants who used homogenised models (ANL-3, IKE-2) obtained
  k_BOC ≈ 1.037–1.043, consistent with the mean of 1.0355 ± 78 pcm.

  Mo in the fuel is the pseudo-fission-product.  Mo in HT9 is structural.
  Both contribute to Mo92 in OpenMC; they are summed before add_nuclide.

  Zr in the benchmark uses natural Zr (Zr90 is listed as a proxy for all
  Zr isotopes in the table; we use add_element('Zr') scaled to the same
  total atom density as listed for Zr90 in the table, which is equivalent
  since the table treats it as a single pseudo-nuclide).
"""

import math
import openmc

# ======================================================================
# TEMPERATURES  [NEA] Table 2.15
# ======================================================================
T_cool   = 432.5 + 273.15   # K
T_struct = T_cool
T_fuel   = 534.0 + 273.15   # K

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

H_active  = 85.82   # cm active fuel (operating, Table 2.16)
H_lower_r = 125.16  # cm lower reflector
H_lower_s = 35.76   # cm lower structure
H_bond_na = 20.06   # cm displaced bond sodium
H_gas_pl  = 101.01  # cm gas plenum
H_upper_s = 112.39  # cm upper structure

z_fuel_bot = -H_active / 2.0
z_fuel_top = +H_active / 2.0
z_lr_bot   = z_fuel_bot - H_lower_r
z_ls_bot   = z_lr_bot   - H_lower_s
z_bond_top = z_fuel_top + H_bond_na
z_pl_top   = z_bond_top + H_gas_pl
z_us_top   = z_pl_top   + H_upper_s

# ======================================================================
# VOLUME FRACTIONS  [NEA] Table 2.20
# ======================================================================
VF_FUEL = 0.3900
VF_NA   = 0.3534
VF_HT9  = 0.2566
# Active core: 0.3900 + 0.3534 + 0.2566 = 1.0000 ✓

# ======================================================================
# PURE-MATERIAL NUMBER DENSITIES  [NEA] Table 2.21  (a/b-cm)
# ======================================================================
ND_NA   = 2.2272e-02   # Na coolant

# HT9 element number densities
ND_HT9_Fe  = 6.9715e-02
ND_HT9_Ni  = 4.2984e-04
ND_HT9_Cr  = 1.0366e-02
ND_HT9_Mn  = 4.5921e-04
ND_HT9_Mo  = 4.9007e-04  # structural Mo in HT9

# ======================================================================
# ZONE-AVERAGED FUEL-PIN NUMBER DENSITIES  [NEA] Tables 2.22 / 2.23
# Five axial zones; simple arithmetic mean used for the homogenised model.
# ======================================================================
import statistics as _st
_mean = lambda lst: sum(lst) / len(lst)

# Inner core
_I = {
    'U234':     _mean([1.1369e-6, 1.0856e-6, 1.0727e-6, 1.1028e-6, 1.1759e-6]),
    'U235':     _mean([3.0421e-5, 2.9338e-5, 2.8961e-5, 3.0070e-5, 3.2571e-5]),
    'U236':     _mean([2.4896e-6, 2.5117e-6, 2.5536e-6, 2.3779e-6, 2.0226e-6]),
    'U238':     _mean([1.9613e-2, 1.9474e-2, 1.9433e-2, 1.9550e-2, 1.9801e-2]),
    'Np237':    _mean([4.6686e-5, 4.6962e-5, 4.6782e-5, 4.7603e-5, 4.8895e-5]),
    'Pu238':    _mean([1.1695e-4, 1.1284e-4, 1.1196e-4, 1.1370e-4, 1.1829e-4]),
    'Pu239':    _mean([2.2076e-3, 2.1814e-3, 2.1754e-3, 2.1813e-3, 2.2011e-3]),
    'Pu240':    _mean([1.3244e-3, 1.2955e-3, 1.2902e-3, 1.2986e-3, 1.3248e-3]),
    'Pu241':    _mean([1.9375e-4, 1.8610e-4, 1.8518e-4, 1.8537e-4, 1.8845e-4]),
    'Pu242':    _mean([2.9277e-4, 2.8911e-4, 2.8818e-4, 2.9038e-4, 2.9569e-4]),
    'Am241':    _mean([1.0791e-4, 1.0465e-4, 1.0353e-4, 1.0686e-4, 1.1421e-4]),
    'Am242_m1': _mean([9.2989e-6, 9.0848e-6, 9.0224e-6, 9.1756e-6, 9.4890e-6]),
    'Am243':    _mean([1.0017e-4, 9.8324e-5, 9.7993e-5, 9.8630e-5, 1.0032e-4]),
    'Cm242':    _mean([5.6250e-6, 5.8208e-6, 5.9476e-6, 5.4901e-6, 4.5416e-6]),
    'Cm243':    _mean([5.4321e-7, 5.0246e-7, 5.0136e-7, 4.8876e-7, 4.8480e-7]),
    'Cm244':    _mean([6.7240e-5, 6.5722e-5, 6.5622e-5, 6.5349e-5, 6.5394e-5]),
    'Cm245':    _mean([1.7397e-5, 1.6743e-5, 1.6663e-5, 1.6696e-5, 1.7026e-5]),
    'Cm246':    _mean([9.2285e-6, 9.1426e-6, 9.1307e-6, 9.1364e-6, 9.1805e-6]),
    'Zr':       7.2802e-3,   # natural Zr (as-listed; use add_element below)
    'Mo_fp':    _mean([9.2873e-4, 1.1464e-3, 1.2031e-3, 1.0625e-3, 7.4065e-4]),
}

# Outer core
_O = {
    'U234':     _mean([1.6317e-6, 1.5766e-6, 1.5638e-6, 1.5894e-6, 1.6552e-6]),
    'U235':     _mean([3.0822e-5, 2.9870e-5, 2.9561e-5, 3.0391e-5, 3.2250e-5]),
    'U236':     _mean([1.7881e-6, 1.8534e-6, 1.8941e-6, 1.7528e-6, 1.4710e-6]),
    'U238':     _mean([1.8244e-2, 1.8144e-2, 1.8115e-2, 1.8191e-2, 1.8359e-2]),
    'Np237':    _mean([9.8244e-5, 9.7300e-5, 9.6775e-5, 9.8481e-5, 1.0175e-4]),
    'Pu238':    _mean([1.6436e-4, 1.6026e-4, 1.5949e-4, 1.6063e-4, 1.6416e-4]),
    'Pu239':    _mean([2.8147e-3, 2.7664e-3, 2.7538e-3, 2.7786e-3, 2.8416e-3]),
    'Pu240':    _mean([1.7467e-3, 1.7191e-3, 1.7135e-3, 1.7231e-3, 1.7508e-3]),
    'Pu241':    _mean([2.8976e-4, 2.8138e-4, 2.8012e-4, 2.8135e-4, 2.8697e-4]),
    'Pu242':    _mean([4.0754e-4, 4.0412e-4, 4.0321e-4, 4.0530e-4, 4.1028e-4]),
    'Am241':    _mean([1.8607e-4, 1.8127e-4, 1.7970e-4, 1.8397e-4, 1.9339e-4]),
    'Am242_m1': _mean([1.2185e-5, 1.2045e-5, 1.2021e-5, 1.2039e-5, 1.2064e-5]),
    'Am243':    _mean([1.3234e-4, 1.3019e-4, 1.2985e-4, 1.3036e-4, 1.3206e-4]),
    'Cm242':    _mean([6.4688e-6, 6.8630e-6, 7.0553e-6, 6.4446e-6, 5.1976e-6]),
    'Cm243':    _mean([6.3471e-7, 6.0893e-7, 6.0901e-7, 5.9753e-7, 5.9372e-7]),
    'Cm244':    _mean([8.0107e-5, 7.8889e-5, 7.8847e-5, 7.8479e-5, 7.8359e-5]),
    'Cm245':    _mean([2.0200e-5, 1.9678e-5, 1.9613e-5, 1.9635e-5, 1.9913e-5]),
    'Cm246':    _mean([1.0443e-5, 1.0371e-5, 1.0361e-5, 1.0367e-5, 1.0410e-5]),
    'Zr':       7.2802e-3,
    'Mo_fp':    _mean([8.1524e-4, 1.0174e-3, 1.0697e-3, 9.4870e-4, 6.6172e-4]),
}

# ======================================================================
# BUILD HOMOGENISED ACTIVE-CORE OPENMC MATERIALS
# ======================================================================

def make_active_core(pin_nd, zone_label, T):
    """
    Create a homogenised active-core OpenMC Material from fuel-pin number
    densities (pin_nd dict) scaled by VF_FUEL, plus Na and HT9 contributions.
    """
    m = openmc.Material(name=f'active_{zone_label}')

    # ---- Fuel nuclides (scaled by volume fraction) ----
    actinides = [
        'U234','U235','U236','U238',
        'Np237',
        'Pu238','Pu239','Pu240','Pu241','Pu242',
        'Am241','Am242_m1','Am243',
        'Cm242','Cm243','Cm244','Cm245','Cm246',
    ]
    for nuc in actinides:
        m.add_nuclide(nuc, pin_nd[nuc] * VF_FUEL)

    # Zr: use add_element so all natural isotopes are included; the benchmark
    # lists a single "Zr" entry — scale to the same total atom density.
    # OpenMC add_element with 'ao' accepts the total element atom density.
    m.add_element('Zr', pin_nd['Zr'] * VF_FUEL)

    # Mo = pseudo-FP Mo (from fuel) + structural Mo (from HT9)
    mo_total = pin_nd['Mo_fp'] * VF_FUEL + ND_HT9_Mo * VF_HT9
    m.add_element('Mo', mo_total)

    # ---- Coolant sodium ----
    m.add_nuclide('Na23', ND_NA * VF_NA)

    # ---- HT9 structural elements (excluding Mo, handled above) ----
    m.add_nuclide('Fe56', ND_HT9_Fe * VF_HT9)   # Fe56 dominant isotope
    m.add_nuclide('Ni58', ND_HT9_Ni * VF_HT9)   # Ni58 dominant isotope
    m.add_nuclide('Cr52', ND_HT9_Cr * VF_HT9)   # Cr52 dominant isotope
    m.add_nuclide('Mn55', ND_HT9_Mn * VF_HT9)   # Mn55 only stable isotope

    # Total atom density
    total = (
        sum(pin_nd[nuc] * VF_FUEL for nuc in actinides)
        + pin_nd['Zr']    * VF_FUEL
        + mo_total
        + ND_NA           * VF_NA
        + (ND_HT9_Fe + ND_HT9_Ni + ND_HT9_Cr + ND_HT9_Mn) * VF_HT9
        # Mo already counted above
    )
    m.set_density('atom/b-cm', total)
    m.temperature = T
    return m

inner_core_active = make_active_core(_I, 'inner', T_fuel)
outer_core_active = make_active_core(_O, 'outer', T_fuel)

# ======================================================================
# NON-FUEL ZONE MATERIALS
# Constructed with explicit add_nuclide/add_element calls from Table 2.21
# volume fractions to avoid any mix_materials rounding issues.
# ======================================================================

def make_na_ht9(vf_na, vf_ht9, name, T):
    """Pure Na + HT9 mixture at given volume fractions (must sum to 1)."""
    assert abs(vf_na + vf_ht9 - 1.0) < 1e-6, f"{name}: vf sum = {vf_na+vf_ht9}"
    m = openmc.Material(name=name)
    m.add_nuclide('Na23', ND_NA   * vf_na)
    m.add_nuclide('Fe56', ND_HT9_Fe * vf_ht9)
    m.add_nuclide('Ni58', ND_HT9_Ni * vf_ht9)
    m.add_nuclide('Cr52', ND_HT9_Cr * vf_ht9)
    m.add_nuclide('Mn55', ND_HT9_Mn * vf_ht9)
    m.add_element('Mo',   ND_HT9_Mo * vf_ht9)
    total = (ND_NA * vf_na
             + (ND_HT9_Fe + ND_HT9_Ni + ND_HT9_Cr + ND_HT9_Mn + ND_HT9_Mo) * vf_ht9)
    m.set_density('atom/b-cm', total)
    m.temperature = T
    return m

# Lower reflector / upper structure: Na 35.34% + HT9 64.66%
lower_refl_mat   = make_na_ht9(0.3534, 0.6466, 'lower_reflector', T_struct)
upper_struct_mat = make_na_ht9(0.3534, 0.6466, 'upper_structure', T_struct)

# Displaced bond sodium zone: Na 74.34% + HT9 25.66%
bond_na_mat = make_na_ht9(0.7434, 0.2566, 'displaced_bond_na', T_struct)

# Gas plenum: Na 35.34% + HT9 25.66% + He 39.00%
# He is negligible in fast spectrum; Na+HT9 renormalized to 1.0:
#   Na  = 35.34/61.00 = 57.93%,  HT9 = 25.66/61.00 = 42.07%
plenum_mat = make_na_ht9(0.5793, 0.4207, 'gas_plenum', T_struct)

# Lower structure: Na 70% + SS-316 30%
# SS-316 approximated by HT9 composition (conservative; both are Fe-Cr-Ni)
# Table 2.21 lower structure number densities used directly:
lower_struct_mat = openmc.Material(name='lower_structure')
lower_struct_mat.add_nuclide('Na23',  1.5591e-02)
lower_struct_mat.add_nuclide('Fe56',  1.5878e-02)
lower_struct_mat.add_nuclide('Ni58',  3.2604e-03)
lower_struct_mat.add_nuclide('Cr52',  3.2355e-03)
lower_struct_mat.add_nuclide('Mn55',  5.0846e-04)
lower_struct_mat.add_element('Mo',    4.3524e-04)
lower_struct_mat.set_density('atom/b-cm',
    1.5591e-2 + 1.5878e-2 + 3.2604e-3 + 3.2355e-3 + 5.0846e-4 + 4.3524e-4)
lower_struct_mat.temperature = T_struct

# Radial reflector: Na 15.50% + HT9 84.50%
radial_refl_mat = make_na_ht9(0.1550, 0.8450, 'radial_reflector', T_struct)

# Radial shield: Na 17.10% + HT9 29.68% + natB4C 53.22% (renorm 0.5323→0.5322)
# B4C from Table 2.21: C 1.9657e-2, B10 1.5018e-2, B11 6.3609e-2 a/b-cm
_VF_shield_Na   = 0.1710
_VF_shield_HT9  = 0.2968
_VF_shield_B4C  = 1.0 - _VF_shield_Na - _VF_shield_HT9  # = 0.5322
radial_shield = openmc.Material(name='radial_shield')
radial_shield.add_nuclide('Na23',  ND_NA        * _VF_shield_Na)
radial_shield.add_nuclide('Fe56',  ND_HT9_Fe    * _VF_shield_HT9)
radial_shield.add_nuclide('Ni58',  ND_HT9_Ni    * _VF_shield_HT9)
radial_shield.add_nuclide('Cr52',  ND_HT9_Cr    * _VF_shield_HT9)
radial_shield.add_nuclide('Mn55',  ND_HT9_Mn    * _VF_shield_HT9)
radial_shield.add_element('Mo',    ND_HT9_Mo    * _VF_shield_HT9)
radial_shield.add_nuclide('B10',   1.5018e-02   * _VF_shield_B4C)
radial_shield.add_nuclide('B11',   6.3609e-02   * _VF_shield_B4C)
radial_shield.add_nuclide('C12',   1.9657e-02   * _VF_shield_B4C)
radial_shield.set_density('atom/b-cm',
    ND_NA * _VF_shield_Na
    + (ND_HT9_Fe + ND_HT9_Ni + ND_HT9_Cr + ND_HT9_Mn + ND_HT9_Mo) * _VF_shield_HT9
    + (1.9657e-2 + 1.5018e-2 + 6.3609e-2) * _VF_shield_B4C)
radial_shield.temperature = T_struct

# Control absorber: Na 28.83% + HT9 20.77% + enrB4C 50.40%
# Enriched B4C from Table 2.21: C 2.0632e-2, B10 5.3642e-2, B11 2.8884e-2
_VF_ctrl_Na  = 0.2883
_VF_ctrl_HT9 = 0.2077
_VF_ctrl_B4C = 1.0 - _VF_ctrl_Na - _VF_ctrl_HT9  # = 0.5040
ctrl_absorber = openmc.Material(name='ctrl_absorber')
ctrl_absorber.add_nuclide('Na23',  ND_NA      * _VF_ctrl_Na)
ctrl_absorber.add_nuclide('Fe56',  ND_HT9_Fe  * _VF_ctrl_HT9)
ctrl_absorber.add_nuclide('Ni58',  ND_HT9_Ni  * _VF_ctrl_HT9)
ctrl_absorber.add_nuclide('Cr52',  ND_HT9_Cr  * _VF_ctrl_HT9)
ctrl_absorber.add_nuclide('Mn55',  ND_HT9_Mn  * _VF_ctrl_HT9)
ctrl_absorber.add_element('Mo',    ND_HT9_Mo  * _VF_ctrl_HT9)
ctrl_absorber.add_nuclide('B10',   5.3642e-02 * _VF_ctrl_B4C)
ctrl_absorber.add_nuclide('B11',   2.8884e-02 * _VF_ctrl_B4C)
ctrl_absorber.add_nuclide('C12',   2.0632e-02 * _VF_ctrl_B4C)
ctrl_absorber.set_density('atom/b-cm',
    ND_NA * _VF_ctrl_Na
    + (ND_HT9_Fe + ND_HT9_Ni + ND_HT9_Cr + ND_HT9_Mn + ND_HT9_Mo) * _VF_ctrl_HT9
    + (2.0632e-2 + 5.3642e-2 + 2.8884e-2) * _VF_ctrl_B4C)
ctrl_absorber.temperature = T_struct

# Empty duct in control: Na 90.74% + HT9 9.26%
ctrl_empty = make_na_ht9(0.9074, 0.0926, 'ctrl_empty_duct', T_struct)

# ======================================================================
# RADIAL BLENDING OF CONTROL ASSEMBLIES INTO FUEL ZONES
# (assembly-count weighted; control absorber zone co-located with active core)
# Inner zone: 85 total (78 fuel + 3 primary-row4 + 4 secondary) => 7/85 ctrl
# Outer zone: 114 total (102 fuel + 12 primary-row7)             => 12/114 ctrl
# ======================================================================
f_ci = 7.0  / 85.0
f_co = 12.0 / 114.0

def blend_two(m_base, f_base, m_ctrl, f_ctrl, name, T):
    """Linear blend of two pre-built OpenMC materials by atom density."""
    # Use mix_materials with 'vo' — both materials were built with
    # set_density('atom/b-cm'), so 'vo' here acts as a fractional weight
    # on the atom densities.  Fractions must sum to 1.
    assert abs(f_base + f_ctrl - 1.0) < 1e-9
    m = openmc.Material.mix_materials([m_base, m_ctrl], [f_base, f_ctrl], 'vo')
    m.name = name
    m.temperature = T
    return m

# Active-core blended — BOEC: rods WITHDRAWN, use ctrl_empty (no B4C) in active zone.
inner_ac = blend_two(inner_core_active, 1-f_ci, ctrl_empty, f_ci,
                     'inner_active_blended', T_fuel)
outer_ac = blend_two(outer_core_active, 1-f_co, ctrl_empty, f_co,
                     'outer_active_blended', T_fuel)

# Non-active axial zones: blend ctrl_empty into structural material
# (control empty duct occupies same fraction of the cross section)
def zblend(base, f_ctrl, ctrl_mat, name):
    return blend_two(base, 1-f_ctrl, ctrl_mat, f_ctrl, name, T_struct)

inner_ls  = zblend(lower_struct_mat, f_ci, lower_struct_mat, 'inner_ls')   # same mat
outer_ls  = zblend(lower_struct_mat, f_co, lower_struct_mat, 'outer_ls')
inner_lr  = zblend(lower_refl_mat,   f_ci, ctrl_empty, 'inner_lr')
outer_lr  = zblend(lower_refl_mat,   f_co, ctrl_empty, 'outer_lr')
inner_bn  = zblend(bond_na_mat,      f_ci, ctrl_empty, 'inner_bn')
outer_bn  = zblend(bond_na_mat,      f_co, ctrl_empty, 'outer_bn')
inner_gp  = zblend(plenum_mat,       f_ci, ctrl_empty, 'inner_gp')
outer_gp  = zblend(plenum_mat,       f_co, ctrl_empty, 'outer_gp')
inner_us  = zblend(upper_struct_mat, f_ci, ctrl_empty, 'inner_us')
outer_us  = zblend(upper_struct_mat, f_co, ctrl_empty, 'outer_us')

# =============================================================================
# MATERIALS FILE
# =============================================================================
materials = openmc.Materials([
    inner_ac, outer_ac,
    inner_ls, outer_ls,
    inner_lr, outer_lr,
    inner_bn, outer_bn,
    inner_gp, outer_gp,
    inner_us, outer_us,
    radial_refl_mat,
    radial_shield,
])
materials.export_to_xml()

# ======================================================================
# GEOMETRY
# ======================================================================
cyl_inner  = openmc.ZCylinder(r=R_inner)
cyl_outer  = openmc.ZCylinder(r=R_outer)
cyl_refl   = openmc.ZCylinder(r=R_refl)
cyl_shield = openmc.ZCylinder(r=R_shield, boundary_type='vacuum')

z0 = openmc.ZPlane(z0=z_ls_bot,   boundary_type='vacuum')
z1 = openmc.ZPlane(z0=z_lr_bot)
z2 = openmc.ZPlane(z0=z_fuel_bot)
z3 = openmc.ZPlane(z0=z_fuel_top)
z4 = openmc.ZPlane(z0=z_bond_top)
z5 = openmc.ZPlane(z0=z_pl_top)
z6 = openmc.ZPlane(z0=z_us_top,   boundary_type='vacuum')

def axial_cells(rad_reg, mats, prefix):
    """Build 6 axial cells for one radial annulus."""
    layers = [(z0,z1,mats[0],'ls'), (z1,z2,mats[1],'lr'),
              (z2,z3,mats[2],'ac'), (z3,z4,mats[3],'bn'),
              (z4,z5,mats[4],'gp'), (z5,z6,mats[5],'us')]
    return [openmc.Cell(name=f'{prefix}_{s}', fill=m,
                        region=rad_reg & +zlo & -zhi)
            for zlo, zhi, m, s in layers]

inner_cells = axial_cells(-cyl_inner,
    [inner_ls, inner_lr, inner_ac, inner_bn, inner_gp, inner_us], 'inner')
outer_cells = axial_cells(+cyl_inner & -cyl_outer,
    [outer_ls, outer_lr, outer_ac, outer_bn, outer_gp, outer_us], 'outer')

cell_refl   = openmc.Cell(name='radial_refl',   fill=radial_refl_mat,
                          region=+cyl_outer & -cyl_refl  & +z0 & -z6)
cell_shield = openmc.Cell(name='radial_shield',  fill=radial_shield,
                          region=+cyl_refl  & -cyl_shield & +z0 & -z6)

universe = openmc.Universe(
    cells=inner_cells + outer_cells + [cell_refl, cell_shield])
geometry = openmc.Geometry(universe)
geometry.export_to_xml()

# =============================================================================
# PLOTS
# =============================================================================
plot_xy = openmc.Plot.from_geometry(geometry)
plot_xy.filename = 'ABR1000_v4_xy'
plot_xy.basis = 'xy'; plot_xy.origin = (0., 0., 0.)
plot_xy.width = (2.2*R_shield, 2.2*R_shield); plot_xy.pixels = (2000, 2000)
plot_xy.color_by = 'material'

plot_yz = openmc.Plot.from_geometry(geometry)
plot_yz.filename = 'ABR1000_v4_yz'
plot_yz.basis = 'yz'; plot_yz.origin = (0., 0., (z_ls_bot+z_us_top)/2.)
plot_yz.width = (2.2*R_shield, 1.05*(z_us_top-z_ls_bot)); plot_yz.pixels = (1500,2000)
plot_yz.color_by = 'material'

openmc.Plots([plot_xy, plot_yz]).export_to_xml()

# =============================================================================
# SETTINGS
# =============================================================================
settings = openmc.Settings()
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-R_outer, -R_outer, z_fuel_bot],
        [ R_outer,  R_outer, z_fuel_top],
        only_fissionable=True
    )
)
settings.batches   = 160
settings.inactive  = 40
settings.particles = 50000

entropy_mesh = openmc.RegularMesh()
entropy_mesh.lower_left  = [-R_outer, -R_outer, z_fuel_bot]
entropy_mesh.upper_right = [ R_outer,  R_outer, z_fuel_top]
entropy_mesh.dimension   = [20, 20, 20]
settings.entropy_mesh = entropy_mesh
settings.temperature  = {'default': T_struct, 'method': 'interpolation'}
settings.export_to_xml()

# =============================================================================
# TALLIES
# =============================================================================
grpstr33 = [
    1.00000E-05, 4.17458E-01, 5.31578E-01, 3.92786E+00, 8.31528E+00,
    1.37096E+01, 2.26033E+01, 3.72665E+01, 6.14421E+01, 1.01301E+02,
    1.67017E+02, 2.75364E+02, 4.53999E+02, 7.48518E+02, 1.23410E+03,
    2.03468E+03, 3.35462E+03, 5.53084E+03, 9.11881E+03, 1.50344E+04,
    2.47875E+04, 4.08677E+04, 6.73794E+04, 1.11090E+05, 1.83156E+05,
    3.01974E+05, 4.97871E+05, 8.20850E+05, 1.35335E+06, 2.23130E+06,
    3.67879E+06, 6.06531E+06, 1.00000E+07, 2.00000E+07
]

inner_ac_cell = inner_cells[2]   # 'inner_ac' axial layer
outer_ac_cell = outer_cells[2]   # 'outer_ac' axial layer

fuel_filt   = openmc.CellFilter([inner_ac_cell, outer_ac_cell])
inner_filt  = openmc.CellFilter([inner_ac_cell])
outer_filt  = openmc.CellFilter([outer_ac_cell])
all_filt    = openmc.CellFilter(inner_cells + outer_cells + [cell_refl, cell_shield])
leak_filt   = openmc.SurfaceFilter([cyl_shield, z0, z6])
neut_filt   = openmc.ParticleFilter(['neutron'])
E_filt      = openmc.EnergyFilter(grpstr33)

tallies = openmc.Tallies()

t = openmc.Tally(name='total_power')
t.filters = [fuel_filt]; t.scores = ['kappa-fission','fission','nu-fission']
tallies.append(t)

t = openmc.Tally(name='inner_power')
t.filters = [inner_filt]; t.scores = ['kappa-fission','fission']
tallies.append(t)

t = openmc.Tally(name='outer_power')
t.filters = [outer_filt]; t.scores = ['kappa-fission','fission']
tallies.append(t)

t = openmc.Tally(name='neutron_balance')
t.filters = [neut_filt, all_filt]
t.scores  = ['absorption','nu-fission','(n,2n)','(n,3n)']
tallies.append(t)

t = openmc.Tally(name='system_leakage')
t.filters = [neut_filt, leak_filt]; t.scores = ['current']
tallies.append(t)

t = openmc.Tally(name='flux_33g')
t.filters = [fuel_filt, E_filt]; t.scores = ['flux','nu-fission']
tallies.append(t)

tallies.export_to_xml()

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "="*72)
print("ABR-1000  Equilibrium Recycled Metal Core  v4 — benchmark compositions")
print("="*72)
print(f"  R_inner/outer/refl/shield  {R_inner:.1f} / {R_outer:.1f} / {R_refl:.1f} / {R_shield:.1f} cm")
print(f"  z_ls_bot → z_us_top        {z_ls_bot:.2f} → {z_us_top:.2f} cm")
print(f"  Active core                {z_fuel_bot:.2f} → {z_fuel_top:.2f} cm  ({H_active} cm)")
print(f"  T_fuel / T_struct          {T_fuel-273.15:.1f} / {T_struct-273.15:.1f} °C")
print(f"  Inner Pu239 (pin, avg)     {_I['Pu239']:.4e} a/b-cm")
print(f"  Outer Pu239 (pin, avg)     {_O['Pu239']:.4e} a/b-cm  (ratio {_O['Pu239']/_I['Pu239']:.3f})")
print(f"  Benchmark k_BOC target     1.0355 ± 0.0078  [NEA] Table 4.3")
print("="*72 + "\n")

openmc.run(threads=6)
