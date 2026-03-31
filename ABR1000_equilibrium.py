"""
ABR-1000 OpenMC Model  —  Equilibrium TRU-Recycled Metal Core
=============================================================
Reference: Grandy & Seidensticker (eds.), "Advanced Burner Reactor 1000 MWth
Reference Concept," ANL-AFCI-202 (ANL-ABR-4), Argonne National Laboratory,
September 2007.

Design basis: EQUILIBRIUM TRU-RECYCLED metal core (U-TRU-10Zr / HT9 / Na)
  Core layout (Fig. II.1-1):
    Inner fuel zone   —  78 assemblies
    Outer fuel zone   — 102 assemblies
    Reflector         — 114 assemblies  (solid HT9 pins)
    Radial shield     —  66 assemblies  (B4C in HT9 tubes)
    Primary control   —  15 assemblies  (row-4: nat-B; row-7: 60%-B10)
    Secondary control —   4 assemblies  (natural B4C)
    Total             — 379 assemblies

  Key equilibrium-recycle parameters (Table II.1-3, metal recycle column):
    TRU feed             : LWR-SF (50 MWd/kg, 10-yr cooling)
    Avg TRU enrichment   : 22.1 wt% in heavy metal
    Number of batches    : 4
    HM inventory at BOEC : 13.2 MT
    Avg discharge burnup : 93 MWd/kg
    Power density        : 303 kW/L
    Thermal power        : 1000 MWt
    Coolant inlet/outlet : 355 deg-C / 510 deg-C

Geometry: homogenized cylindrical approximation of the hex assembly lattice.
  Effective radii are derived from hex cell area and assembly counts.
  A full pin-by-pin hex-lattice model requires explicit universe definitions.

  Per-zone TRU enrichment note: the report gives only the CORE-AVERAGE of
  22.1 wt% (Table II.1-3). The zone-specific split is determined iteratively
  in the ANL neutronics code but is not tabulated; both zones use 22.1% here.
"""

# %matplotlib inline
import math
import openmc

# =============================================================================
# SECTION 1 — BASE NUCLIDE / ELEMENT MATERIALS
# =============================================================================

# -- Uranium isotopes ----------------------------------------------------------
U235 = openmc.Material(name='U235')
U235.add_nuclide('U235', 1.0)
U235.set_density('g/cm3', 19.1)

U238 = openmc.Material(name='U238')
U238.add_nuclide('U238', 1.0)
U238.set_density('g/cm3', 19.1)

# -- Transuranics (individual nuclides) ----------------------------------------
Np237 = openmc.Material(name='Np237')
Np237.add_nuclide('Np237', 1.0)
Np237.set_density('g/cm3', 20.5)

Pu238 = openmc.Material(name='Pu238')
Pu238.add_nuclide('Pu238', 1.0)
Pu238.set_density('g/cm3', 19.8)

Pu239 = openmc.Material(name='Pu239')
Pu239.add_nuclide('Pu239', 1.0)
Pu239.set_density('g/cm3', 19.8)

Pu240 = openmc.Material(name='Pu240')
Pu240.add_nuclide('Pu240', 1.0)
Pu240.set_density('g/cm3', 19.8)

Pu241 = openmc.Material(name='Pu241')
Pu241.add_nuclide('Pu241', 1.0)
Pu241.set_density('g/cm3', 19.8)

Pu242 = openmc.Material(name='Pu242')
Pu242.add_nuclide('Pu242', 1.0)
Pu242.set_density('g/cm3', 19.8)

Am241 = openmc.Material(name='Am241')
Am241.add_nuclide('Am241', 1.0)
Am241.set_density('g/cm3', 13.7)

Am242m = openmc.Material(name='Am242m')
Am242m.add_nuclide('Am242_m1', 1.0)  # metastable — correct OpenMC nuclide tag
Am242m.set_density('g/cm3', 13.7)

Am243 = openmc.Material(name='Am243')
Am243.add_nuclide('Am243', 1.0)
Am243.set_density('g/cm3', 13.7)

Cm243 = openmc.Material(name='Cm243')
Cm243.add_nuclide('Cm243', 1.0)
Cm243.set_density('g/cm3', 13.5)

Cm244 = openmc.Material(name='Cm244')
Cm244.add_nuclide('Cm244', 1.0)
Cm244.set_density('g/cm3', 13.5)

Cm245 = openmc.Material(name='Cm245')
Cm245.add_nuclide('Cm245', 1.0)
Cm245.set_density('g/cm3', 13.5)

# -- Zirconium alloying element ------------------------------------------------
Zr = openmc.Material(name='Zr')
Zr.add_element('Zr', 1.0)
Zr.set_density('g/cm3', 6.52)

# =============================================================================
# SECTION 2 — TRU VECTOR AND FUEL ALLOY
# =============================================================================

# -- Depleted uranium: 0.2 wt% U-235 / 99.8 wt% U-238 ------------------------
DU = openmc.Material.mix_materials(
    [U235,  U238],
    [0.002, 0.998],
    'wo'
)
DU.name = 'depleted_U'

# -- LWR-SF TRU vector — Table II.1-1 (wt%) -----------------------------------
# Assumed: 50 MWd/kg burnup, 10-year post-irradiation cooling.
# This is the TRU feed for the EQUILIBRIUM RECYCLED core.
TRU = openmc.Material.mix_materials(
    [Np237, Pu238, Pu239, Pu240,  Pu241, Pu242,
     Am241, Am242m, Am243, Cm243, Cm244, Cm245],
    [0.0472, 0.0218, 0.4734, 0.2282, 0.0842, 0.0684,
     0.0561, 0.0001,  0.0156, 0.0000, 0.0046, 0.0004],
    'wo'
)
TRU.name = 'TRU_LWR_SF'

# -- U-TRU-10Zr metal fuel alloy -----------------------------------------------
# Composition from Table II.1-3 (metal recycle column):
#   Average TRU enrichment eps = 22.1 wt% in heavy metal (HM)
#   Zr = 10 wt% of total ternary alloy  (Table II.1-2, pin material)
#   HM fraction of total  = 1.0 - 0.10 = 0.90
#   TRU fraction of total = eps x HM_frac = 0.221 x 0.90 = 0.1989
#   DU  fraction of total = (1-eps) x HM_frac = 0.779 x 0.90 = 0.7011
_eps   = 0.221             # TRU/HM wt fraction (Table II.1-3 recycle average)
_w_Zr  = 0.10              # Zr wt fraction in ternary alloy (Table II.1-2)
_w_HM  = 1.0 - _w_Zr      # heavy metal fraction = 0.90
_w_TRU = _eps * _w_HM     # = 0.1989
_w_DU  = (1.0 - _eps) * _w_HM  # = 0.7011

fuel_alloy = openmc.Material.mix_materials(
    [DU,    TRU,    Zr],
    [_w_DU, _w_TRU, _w_Zr],
    'wo'
)
fuel_alloy.name = 'U_TRU_10Zr_recycle'

# =============================================================================
# SECTION 3 — STRUCTURAL / COOLANT / ABSORBER MATERIALS
# =============================================================================

# -- HT9 ferritic-martensitic steel (cladding AND duct for all assembly types) -
# Used as cladding, wire wrap, and hexagonal duct material (Table II.1-2).
HT9 = openmc.Material(name='HT9')
HT9.add_element('Fe', 0.854, 'wo')
HT9.add_element('Cr', 0.120, 'wo')
HT9.add_element('Mo', 0.010, 'wo')
HT9.add_element('Mn', 0.006, 'wo')
HT9.add_element('Si', 0.004, 'wo')
HT9.add_element('C',  0.002, 'wo')
HT9.add_element('Ni', 0.004, 'wo')
HT9.set_density('g/cm3', 7.70)

# -- Sodium coolant at average operating temperature ---------------------------
# Core inlet: 355 deg-C, outlet: 510 deg-C  =>  T_avg ~ 432 deg-C (Table II.5-1)
# Density from Fink-Leibowitz correlation at ~432 deg-C ~ 0.860 g/cm3.
sodium = openmc.Material(name='Sodium')
sodium.add_element('Na', 1.0)
sodium.set_density('g/cm3', 0.860)

# -- B4C with NATURAL boron (19.9 atom% B-10) ----------------------------------
# Used for: radial shield, secondary control, row-4 primary control
# (Table II.1-2, note a; and body text describing control assemblies).
B4C_nat = openmc.Material(name='B4C_natural_B')
B4C_nat.add_nuclide('B10', 4.0 * 0.199, 'ao')  # 19.9 at% B-10 in natural B
B4C_nat.add_nuclide('B11', 4.0 * 0.801, 'ao')
B4C_nat.add_element('C',   1.0, 'ao')
B4C_nat.set_density('g/cm3', 2.52)

# -- B4C with 60%-enriched B-10 ------------------------------------------------
# Used for: row-7 primary control assemblies only (Table II.1-2, note b).
B4C_60 = openmc.Material(name='B4C_60pct_B10')
B4C_60.add_nuclide('B10', 4.0 * 0.60, 'ao')
B4C_60.add_nuclide('B11', 4.0 * 0.40, 'ao')
B4C_60.add_element('C',   1.0, 'ao')
B4C_60.set_density('g/cm3', 2.52)

# =============================================================================
# SECTION 4 — HOMOGENIZED ASSEMBLY MATERIALS
#
# Volume fractions at fabrication from Table II.1-2 (metal fuel core):
#
#  Assembly type | Fuel/Abs | Bond Na | Structure | Coolant Na
#  -------------+----------+---------+-----------+-----------
#  Fuel (metal) |  29.2 %  |  9.8 %  |   25.7 %  |  35.3 %
#  Reflector    |   ---    |  ---    |   84.5 %  |  15.5 %
#  Shield       |  43.1 %  |  ---    |   29.7 %  |  27.2 %  (by difference)
#  Control      |  42.8 %  |  7.6 %* |   20.8 %  |  28.8 %
#
#  * Control pin bond is He (negligible density). The 7.6% He bond volume is
#    merged with the 28.8% Na coolant => total Na = 36.4% for the mix call.
#
# Mixing mode 'vo' (volume fraction): rho_mix = SUM(vf_i * rho_i).
# =============================================================================

# -- Inner fuel assembly (78 assemblies) ---------------------------------------
# Volume fractions: fuel=29.2%, Na(bond+cool)=9.8+35.3=45.1%, HT9=25.7%
inner_fuel = openmc.Material.mix_materials(
    [fuel_alloy, sodium,  HT9],
    [0.292,      0.451,   0.257],
    'vo'
)
inner_fuel.name = 'inner_fuel_asm'

# -- Outer fuel assembly (102 assemblies) --------------------------------------
# Same intra-assembly design as inner; TRU enrichment zoning is at the
# assembly level (see module docstring note — average 22.1% applied to both).
outer_fuel = openmc.Material.mix_materials(
    [fuel_alloy, sodium,  HT9],
    [0.292,      0.451,   0.257],
    'vo'
)
outer_fuel.name = 'outer_fuel_asm'

# -- Reflector assembly (114 assemblies — 91 solid HT9 pins + duct) ------------
# HT9: 75.3% (pins) + 9.2% (duct) = 84.5%; Na coolant = 15.5%
reflector = openmc.Material.mix_materials(
    [HT9,  sodium],
    [0.845, 0.155],
    'vo'
)
reflector.name = 'reflector_asm'

# -- Radial shield assembly (66 assemblies) ------------------------------------
# 19 thick HT9 tubes containing natural B4C pellets (smear density 81%).
# B4C=43.1%, HT9=29.7%, Na=27.2% (by difference from Table II.1-2).
radial_shield = openmc.Material.mix_materials(
    [B4C_nat, HT9,   sodium],
    [0.431,   0.297,  0.272],
    'vo'
)
radial_shield.name = 'radial_shield_asm'

# -- Lower axial shield material (below active core, fuel-assembly footprint) --
# The lower 124.5 cm of each fuel pin contains a solid HT9 plug (lower shield).
# In the homogenized model the fuel slug volume is replaced by HT9.
# Approximate: HT9 (slug replaced + clad + duct) ~ 29.2+25.7=54.9%, Na=45.1%.
lower_shield_mat = openmc.Material.mix_materials(
    [HT9,  sodium],
    [0.549, 0.451],
    'vo'
)
lower_shield_mat.name = 'lower_axial_shield'

# -- Primary control — row 4 (3 assemblies, natural boron) ---------------------
primary_ctrl_r4 = openmc.Material.mix_materials(
    [B4C_nat, HT9,   sodium],
    [0.428,   0.208,  0.364],
    'vo'
)
primary_ctrl_r4.name = 'primary_ctrl_row4_natB'

# -- Primary control — row 7 (12 assemblies, 60% B-10 enrichment) --------------
primary_ctrl_r7 = openmc.Material.mix_materials(
    [B4C_60, HT9,   sodium],
    [0.428,  0.208,  0.364],
    'vo'
)
primary_ctrl_r7.name = 'primary_ctrl_row7_60B10'

# -- Secondary control (4 assemblies, natural boron) ---------------------------
secondary_ctrl = openmc.Material.mix_materials(
    [B4C_nat, HT9,   sodium],
    [0.428,   0.208,  0.364],
    'vo'
)
secondary_ctrl.name = 'secondary_ctrl_natB'

# =============================================================================
# MATERIALS FILE
# =============================================================================
materials_file = openmc.Materials([
    inner_fuel, outer_fuel,
    reflector, radial_shield, lower_shield_mat,
    primary_ctrl_r4, primary_ctrl_r7, secondary_ctrl,
    sodium
])
materials_file.export_to_xml()

# =============================================================================
# SECTION 5 — GEOMETRY
#
# Radial boundaries derived from hex assembly pitch and assembly counts
# ---------------------------------------------------------------------------
# Hex cell area:  A_hex = (sqrt(3)/2) x pitch^2
#   pitch = 16.142 cm  =>  A_hex = 225.655 cm^2
#
# Equivalent cylinder radius for N assemblies:  r = sqrt(N x A_hex / pi)
#
# Grouping (Fig. II.1-1 + Table II.1-2):
#   Inner zone  : 78 fuel + 4 sec-ctrl + 3 pri-ctrl(row-4)  =  85 => r ~ 78.1 cm
#   + Outer zone: 102 fuel + 12 pri-ctrl(row-7)             = 199 => r ~119.6 cm
#   + Reflector : 114 assemblies                            = 313 => r ~149.9 cm
#   + Rad shield:  66 assemblies                            = 379 => r ~165.0 cm
#
# Power density self-check:
#   V_active = 180 fuel asm x 225.655 cm^2 x 81.3 cm = 3 302 200 cm^3 = 3302 L
#   P_density = 1000 MWt / 3302 L = 302.8 kW/L   (ANL Table II.1-3: 303 kW/L) ✓
#
# Axial structure (Table II.1-2, metal fuel pin dimensions, centred at z=0):
#   Lower HT9 shield pin : 124.5 cm   (z = -165.15 to -40.65)
#   Active fuel zone     :  81.3 cm   (z =  -40.65 to +40.65)
#   Fission gas plenum   : 124.5 cm   (z =  +40.65 to +165.15)
# =============================================================================

pitch = 16.142                                 # cm (Table II.1-2)
A_hex = (math.sqrt(3) / 2.0) * pitch ** 2    # = 225.655 cm^2 per hex cell

def r_eff(n_asm):
    """Effective cylinder radius (cm) for n_asm hex assemblies."""
    return math.sqrt(n_asm * A_hex / math.pi)

# -- Radial boundary cylinders -------------------------------------------------
R_inner  = r_eff(85)    # ~  78.1 cm  inner zone + sec/row-4 pri ctrl
R_outer  = r_eff(199)   # ~ 119.6 cm  + outer zone + row-7 pri ctrl
R_refl   = r_eff(313)   # ~ 149.9 cm  + reflector assemblies
R_shield = r_eff(379)   # ~ 165.0 cm  + radial shield assemblies

cyl_inner  = openmc.ZCylinder(r=R_inner)
cyl_outer  = openmc.ZCylinder(r=R_outer)
cyl_refl   = openmc.ZCylinder(r=R_refl)
cyl_shield = openmc.ZCylinder(r=R_shield, boundary_type='vacuum')

# -- Axial boundary planes -----------------------------------------------------
H_active  = 81.3    # cm — active fuel height          (Table II.1-2)
H_lower   = 124.5   # cm — lower HT9 shield pin length (Table II.1-2)
H_plenum  = 124.5   # cm — fission gas plenum height   (Table II.1-2)

z_bot_active = -H_active / 2.0              # -40.65 cm
z_top_active = +H_active / 2.0              # +40.65 cm
z_bot        =  z_bot_active - H_lower      # -165.15 cm
z_top        =  z_top_active + H_plenum     # +165.15 cm

plane_bot     = openmc.ZPlane(z0=z_bot,        boundary_type='vacuum')
plane_bot_act = openmc.ZPlane(z0=z_bot_active)
plane_top_act = openmc.ZPlane(z0=z_top_active)
plane_top     = openmc.ZPlane(z0=z_top,        boundary_type='vacuum')

# Axial slab shorthand
active_slab = +plane_bot_act & -plane_top_act
lower_slab  = +plane_bot     & -plane_bot_act
upper_slab  = +plane_top_act & -plane_top

# ==================== ACTIVE CORE AXIAL SLICE ====================

cell_inner_act = openmc.Cell(
    name   = 'inner_core_active',
    fill   = inner_fuel,
    region = -cyl_inner & active_slab
)

cell_outer_act = openmc.Cell(
    name   = 'outer_core_active',
    fill   = outer_fuel,
    region = +cyl_inner & -cyl_outer & active_slab
)

cell_refl_act = openmc.Cell(
    name   = 'radial_reflector_active',
    fill   = reflector,
    region = +cyl_outer & -cyl_refl & active_slab
)

cell_shield_act = openmc.Cell(
    name   = 'radial_shield_active',
    fill   = radial_shield,
    region = +cyl_refl & -cyl_shield & active_slab
)

# ============= LOWER AXIAL ZONE (HT9 lower shield pins) =============
# Within fuel-assembly radius: lower HT9 shield plugs in sodium flow.
# Radially beyond: reflector/shield duct structures continue downward.

cell_lower_fuel_zone = openmc.Cell(
    name   = 'lower_axial_fuel_zone',
    fill   = lower_shield_mat,
    region = -cyl_outer & lower_slab
)

cell_lower_outer_ring = openmc.Cell(
    name   = 'lower_axial_outer_ring',
    fill   = reflector,       # reflector/shield duct material
    region = +cyl_outer & -cyl_shield & lower_slab
)

# =========== UPPER AXIAL ZONE (fission gas plenum + sodium) ===========
# Gas plenum is sealed inside pins; above the core the open sodium
# pool dominates. Modelled as pure sodium flowing upward.

cell_upper = openmc.Cell(
    name   = 'upper_axial_plenum',
    fill   = sodium,
    region = -cyl_shield & upper_slab
)

# -- Assemble universe and geometry --------------------------------------------
main_universe = openmc.Universe(cells=[
    cell_inner_act, cell_outer_act,
    cell_refl_act,  cell_shield_act,
    cell_lower_fuel_zone, cell_lower_outer_ring,
    cell_upper
])

geometry = openmc.Geometry(main_universe)
geometry.export_to_xml()

# =============================================================================
# SECTION 6 — SETTINGS
# =============================================================================
# Restrict fission source to the active fissile core region.
src_box = openmc.stats.Box(
    lower_left  = [-R_outer, -R_outer, z_bot_active],
    upper_right = [ R_outer,  R_outer, z_top_active]
)
source = openmc.IndependentSource(
    space       = src_box,
    constraints = {'fissionable': True}
)

settings = openmc.Settings()
settings.source    = source
settings.batches   = 110     # 10 inactive + 100 active
settings.inactive  = 10
settings.particles = 5000    # increase to >=50 000 for production k-eff runs
settings.export_to_xml()

# =============================================================================
# SECTION 7 — TALLIES  (uncomment as needed)
# =============================================================================
# Flux and fission rate tally across active fuel zones:
#
# flux_tally = openmc.Tally(name='core_flux_fission')
# flux_tally.scores  = ['flux', 'fission', 'nu-fission']
# flux_tally.filters = [openmc.CellFilter([cell_inner_act, cell_outer_act])]
# tallies = openmc.Tallies([flux_tally])
# tallies.export_to_xml()

# =============================================================================
# SECTION 8 — RUN
# =============================================================================
openmc.run()

# =============================================================================
# DESIGN PARAMETER SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("ABR-1000  Equilibrium Recycled Metal Core  (ANL-AFCI-202)")
print("=" * 70)
print(f"  Assembly pitch                    {pitch:.3f} cm   (Table II.1-2)")
print(f"  Hex cell area                     {A_hex:.2f} cm^2")
print(f"  Pins per fuel assembly            271             (Table II.1-2)")
print(f"  Fuel pin OD / clad thickness      0.755 / 0.056 cm")
print(f"  Fuel slug OD (75% smear)          0.557 cm")
print(f"  Wire wrap diameter                0.131 cm")
print(f"  P/D ratio                         1.180")
print(f"  Active core height                {H_active:.1f} cm")
print(f"  Lower HT9 shield length           {H_lower:.1f} cm")
print(f"  Fission gas plenum length         {H_plenum:.1f} cm")
print(f"  Fuel smear density                75 % TD")
print(f"  TRU enrichment (avg, recycle)     22.1 wt% in HM  (Table II.1-3)")
print(f"  Zr content                        10.0 wt% of alloy")
print(f"  TRU feed                          LWR-SF (Table II.1-1)")
print()
print(f"  R_inner  ( 85 asm)   {R_inner:.2f} cm")
print(f"  R_outer  (199 asm)   {R_outer:.2f} cm")
print(f"  R_refl   (313 asm)   {R_refl:.2f} cm")
print(f"  R_shield (379 asm)   {R_shield:.2f} cm")
print()
V_act = 180 * A_hex * H_active
print(f"  Active core volume   {V_act/1000:.1f} L")
print(f"  Power density check  {1e6/V_act:.1f} kW/L  (ANL: 303 kW/L)")
print(f"  Coolant inlet/outlet 355 / 510 deg-C")
print(f"  Na density (avg T)   0.860 g/cm3  (~432 deg-C average)")
print("=" * 70 + "\n")
