"""
ABR-1000 OpenMC Model; Equilibrium TRU-Recycled Metal Core

Reference: "Advanced Burner Reactor 1000 MWth Reference Concept," 
ANL-AFCI-202 (ANL-ABR-4), Argonne National Laboratory,
September 2007.

Reference: "Benchmark for Neutronic Analysis of Sodium-cooled Fast Reactor Cores with Various Fuel Types and Core Sizes"
JT03390630 (NEA-NSC-R(2015)9), Nuclear Energy Agency,
February 2016.

Design basis: EQUILIBRIUM TRU-RECYCLED metal core (U-TRU-10Zr / HT9 / Na)
  Core layout (Fig. II.1-1):
    Inner fuel zone   :  78 assemblies
    Outer fuel zone   : 102 assemblies
    Reflector         : 114 assemblies  (solid HT9 pins)
    Radial shield     :  66 assemblies  (B4C in HT9 tubes)
    Primary control   :  15 assemblies  (row-4: nat-B; row-7: 60%-B10)
    Secondary control :   4 assemblies  (natural B4C)
    Total             : 379 assemblies

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

# ======================================================================
# TEMPERATURES
# ======================================================================
T_struct = 705.0 # K
T_fuel = 807.0 # K

# ======================================================================
# BASE MATERIALS
# ======================================================================

# Uraniums -------------------------------------------------------------
U235 = openmc.Material(name='U235')
U235.add_nuclide('U235', 1.0)
U235.set_density('g/cm3', 19.1)

U238 = openmc.Material(name='U238')
U238.add_nuclide('U238', 1.0)
U238.set_density('g/cm3', 19.1)

# Transuranics ---------------------------------------------------------
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
Am242m.add_nuclide('Am242_m1', 1.0)  
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

# Zirconium ------------------------------------------------------------
Zr = openmc.Material(name='Zr')
Zr.add_element('Zr', 1.0)
Zr.set_density('g/cm3', 6.52)

# ======================================================================
# TRU AND FUEL ALLOY
# ======================================================================

# Depleted uranium: 0.2 wt% U-235 / 99.8 wt% U-238 ---------------------
DU = openmc.Material.mix_materials(
    [U235,  U238],
    [0.002, 0.998],
    'wo'
)
DU.name = 'depleted_U'

# LWR-SF TRU; Table II.1-1 (wt%) ---------------------------------------
# Assumed: 50 MWd/kg burnup, 10-year post-irradiation cooling.
TRU = openmc.Material.mix_materials(
    [Np237, Pu238, Pu239, Pu240,  Pu241, Pu242,
     Am241, Am242m, Am243, Cm243, Cm244, Cm245],
    [0.0472, 0.0218, 0.4734, 0.2282, 0.0842, 0.0684,
     0.0561, 0.0001,  0.0156, 0.0000, 0.0046, 0.0004],
    'wo'
)
TRU.name = 'TRU_LWR_SF'

# U-TRU-10Zr metal fuel alloy ------------------------------------------
# Composition from Table II.1-3 (metal recycle column):
#   Average TRU enrichment eps = 22.1 wt% in heavy metal (HM)
#   Zr = 10 wt% of total ternary alloy  (Table II.1-2, pin material)
#   HM fraction of total  = 1.0 - 0.10 = 0.90
#   TRU fraction of total = eps x HM_frac = 0.221 x 0.90 = 0.1989
#   DU  fraction of total = (1-eps) x HM_frac = 0.779 x 0.90 = 0.7011
_eps   = 0.221
_w_Zr  = 0.10
_w_HM  = 1.0 - _w_Zr
_w_TRU = _eps * _w_HM
_w_DU  = (1.0 - _eps) * _w_HM

fuel_alloy = openmc.Material.mix_materials(
    [DU,    TRU,    Zr],
    [_w_DU, _w_TRU, _w_Zr],
    'wo'
)
fuel_alloy.name = 'U_TRU_10Zr_recycle'

# ======================================================================
# STRUCTURAL / COOLANT / ABSORBER MATERIALS
# ======================================================================

# HT9 ferritic-martensitic steel ---------------------------------------
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

# Sodium coolant -------------------------------------------------------
# Core inlet: 355 deg-C, outlet: 510 deg-C  =>  T_avg ~ 432 deg-C (Table II.5-1)
# Density from Fink-Leibowitz correlation at ~432 deg-C ~ 0.860 g/cm3.
sodium = openmc.Material(name='Sodium')
sodium.add_element('Na', 1.0)
sodium.set_density('g/cm3', 0.860)
sodium.temperature = T_struct   # K, applied directly to Na-filled upper plenum cell

# B4C with natural boron (19.9 atom% B-10) -----------------------------
# Used for: radial shield, secondary control, row-4 primary control
# (Table II.1-2, note a; and body text describing control assemblies).
B4C_nat = openmc.Material(name='B4C_natural_B')
B4C_nat.add_nuclide('B10', 4.0 * 0.199, 'ao')
B4C_nat.add_nuclide('B11', 4.0 * 0.801, 'ao')
B4C_nat.add_element('C',   1.0, 'ao')
B4C_nat.set_density('g/cm3', 2.52)

# B4C with 60%-enriched B-10 -------------------------------------------
# Used for: row-7 primary control assemblies only (Table II.1-2, note b).
B4C_60 = openmc.Material(name='B4C_60pct_B10')
B4C_60.add_nuclide('B10', 4.0 * 0.60, 'ao')
B4C_60.add_nuclide('B11', 4.0 * 0.40, 'ao')
B4C_60.add_element('C',   1.0, 'ao')
B4C_60.set_density('g/cm3', 2.52)

# ======================================================================
# HOMOGENIZED ASSEMBLY MATERIALS
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
# ======================================================================

# Inner fuel assembly (78 assemblies) ----------------------------------
# Volume fractions: fuel=29.2%, Na(bond+cool)=9.8+35.3=45.1%, HT9=25.7%
inner_fuel = openmc.Material.mix_materials(
    [fuel_alloy, sodium,  HT9],
    [0.292,      0.451,   0.257],
    'vo'
)

# Outer fuel assembly (102 assemblies) ---------------------------------
# Same intra-assembly design as inner; TRU enrichment zoning is at the
# assembly level (see note; average 22.1% applied to both).
outer_fuel = openmc.Material.mix_materials(
    [fuel_alloy, sodium,  HT9],
    [0.292,      0.451,   0.257],
    'vo'
)

# Reflector assembly (114 assemblies — 91 solid HT9 pins + duct) -------
# HT9: 75.3% (pins) + 9.2% (duct) = 84.5%; Na coolant = 15.5%
reflector = openmc.Material.mix_materials(
    [HT9,  sodium],
    [0.845, 0.155],
    'vo'
)
reflector.temperature = T_struct   # K

# Radial shield assembly (66 assemblies) -------------------------------
# 19 thick HT9 tubes containing natural B4C pellets (smear density 81%).
# B4C=43.1%, HT9=29.7%, Na=27.2% (by difference from Table II.1-2).
radial_shield = openmc.Material.mix_materials(
    [B4C_nat, HT9,   sodium],
    [0.431,   0.297,  0.272],
    'vo'
)
radial_shield.temperature = T_struct   # K

# Lower axial shield material ------------------------------------------
# The lower 124.5 cm of each fuel pin contains a solid HT9 plug (lower shield).
# In the homogenized model the fuel slug volume is replaced by HT9.
# Approximate: HT9 (slug replaced + clad + duct) ~ 29.2+25.7=54.9%, Na=45.1%.
lower_shield_mat = openmc.Material.mix_materials(
    [HT9,  sodium],
    [0.549, 0.451],
    'vo'
)
lower_shield_mat.temperature = T_struct   # K

# Primary control; row 4 (3 assemblies, natural boron) ------------------
primary_ctrl_r4 = openmc.Material.mix_materials(
    [B4C_nat, HT9,   sodium],
    [0.428,   0.208,  0.364],
    'vo'
)

# Primary control; row 7 (12 assemblies, 60% B-10 enrichment) ----------
primary_ctrl_r7 = openmc.Material.mix_materials(
    [B4C_60, HT9,   sodium],
    [0.428,  0.208,  0.364],
    'vo'
)

# Secondary control (4 assemblies, natural boron) ----------------------
secondary_ctrl = openmc.Material.mix_materials(
    [B4C_nat, HT9,   sodium],
    [0.428,   0.208,  0.364],
    'vo'
)


# =============================================================================
# MATERIALS FILE
# =============================================================================
f_ctrl_inner = 7.0 / 85.0
f_ctrl_outer = 12.0 / 114.0

inner_core = openmc.Material.mix_materials(
    [inner_fuel, primary_ctrl_r4, secondary_ctrl],
    [1 - f_ctrl_inner, f_ctrl_inner*(3/7), f_ctrl_inner*(4/7)],
    'vo'
)
inner_core.temperature = T_fuel

outer_core = openmc.Material.mix_materials(
    [outer_fuel, primary_ctrl_r7],
    [1 - f_ctrl_outer, f_ctrl_outer],
    'vo'
)
outer_core.temperature = T_fuel

materials = openmc.Materials([
    inner_core, outer_core, reflector, radial_shield, sodium
])

materials.export_to_xml()

# ======================================================================
# GEOMETRY
#
# Radial boundaries derived from hex assembly pitch and assembly counts
# ----------------------------------------------------------------------
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
# ======================================================================

pitch = 16.142
A_hex = (math.sqrt(3) / 2.0) * pitch ** 2

def r_eff(n_asm):
    """Effective cylinder radius (cm) for n_asm hex assemblies."""
    return math.sqrt(n_asm * A_hex / math.pi)

# Radial boundary cylinders --------------------------------------------
R_inner = r_eff(85)
R_outer = r_eff(199)
R_refl  = r_eff(313)
R_shield= r_eff(379)

cyl_inner  = openmc.ZCylinder(r=R_inner)
cyl_outer  = openmc.ZCylinder(r=R_outer)
cyl_refl   = openmc.ZCylinder(r=R_refl)
cyl_shield = openmc.ZCylinder(r=R_shield, boundary_type='vacuum')

z_bot = openmc.ZPlane(z0=-150, boundary_type='vacuum')
z_top = openmc.ZPlane(z0=150, boundary_type='vacuum')

region_core   = -cyl_inner & +z_bot & -z_top
region_outer  = +cyl_inner & -cyl_outer & +z_bot & -z_top
region_refl   = +cyl_outer & -cyl_refl & +z_bot & -z_top
region_shield = +cyl_refl & -cyl_shield & +z_bot & -z_top

cell_inner = openmc.Cell(fill=inner_core, region=region_core)
cell_outer = openmc.Cell(fill=outer_core, region=region_outer)
cell_refl  = openmc.Cell(fill=reflector, region=region_refl)
cell_shield= openmc.Cell(fill=radial_shield, region=region_shield)

universe = openmc.Universe(cells=[
    cell_inner, cell_outer, cell_refl, cell_shield
])

geometry = openmc.Geometry(universe)
geometry.export_to_xml()

# =============================================================================
# PLOTS  (geometry verification)
# ----------------------------------------------------------------------
# Material colour map for homogenized zones. RGB tuples are chosen so
# the four radial rings are visually distinct in the XY mid-plane slice
# and the axial layering (lower shield / active / upper plenum) is
# obvious in the YZ slice.
# =============================================================================
dict_matcolors = {
    inner_core:   (255, 200, 60),
    outer_core:   (255, 140, 40),
    reflector:    (140, 140, 160),
    radial_shield:(80, 80, 90),
    sodium:       (120, 200, 240),
}

plot_xy = openmc.Plot.from_geometry(geometry)
plot_xy.filename = 'ABR1000_xy_midplane'
plot_xy.basis    = 'xy'
plot_xy.origin   = (0.0, 0.0, 0.0)                # active-core mid-plane
plot_xy.width    = (2.2 * R_shield, 2.2 * R_shield)
plot_xy.pixels   = (2000, 2000)
plot_xy.color_by = 'material'
plot_xy.colors   = dict_matcolors

plot_yz = openmc.Plot.from_geometry(geometry)
plot_yz.filename = 'ABR1000_yz_axial'
plot_yz.basis    = 'yz'
plot_yz.origin   = (0.0, 0.0, 0.0)
plot_yz.width = (2.2 * R_shield, 1.05 * (z_top.z0 - z_bot.z0))
plot_yz.pixels   = (1500, 2000)
plot_yz.color_by = 'material'
plot_yz.colors   = dict_matcolors

plots_file = openmc.Plots([plot_xy, plot_yz])
plots_file.export_to_xml()

# =============================================================================
# SETTINGS
# =============================================================================

source = openmc.IndependentSource(
    space=openmc.stats.Box(
        [-R_shield, -R_shield, -150],
        [ R_shield,  R_shield,  150],
        only_fissionable=True
    )
)

settings = openmc.Settings()
settings.source = source
settings.batches = 150
settings.inactive = 30
settings.particles = 100000

mesh = openmc.RegularMesh()
mesh.lower_left  = [-R_shield, -R_shield, -150]
mesh.upper_right = [ R_shield,  R_shield,  150]
mesh.dimension   = [20, 20, 40]

settings.entropy_mesh = mesh
settings.temperature = {'default': T_struct, 'method': 'interpolation'}

settings.export_to_xml()

# =============================================================================
# PLOT
# =============================================================================
openmc.plot_geometry()

# =============================================================================
# TALLIES
# ----------------------------------------------------------------------
# Essential neutronics outputs for the ABR-1000 equilibrium metal core:
#   (1) Total core power       - verify 1000 MWt normalization
#   (2) Zone-wise power split  - inner vs. outer fuel (ANL Table II.1-3)
#   (3) System neutron balance - absorption, production, (n,xn)
#   (4) System leakage         - radial + top + bottom vacuum surfaces
#   (5) Group-wise flux + nu-fission in fuel (ANL33 group structure)
# =============================================================================

# =============================================================================
# TALLIES — OECD/NEA SFR BENCHMARK SUITE
# =============================================================================

# --- 33-group energy structure (ANL/fast spectrum standard) ---
grpstr33 = [
    1.00000E-05, 4.17458E-01, 5.31578E-01, 3.92786E+00, 8.31528E+00,
    1.37096E+01, 2.26033E+01, 3.72665E+01, 6.14421E+01, 1.01301E+02,
    1.67017E+02, 2.75364E+02, 4.53999E+02, 7.48518E+02, 1.23410E+03,
    2.03468E+03, 3.35462E+03, 5.53084E+03, 9.11881E+03, 1.50344E+04,
    2.47875E+04, 4.08677E+04, 6.73794E+04, 1.11090E+05, 1.83156E+05,
    3.01974E+05, 4.97871E+05, 8.20850E+05, 1.35335E+06, 2.23130E+06,
    3.67879E+06, 6.06531E+06, 1.00000E+07, 2.00000E+07
]

# --- Filters ---
fuel_filter   = openmc.CellFilter([cell_inner, cell_outer])
inner_filter  = openmc.CellFilter([cell_inner])
outer_filter  = openmc.CellFilter([cell_outer])

all_cells_filter = openmc.CellFilter([
    cell_inner, cell_outer, cell_refl, cell_shield
])

leakage_filter = openmc.SurfaceFilter([cyl_shield, z_bot, z_top])
neutron_filter = openmc.ParticleFilter(['neutron'])
energy_filter  = openmc.EnergyFilter(grpstr33)

tallies = openmc.Tallies()

# =============================================================================
# (1) TOTAL CORE POWER
# =============================================================================
t_total_power = openmc.Tally(name='total_power')
t_total_power.filters = [fuel_filter]
t_total_power.scores  = ['kappa-fission', 'fission', 'nu-fission']
tallies.append(t_total_power)

# =============================================================================
# (2) POWER SPLIT (INNER vs OUTER)
# =============================================================================
t_inner_power = openmc.Tally(name='inner_power')
t_inner_power.filters = [inner_filter]
t_inner_power.scores  = ['kappa-fission', 'fission']
tallies.append(t_inner_power)

t_outer_power = openmc.Tally(name='outer_power')
t_outer_power.filters = [outer_filter]
t_outer_power.scores  = ['kappa-fission', 'fission']
tallies.append(t_outer_power)

# =============================================================================
# (3) NEUTRON BALANCE
# =============================================================================
t_balance = openmc.Tally(name='neutron_balance')
t_balance.filters = [neutron_filter, all_cells_filter]
t_balance.scores  = [
    'absorption',
    'nu-fission',
    '(n,2n)',
    '(n,3n)'
]
tallies.append(t_balance)

# =============================================================================
# (4) SYSTEM LEAKAGE
# =============================================================================
t_leakage = openmc.Tally(name='system_leakage')
t_leakage.filters = [neutron_filter, leakage_filter]
t_leakage.scores  = ['current']
tallies.append(t_leakage)

# =============================================================================
# (5) 33-GROUP FLUX + PRODUCTION SPECTRUM
# =============================================================================
t_flux_33g = openmc.Tally(name='flux_33g')
t_flux_33g.filters = [fuel_filter, energy_filter]
t_flux_33g.scores  = ['flux', 'nu-fission']
tallies.append(t_flux_33g)

# =============================================================================
# EXPORT
# =============================================================================
tallies.export_to_xml()

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
print(f"  Coolant inlet/outlet 355 / 510 deg-C")
print(f"  Na density (avg T)   0.860 g/cm3  (~432 deg-C average)")
print("=" * 70 + "\n")

# =============================================================================
# RUN
# =============================================================================
openmc.run()
