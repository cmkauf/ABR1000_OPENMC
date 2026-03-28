%matplotlib inline
import openmc

# =============================================================================
# MATERIALS — Pure nuclides / elements
# =============================================================================

# --- Uranium isotopes ---
U235 = openmc.Material(name='U235')
U235.add_nuclide('U235', 1.0)
U235.set_density('g/cm3', 19.1)

U238 = openmc.Material(name='U238')
U238.add_nuclide('U238', 1.0)
U238.set_density('g/cm3', 19.1)

# --- Transuranics ---
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

# FIX: metastable Am-242 uses nuclide name 'Am242_m1', not 'Am242'
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

# --- Other elements ---
Zr = openmc.Material(name='Zr')
Zr.add_element('Zr', 1.0)
Zr.set_density('g/cm3', 6.52)   # FIX: correct density for Zr metal is ~6.52, not 13.5

# =============================================================================
# MATERIAL MIXTURES
# =============================================================================

# Depleted uranium (0.2 wt% U235, 99.8 wt% U238) — used as diluent in fuel
# FIX: original U mix had empty lists
U = openmc.Material.mix_materials(
    [U235,  U238],
    [0.002, 0.998],
    'wo'
)

# TRU vector — recycled core (weight fractions must sum to 1.0)
TRU = openmc.Material.mix_materials(
    [Np237, Pu238, Pu239, Pu240, Pu241, Pu242, Am241, Am242m, Am243, Cm243, Cm244, Cm245],
    [0.0472, 0.0218, 0.4734, 0.2282, 0.0842, 0.0684, 0.0561, 0.0001, 0.0156, 0.0000, 0.0046, 0.0004],
    'wo'
)

# Weapons-grade Pu — startup core
# FIX 1: two missing commas between arguments
# FIX 2: fractions were percentages; divide by 100 so they sum to 1.0
WGPu = openmc.Material.mix_materials(
    [Pu238, Pu239, Pu240,  Pu241, Pu242],
    [0.0001, 0.9381, 0.0581, 0.0035, 0.0002],
    'wo'
)

# --- Fuel zone compositions (U-TRU-Zr metal fuel) ---
# Inner core: 15 wt% TRU, 10 wt% Zr, 75 wt% depleted U
inner_fuel = openmc.Material.mix_materials(
    [U,    TRU,  Zr],
    [0.75, 0.15, 0.10],
    'wo'
)
inner_fuel.name = 'inner_fuel'

# Outer core: 20 wt% TRU, 10 wt% Zr, 70 wt% depleted U
outer_fuel = openmc.Material.mix_materials(
    [U,    TRU,  Zr],
    [0.70, 0.20, 0.10],
    'wo'
)
outer_fuel.name = 'outer_fuel'

# --- Structural / coolant materials ---

# HT-9 ferritic-martensitic steel cladding (Fe-12Cr-1Mo)
cladding = openmc.Material(name='HT9_cladding')
cladding.add_element('Fe', 0.854, 'wo')
cladding.add_element('Cr', 0.120, 'wo')
cladding.add_element('Mo', 0.010, 'wo')
cladding.add_element('Mn', 0.006, 'wo')
cladding.add_element('Si', 0.004, 'wo')
cladding.add_element('C',  0.002, 'wo')
cladding.add_element('Ni', 0.004, 'wo')
cladding.set_density('g/cm3', 7.7)

# Sodium coolant (density at ~500 °C operating temperature)
sodium = openmc.Material(name='Sodium')
sodium.add_element('Na', 1.0)
sodium.set_density('g/cm3', 0.850)

# Reflector: depleted uranium (same as U mix above, but we clone for clarity)
reflector = openmc.Material.mix_materials(
    [U235,  U238],
    [0.002, 0.998],
    'wo'
)
reflector.name = 'reflector'

# Radiation shield: B4C + SS316 (50/50 wt%) — simplified
# B4C component
B4C = openmc.Material(name='B4C')
B4C.add_element('B', 4.0, 'ao')
B4C.add_element('C', 1.0, 'ao')
B4C.set_density('g/cm3', 2.52)

# SS316 for shield matrix
SS316 = openmc.Material(name='SS316')
SS316.add_element('Fe', 0.655, 'wo')
SS316.add_element('Cr', 0.170, 'wo')
SS316.add_element('Ni', 0.120, 'wo')
SS316.add_element('Mo', 0.025, 'wo')
SS316.add_element('Mn', 0.020, 'wo')
SS316.add_element('Si', 0.010, 'wo')
SS316.set_density('g/cm3', 8.0)

shield = openmc.Material.mix_materials(
    [B4C,  SS316],
    [0.50, 0.50],
    'wo'
)
shield.name = 'shield'

# Primary control: 90 wt% B4C (enriched 90% B-10), 10 wt% SS316
B4C_enriched = openmc.Material(name='B4C_enriched')
B4C_enriched.add_nuclide('B10', 0.9 * 4, 'ao')   # 90% enriched in B-10
B4C_enriched.add_nuclide('B11', 0.1 * 4, 'ao')
B4C_enriched.add_element('C', 1.0, 'ao')
B4C_enriched.set_density('g/cm3', 2.52)

primary_control = openmc.Material.mix_materials(
    [B4C_enriched, SS316],
    [0.90,         0.10],
    'wo'
)
primary_control.name = 'primary_control'

# Secondary control: natural B4C + SS316 (less enriched)
secondary_control = openmc.Material.mix_materials(
    [B4C,  SS316],
    [0.80, 0.20],
    'wo'
)
secondary_control.name = 'secondary_control'

# =============================================================================
# MATERIALS FILE
# FIX: original referenced undefined 'sodium' and 'clad' (should be 'cladding')
# =============================================================================
materials_file = openmc.Materials([
    inner_fuel, outer_fuel, cladding, sodium,
    reflector, shield, primary_control, secondary_control
])
materials_file.export_to_xml()

# =============================================================================
# GEOMETRY — Simplified cylindrical SFR core
#
# Radial zones (outward):
#   Inner core fuel  | Outer core fuel  | Reflector  | Shield
# Axial zones (upward):
#   Lower shield | Lower reflector | Active fuel | Upper reflector | Upper shield
#
# For a full hexagonal pin/assembly lattice see OpenMC hex-lattice examples.
# =============================================================================

# --- Radial boundaries (cm) ---
r_inner_core  = openmc.ZCylinder(r=60.0)    # Inner fuel zone outer radius
r_outer_core  = openmc.ZCylinder(r=90.0)    # Outer fuel zone outer radius
r_reflector   = openmc.ZCylinder(r=110.0)   # Reflector outer radius
r_shield      = openmc.ZCylinder(r=130.0,   boundary_type='vacuum')

# --- Axial boundaries (cm) ---
z_bot_shield   = openmc.ZPlane(z0=-110.0,  boundary_type='vacuum')
z_bot_refl     = openmc.ZPlane(z0=-90.0)
z_bot_fuel     = openmc.ZPlane(z0=-60.0)
z_top_fuel     = openmc.ZPlane(z0= 60.0)
z_top_refl     = openmc.ZPlane(z0= 90.0)
z_top_shield   = openmc.ZPlane(z0= 110.0,  boundary_type='vacuum')

# --- Inner core fuel region ---
inner_region = (
    -r_inner_core &
    +z_bot_fuel   & -z_top_fuel
)
inner_cell = openmc.Cell(name='inner_fuel', fill=inner_fuel, region=inner_region)

# --- Outer core fuel region ---
outer_region = (
    +r_inner_core & -r_outer_core &
    +z_bot_fuel   & -z_top_fuel
)
outer_cell = openmc.Cell(name='outer_fuel', fill=outer_fuel, region=outer_region)

# --- Radial reflector (axial fuel span) ---
radial_refl_region = (
    +r_outer_core & -r_reflector &
    +z_bot_fuel   & -z_top_fuel
)
radial_refl_cell = openmc.Cell(name='radial_reflector', fill=reflector,
                                region=radial_refl_region)

# --- Radial shield (axial fuel span) ---
radial_shield_region = (
    +r_reflector & -r_shield &
    +z_bot_fuel  & -z_top_fuel
)
radial_shield_cell = openmc.Cell(name='radial_shield', fill=shield,
                                  region=radial_shield_region)

# --- Axial reflectors (top and bottom, within core radius) ---
bot_refl_region = -r_outer_core & +z_bot_refl & -z_bot_fuel
bot_refl_cell   = openmc.Cell(name='bot_reflector', fill=reflector, region=bot_refl_region)

top_refl_region = -r_outer_core & +z_top_fuel & -z_top_refl
top_refl_cell   = openmc.Cell(name='top_reflector', fill=reflector, region=top_refl_region)

# --- Axial shields (top and bottom) ---
bot_shield_region = -r_shield & +z_bot_shield & -z_bot_refl
bot_shield_cell   = openmc.Cell(name='bot_shield', fill=shield, region=bot_shield_region)

top_shield_region = -r_shield & +z_top_refl & -z_top_shield
top_shield_cell   = openmc.Cell(name='top_shield', fill=shield, region=top_shield_region)

# --- Sodium fill for all remaining space inside the vacuum boundary ---
# (coolant fills gaps not occupied by solid regions above)
sodium_region = (
    -r_shield &
    +z_bot_shield & -z_top_shield &
    ~inner_region &
    ~outer_region &
    ~radial_refl_region &
    ~radial_shield_region &
    ~bot_refl_region &
    ~top_refl_region &
    ~bot_shield_region &
    ~top_shield_region
)
sodium_cell = openmc.Cell(name='sodium', fill=sodium, region=sodium_region)

# --- Universe and geometry ---
main_u = openmc.Universe(cells=[
    inner_cell, outer_cell,
    radial_refl_cell, radial_shield_cell,
    bot_refl_cell, top_refl_cell,
    bot_shield_cell, top_shield_cell,
    sodium_cell
])

# FIX: original passed undefined 'main_u'; now it is defined above
geom = openmc.Geometry(main_u)
geom.export_to_xml()

# =============================================================================
# SETTINGS
# =============================================================================
lower_left  = [-130, -130, -110]
upper_right = [ 130,  130,  110]
uniform_dist = openmc.stats.Box(lower_left, upper_right)
src = openmc.IndependentSource(
    space=uniform_dist,
    constraints={"fissionable": True}
)

settings = openmc.Settings()
settings.source   = src
settings.batches  = 100
settings.inactive = 10
settings.particles = 1000
settings.export_to_xml()

# =============================================================================
# EXECUTION
# =============================================================================
openmc.run()
