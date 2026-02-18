%matplotlib inline
import openmc

## Materials
# Uraniums
U235 = openmc.Material(name='U235')
U235.add_nuclide('U235', 1.0)
U235.set_density('g/cm3', 19.1)

U238 = openmc.Material(name='U238')
U238.add_nuclide('U238', 1.0)
U238.set_density('g/cm3', 19.1)

# Transuranics
Np237 = openmc.Material(name='Np237')
Np237.add_nuclide('Np237', 1.0)
Np237.set_density('g/cm3', 20.5)

Pu238 = openmc.Material(name='Pu238')
Pu238.add_nuclide('Pu238', 1.0)
Pu238.set_density('g/cm3', 19.8)

Pu239 = openmc.Material(name='Pu239')
Pu239.add_nuclide('Pu239', 1.0)
Pu239.set_density('g/cm3',19.8 )

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

Am242m = openmc.Material(name='Am242')
Am242m.add_nuclide('Am242', 1.0)
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

# Other
Zr = openmc.Material(name='Zr')
Zr.add_element('Zr', 1.0)
Zr.set_density('g/cm3', 13.5)

## Material mixtures
# Materials
U = openmc.Material.mix_materials(
  [],
  [],
  'wo'
)

TRU = openmc.Material.mix_materials(
  [Np237, Pu238, Pu239, Pu240, Pu241, Pu242, Am241, Am242m, Am243, Cm243, Cm244, Cm245],
  [0.0472, 0.0218, 0.4734, 0.2282, 0.0842, 0.0684, 0.0561, 0.0001, 0.0156, 0.0000, 0.0046, 0.0004],
  'wo'
) # Recycled Core

WGPu = openmc.Material.mix_materials(
  [Pu238, Pu239, Pu240, Pu241, Pu242]
  [0.01, 93.81, 5.81, 0.35, 0.02]
  'wo'
) # Startup Core

# Components
inner = openmc.Material.mix_materials(
  [],
  [],
  'wo')

outer = openmc.Material.mix_materials(
  [],
  [],
  'wo')

cladding = openmc.Material.mix_materials(
  [],
  [],
  'wo')

reflector = openmc.Material.mix_materials(
  [],
  [],
  'wo')

shield = openmc.Material.mix_materials(
  [],
  [],
  'wo')

primary_control = openmc.Material.mix_materials(
  [],
  [],
  'wo')

secondary_control = openmc.Material.mix_materials(
  [],
  [],
  'wo')

## Materials file
# Instantiate a Materials collection and export to xml
materials_file = openmc.Materials([inner, outer, sodium, clad])
materials_file.export_to_xml()

## Assembly Geometry

## Inner assemblies

## Outer assemblies

## Primary Control

## Secondary Control

## Reflector

## Shield

## Core Geometry

## Geometry file
geom = openmc.Geometry(main_u)
geom.export_to_xml()

## Simulation parameters
lower_left = [-300, -300, -50]
upper_right = [300, 300, 50]
uniform_dist = openmc.stats.Box(lower_left, upper_right)
src = openmc.IndependentSource(space=uniform_dist, constraints={"fissionable": True})

settings = openmc.Settings()
settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 1000

settings.export_to_xml()

## EXECUTION
openmc.run()
