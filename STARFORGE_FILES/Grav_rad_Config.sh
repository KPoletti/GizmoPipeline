HYDRO_MESHLESS_FINITE_MASS
BOX_SPATIAL_DIMENSION=3
BOX_PERIODIC
GRAVITY_NOT_PERIODIC
EOS_GAMMA=(5.0/3.0)
SINGLE_STAR_STARFORGE_DEFAULTS
COOLING # enables detailed cooling physics, see the FIRE-3 paper for details
METALS
# Radiation
RT_M1 # Simulation method

# Radiation Bands
RT_LYMAN_WERNER                        # specific lyman-werner [narrow H2 dissociating] band 11.2 to 13.6 eV, 
RT_PHOTOELECTRIC                       # far-uv (8-13.6eV): track photo-electric heating photons + their dust interactions
RT_NUV                                 # near-UV: 1550-3600 Angstrom (where direct stellar emission dominates) 8.00 to 3.44 eV
RT_OPTICAL_NIR                         # optical+near-ir: 3600 Angstrom-3 micron (where direct stellar emission dominates) 3.44eV to 0.41 eV
RT_ISRF_BACKGROUND

# I am unsure if this does anything
SINGLE_STAR_FB_RAD # enables explicit radiative transfer, include emission/absorption of gas and dust, and stellar emission
ADAPTIVE_TREEFORCE_UPDATE=0.0625 # Enables gravity optimization for gas cells introduced in
# RT_SPEEDOFLIGHT_REDUCTION=0.01

# Outputs
OUTPUT_RT_RAD_FLUX
RT_RAD_PRESSURE_OUTPUT
OUTPUT_TIMESTEP
OUTPUT_IN_DOUBLEPRECISION      # snapshot files will be written in double precision
OUTPUT_TEMPERATURE
STOP_WHEN_BELOW_MINTIMESTEP

