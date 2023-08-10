from math import pi as PI
import numpy as np

# module
from src.format_data import formatData

# NOTE: Maintain variable names when modifying parameters/functions

# =========== wing twist functions ===========
                                                        # units     # description:
                                                        #           local wing chord divided by wing span - defines shape of the wing. 
def chord_ratio(x):					
	return (1 - x**2)**0.5 * (4/(PI*12))
                                                        # deg       glide twist angle at a given dimensionless spanwise position. x=1 is tip
def twist_glide(x):
	return 5
		
def twist(x):
	return np.arctan( 0.6 * np.abs(x) ) * 180 / PI
                                                        # deg       upstroke twist angle at a given dimensionless spanwise position. 
def twist_up(x):
	return twist_glide(x) + twist(x)
                                                        # deg       downstroke twist angle at a given dimensionless spanwise position.	
def twist_down(x):
	return twist_glide(x) - twist(x)

# =========== design parameters ==============
                                # units     # description
# integration steps
n_step      = 201               #           number of wingspan subdivisions to calculate lift. Larger number slows down calculations
# flight parameters
amplitude   = 55                # deg       flapping amplitude
dihedral    = 2.5               # deg       average dihedral angle
ld_ratio    = 4                 # %         lift to drag ratio of entire aircraft
gravity     = 9.81              # m/^2      acceleration due to gravity 
mass_total  = 300e-3            # kg        total mass of the aircraft
mass_wing   = 52e-3             # kg        mass of both half wingspans
air_density = 1.204             # kg/m^3    density of atmosphere
area        = 0.110             # m^2       area of entire wing
# motor properties
motor_kv    = 3100              # rpm/volt  rpm for every volt fed to the motor
voltage     = 3.7*3             # volt      voltage source of motor
i_stall     = 6                 # amp       stall current of motor
i_noload    = 0.5               # amp       no-load current of motor
# spring geometry
xoffset     = -10e-3            # m         horizontal offset of wing shoulder
yoffset     = 80e-3             # m         vertical offset of wing shoulder
rspring     = 50e-3             # m         radius of wing shoulder

# ========== solve design ========== [do not edit]
formatData( [n_step,   chord_ratio, twist_glide, twist_up, twist_down                               \
            ,ld_ratio, amplitude,   dihedral,    gravity,  mass_total, mass_wing, air_density, area \
            ,voltage,  i_stall,     i_noload,    motor_kv                                          \
            ,xoffset,  yoffset,     rspring] )
