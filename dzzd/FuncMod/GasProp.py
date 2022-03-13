# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 07:24:29 2019

@author: dell
"""
import numpy
from scipy import interpolate
import math
import scipy

def GasProp(MF_in, FuelF):
    global frac
    mass_O2 = MF_in * 0.2314 - FuelF *(4 * 15.9994 / 16.0426 * c1 \
				+ 7 * 15.9994/30.0694 * c2 \
				+ 10 * 15.9994/44.0962 * c3 \
                +13 * 15.9994/58.1222 * c4 \
                + 13 * 15.9994/58.1222 * ic4 \
                +16 * 15.9994/72.15 * c5 \
                + 16 * 15.9994/72.15 * ic5 \
                +19 * 15.9994/86.17 * c6)
    mass_N2	= MF_in * 0.7552 + FuelF * N2
    mass_CO2 = MF_in * 0.0005 + FuelF *( co2 \
					+ 44.0098/16.0426 * c1 \
					+ 2 * 44.0098/30.0694 * c2 \
					+ 3* 44.0098/44.0962 * c3 \
					+4 * 44.0098/58.122 * c4 \
                    + 4 * 44.0098/58.122 * ic4 \
                    +5 * 44.0098/72.15 * c5 \
                    + 5 * 44.0098/72.15 * ic5 \
                    +6 * 44.0098/86.17 * c6)
    mass_H2O = FuelF *( 2 * 18.0152/16.0426 * c1 \
					+ 3 * 18.0152/30.0694 * c2 \
					+ 4 * 18.0152/44.0962 * c3 \
                    +5 * 18.0152/58.122 * c4 \
                    + 5 * 18.0152/58.1222 * ic4 \
                    +6 * 18.0152/72.15 * c5 \
                    + 6 * 18.0152/72.15 * ic5 \
                    +7 * 18.0152/86.17 * c6)
    mass_Ar = MF_in * 0.0129


    TotalMass = mass_O2 + mass_N2 + mass_CO2 + mass_H2O + mass_Ar
    mass_O2 = mass_O2 / TotalMass
    mass_N2 = mass_N2 / TotalMass
    mass_CO2 = mass_CO2 / TotalMass
    mass_H2O = mass_H2O / TotalMass
    mass_Ar = mass_Ar / TotalMass


    fluid = CoolProp.AbstractState('REFPROP','oxygen&nitrogen&co2&water&argon')
    fluid.set_mass_fractions([mass_O2, mass_N2, mass_CO2, mass_H2O, mass_Ar])
    frac = fluid.get_mole_fractions()

