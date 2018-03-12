import math
import tarfile
import os
import sys
import argparse
from base_input import base_string

def merge_comps(compA, compB):
    """
    """
    for isotope in compB:
        if isotope in compA.keys():
            compA[isotope] += compB[isotope]
        else:
            compA.update({isotope : compB[isotope]})

    return compA

class HomogeneousInput:
    """Class to write MCNP input files for pin cell modeling in an infinite
    lattice (in the x-y directions).
    """
    
    # base template string modified by the methods below
    base_string = base_string

    def __init__(self, radius, length):
        """Initialize parameters.
        """
        self.L = length
        self.r = radius
        self.vol = self.r**2 * math.pi * self.L
        self.refl_t = 0.5 * self.r
        self.z_ref_t = self.L * 0.5
        self.vol_refl = ((self.r + self.refl_t)**2 - self.r**2)*math.pi * self.L
    
    def homogenize_fuel_clad_cool(self, frac_fuel, enrich=0.9):
        """
        """
        ############ Channel Dimensions ########
        # Radius = 0.5 cm, clad_thick = 0.0031 #
        clad_to_cool = 0.127844                #
        frac_cool = 1 - frac_fuel              #
        frac_clad = frac_cool * clad_to_cool   #
        uran_frac = 0.6 # UN in tungsten       #
        ########################################
        # densities
        rho_W = 19.3
        rho_cool = 87.8
        rho_UN = 11.3 # g/cc
        rho_In = 8.19
        
        mass_fuel = frac_fuel * uran_frac * rho_UN
        mass_matr = frac_fuel * (1 - uran_frac) * rho_W
        mass_cool = frac_cool * rho_cool
        mass_clad = frac_clad * rho_In
        mass_total = mass_fuel + mass_matr + mass_cool + mass_clad

        # fuel mass fractions
        U_fraction = {'U' : 0.94441, 'N' : 0.05559}
        Uran_comp =  {'235' : enrich, '238' : 1 - enrich}
        
        fuel_comp = {7015 : U_fraction['N'],
                     92235 : U_fraction['U'] * Uran_comp['235'],
                     92238 : U_fraction['U'] * Uran_comp['238']
                    }

        matr_comp = {74180 : 1.1746e-03,
                     74182 : 2.6227e-01,
                     74183 : 1.4241e-01,
                     74184 : 3.0658e-01,
                     74186 : 2.8757e-01
                    }

        clad_comp = {6000   : 0.00073,
                     13027  : 0.005,
                     14000  : 0.00318,
                     15031  : 0.00014,
                     16000  : 0.00014,
                     22000  : 0.009,
                     24000  : 0.19,
                     25055  : 0.00318,
                     26000  : 0.17,
                     28000  : 0.525,
                     27059  : 0.0091,
                     29000  : 0.00273,
                     41093  : 0.05125,
                     42000  : 0.0305
                    }

        cool_comp = {6000 : 0.272912,
                     8016 : 0.727088
                    }


        components = [fuel_comp, matr_comp, clad_comp, cool_comp]
        
        # get isotopic masses for fuel, matrix, coolant, cladding
        print(fuel_comp)
        fuel_comp = {isotope : fuel_comp[isotope] 
                     * mass_fuel for isotope in fuel_comp}
        print(fuel_comp)
        matr_comp = {isotope : matr_comp[isotope] 
                     * mass_matr for isotope in matr_comp}
        clad_comp = {isotope : clad_comp[isotope]
                     * mass_clad for isotope in clad_comp}
        cool_comp = {isotope : cool_comp[isotope]
                     * mass_cool for isotope in cool_comp}

        homog_comp = {}
        for frac in components:
            homog_comp = merge_comps(homog_comp, frac)
        
        # write the mcnp string
        self.fuel_string = "\
m1\n"
        for isotope in sorted(homog_comp.keys()):
            self.fuel_string += "    {0} -{1}\n".format(isotope,
                    round(homog_comp[isotope] / mass_total, 7))

        print(self.fuel_string)

        

    def write_input(self):
        """ Write MCNP6 input files.
        This function writes the MCNP6 input files for the leakage experiment using
        the template input string. It writes a bare and reflected core input file
        for each core radius.
        """
        templ = self.base_string
        file_string = templ.substitute(r_core = self.r,
                                       core_z = self.L,
                                       r_refl = self.r + self.refl_t,
                                       refl_min = -self.z_ref_t,
                                       refl_max = self.L + self.z_ref_t,
                                       fuel_vol = self.vol,
                                       refl_vol = self.vol_refl)
        # write the file
        filename = 'r_{0}_{1}.i'.format(self.r, self.L)
        ifile = open(filename, 'w')
        ifile.write(file_string)
        ifile.close()

        return filename

if __name__ == '__main__':
   test = HomogeneousInput(15,15)
   test.homogenize_fuel_clad_cool(0.6)
