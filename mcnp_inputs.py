import math
import tarfile
import os
import sys
import argparse
from base_input import base_string

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
        ########################################
        
        # fuel mass fractions
        U_fraction = {'U' : 0.94441, 'N' : 0.05559}
        Uran_comp =  {'235' : enrich, '238' : 1 - enrich}
        
        fuel_comp = {7015 : U_fraction['N'],
                     92235 : U_fraction['U'] * Uran_comp['235'],
                     92238 : U_fraction['U'] * Uran_comp['238']
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
        


        fuel_mat = "\
m1
    7015 -{0}
    92235 -{1}
    92238 -{2}".format(fuel_comp['N'], 
                       fuel_comp['U'] * Uran_comp['235'],
                       fuel_comp['U'] * Uran_comp['238'])


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
    
    tarball = tarfile.open('radius_sweep.tar.gz', 'w|gz')
    input_list = open('input_list.txt', 'w')
     
    for radius in range(5, 45, 1):
        inp = HomogeneousInput(radius, 15)
        print(radius)
        infilename = inp.write_input()
        input_list.write(infilename+'\n')
        tarball.add(infilename)

    input_list.close()
    tarball.add('input_list.txt')
