import math
from string import Template
import material_data as md

def merge_comps(compA, compB):
    """Merge two compositions for homogenization. Combine like isotopes when
    necessary. This function takes two input compositions, and returns the union
    of the two compositions.

    Arguments:
    ----------
        compA (dict): first composition
        compB (dict): second composition

    Returns:
    --------
        compA (dict): merged comp. of A and B
    """
    for isotope in compB:
        if isotope in compA:
            compA[isotope] += compB[isotope]
        else:
            compA.update({isotope : compB[isotope]})

    return compA


class HomogeneousInput:
    """Write Homogeneous Input File.
    Class to write homogeneous MCNP burnup input files.
    """
    
    def __init__(self, radius, volfrac_fuel, length):
        """Initialize geometric reactor parameters.

        Initialized Attributes:
        -----------------------
            z (float): reactor core height
            r (float): reactor core radius
            frac_fuel (float): fuel to coolant channel fraction
        """
        self.z = length
        self.r = radius
        self.frac_fuel = volfrac_fuel
        self.refl_t = 15 # 15 cm reflector

    def calculate_volume(self):
        """Get volumes for core and reflector regions.
        
        Modified Attributes:
        --------------------
            vol (float): homogenized core volume
            vol_refl (float): reflector volume
        """
        self.core_vol = self.r**2 * math.pi * self.z
        self.refl_vol = ((self.r + self.refl_t)**2 - self.r**2)*math.pi * self.z

    
    def homog_core(self, enrich=0.9):
        """Homogenize the fuel, clad, and coolant.
        
        Arguments:
        ----------
            enrich (float) (optional): uranium enrichment
        
        Modified Attributes:
        --------------------
            rho (float): fuel density
            fuel_string (str): mcnp-style fuel string
        """
        ############ Channel Dimensions ################
        # Radius = 0.5 cm, clad_thick = 0.0031         #
        clad_to_cool = 0.127844                        #
        frac_clad = (1 - self.frac_fuel) * clad_to_cool# 
        frac_cool = 1 - frac_clad - self.frac_fuel     #
        uran_frac = 0.6 # UN in tungsten               #
        ################################################

        # densities
        rho_W = 19.3
        rho_cool = 0.087
        rho_UN = 11.3 # g/cc
        rho_In = 8.19

        # volume-weighted densities
        mass_fuel = self.frac_fuel * uran_frac * rho_UN
        mass_matr = self.frac_fuel * (1 - uran_frac) * rho_W
        mass_cool = frac_cool * rho_cool
        mass_clad = frac_clad * rho_In
        # smeared density
        self.rho = round(mass_fuel + mass_matr + mass_cool + mass_clad, 5)
        
        # fuel mass fractions
        Uran_comp =  {'235' : enrich, '238' : 1 - enrich}
        
        # Uranium Nitride
        fuel_comp = {7015 : md.UN[7015],
                     92235 : md.UN[92000] * Uran_comp['235'],
                     92238 : md.UN[92000] * Uran_comp['238'] }

        
        # get isotopic mass fractions for fuel, matrix, coolant, cladding
        fuel_mfrac = {isotope : fuel_comp[isotope] 
                     * (mass_fuel / self.rho) for isotope in fuel_comp}
        matr_mfrac = {isotope : md.W[isotope] 
                     * (mass_matr / self.rho) for isotope in md.W}
        clad_mfrac = {isotope : md.In[isotope]
                     * (mass_clad / self.rho) for isotope in md.In}
        cool_mfrac = {isotope : md.CO2[isotope]
                     * (mass_cool / self.rho) for isotope in md.CO2}

        homog_comp = {}
        components = [fuel_mfrac, matr_mfrac, clad_mfrac, cool_mfrac]
        # homogenize the material by merging components
        for frac in components:
            homog_comp = merge_comps(homog_comp, frac)
        
        # write the mcnp string
        self.fuel_string = "m1\n"
        endline = '\n'
        # loop through isotopes and write mcnp-style mass fractions
        for idx, isotope in enumerate(sorted(homog_comp.keys())):
            # no endline character for last isotope
            if idx == len(homog_comp.keys()) - 1:
                endline = ''
            self.fuel_string += "     {0} -{1}{2}".format(isotope,
                    round(homog_comp[isotope], 7), endline)
        
    def write_input(self):
        """ Write MCNP6 input files.
        This function writes the MCNP6 input files for the leakage experiment using
        the template input string. It writes a bare and reflected core input file
        for each core radius.
        """
        self.calculate_volume()
        self.homog_core()
        input_tmpl = open('base_input.txt')
        templ = Template(input_tmpl.read())
        file_string = templ.substitute(cool_frac = self.frac_fuel,
                                       r_core = self.r,
                                       core_z = self.z,
                                       r_refl = self.r + self.refl_t,
                                       refl_min = -self.refl_t,
                                       refl_max = self.z + self.refl_t,
                                       fuel_string = self.fuel_string,
                                       fuel_rho = self.rho,
                                       fuel_vol = self.core_vol,
                                       refl_vol = self.core_vol)
        # write the file
        filename = 'r_{0}_{1}.i'.format(round(self.frac_fuel, 3), 
                                                 round(self.r, 3))
        ifile = open(filename, 'w')
        ifile.write(file_string)
        ifile.close()

        return filename

