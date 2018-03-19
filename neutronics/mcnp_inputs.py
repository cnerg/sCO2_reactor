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
    
    def __init__(self, radius, length, thick_refl=15):
        """Initialize geometric reactor parameters.

        Initialized Attributes:
        -----------------------
            z (float): reactor core height
            r (float): reactor core radius
            frac_fuel (float): fuel to coolant channel fraction
        """
        self.z = length
        self.r = radius
        self.refl_t = thick_refl

    def calc_vol_vfrac(self, r_cool, PD, c):
        """Get volumes for core and reflector regions. Calculate the volume
        fraction of core components. The volume fractions are used to homogenize
        the core material.
        
        Modified Attributes:
        --------------------
            core_vol (float): homogenized core volume
            refl_vol (float): reflector volume
            vfrac_cool (float): coolant volume fraction
            vfrac_cermet (float): cermet matrix volume fraction
            vfrac_clad (float): cladding volume fraction
        """
        # core and reflector volume required for depletion calculation
        self.core_vol = self.r**2 * math.pi * self.z
        self.refl_vol = ((self.r + self.refl_t)**2 - self.r**2)*math.pi * self.z
        
        pitch = 2*r_cool*PD
        # calculate 'volumes' for fixed length
        v_cool = (r_cool ** 2 * math.pi)
        # clad volume fraction
        v_clad = ((r_cool + c)**2 - r_cool**2)*math.pi
        # fuel volume fraction
        v_cermet = (math.sqrt(3)*pitch**2 / 2.0) - (r_cool + c) ** 2 * math.pi 

        cell_vol = v_cool + v_clad + v_cermet
        # calculate vfracs from total cell volume
        self.vfrac_cool = v_cool / cell_vol
        self.vfrac_clad = v_clad / cell_vol
        self.vfrac_cermet = v_cermet / cell_vol
        
        print(self.vfrac_cool)
        print(self.vfrac_clad)
        print(self.vfrac_cermet)
        
    def homog_core(self, enrich=0.9, r_cool=0.5, PD=1.48, rho_cool=0.087, c=0.031):
        """Homogenize the fuel, clad, and coolant.
        
        Arguments:
        ----------
            enrich (float) (optional): uranium enrichment
        
        Modified Attributes:
        --------------------
            rho (float): fuel density
        
        Returns:
        --------
            homog_comp (dict): homogenized, isotopic, core composition (mass %)
        """
        # get volumes, volume fractions
        self.calc_vol_vfrac(r_cool, PD, c)
        # volume-weighted densities/masses
        mfrac_fuel = self.vfrac_cermet * md.vfrac_UN * md.rho_UN
        mfrac_matr = self.vfrac_cermet * (1 - md.vfrac_UN) * md.rho_W
        mfrac_cool = self.vfrac_cool * rho_cool
        mfrac_clad = self.vfrac_clad * md.rho_In
        # smeared density
        self.rho = round(mfrac_fuel + mfrac_matr + mfrac_cool + mfrac_clad, 5)
        
        # get UN composition
        fuel_comp = md.enrich_fuel(enrich)

        components = {
        # get isotopic mass fractions for fuel, matrix, coolant, cladding
            'normed_fuel_mfrac' : {iso : comp*mfrac_fuel / self.rho 
                                 for iso, comp in fuel_comp.items()},
            'normed_matr_mfrac' : {iso : comp*mfrac_matr / self.rho 
                                 for iso, comp in md.W.items()},
            'normed_clad_mfrac' : {iso : comp*mfrac_clad / self.rho 
                                 for iso, comp in md.In.items()},
            'normed_cool_mfrac' : {iso : comp*mfrac_cool / self.rho 
                                 for iso, comp in md.CO2.items()}
                     }
        
        
        # homogenize the material by merging components
        homog_comp = {}
        for frac in components:
            homog_comp = merge_comps(homog_comp, components[frac])
        
        return homog_comp

    def write_mat_string(self, homog_comp):
        """Write the fuel composition in MCNP-format.

        Arguments:
        ----------
            homog_comp (dict): normalized isotopic mass fractions
        
        Modified Attributes:
        --------------------
            fuel_string (str): MCNP-style material card
        """
        # Initialize material card string
        self.fuel_string = "m1\n"

        # loop through isotopes and write mcnp-style mass fractions
        for idx, isotope in enumerate(sorted(homog_comp.keys())):
            self.fuel_string += "     {0} -{1}".format(isotope,
                    round(homog_comp[isotope], 7))
            # no endline character for last isotope
            if idx < len(homog_comp.keys()) - 1:
                self.fuel_string += '\n'
                
        
    def write_input(self):
        """ Write MCNP6 input files.
        This function writes the MCNP6 input files for the leakage experiment using
        the template input string. It writes a bare and reflected core input file
        for each core radius.
        """
        # homogenize material and write MCNP card
        core_comp = self.homog_core()
        self.write_mat_string(core_comp)
        # load template, substitute parameters and write input file
        input_tmpl = open('base_input.txt')
        templ = Template(input_tmpl.read())
        file_string = templ.substitute(cool_frac = self.vfrac_cermet,
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
        filename = 'r_{0}_{1}.i'.format(round(self.vfrac_cermet, 3), 
                                                 round(self.r, 3))
        ifile = open(filename, 'w')
        ifile.write(file_string)
        ifile.close()

        return filename
