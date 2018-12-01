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

def write_mat_string(mat_num, comp, XS, sab=None):
    """Write the fuel composition in MCNP-format.

    Arguments:
    ----------
        homog_comp (dict): normalized isotopic mass fractions
    
    Modified Attributes:
    --------------------
        fuel_string (str): MCNP-style material card
    """
    # Initialize material card string
    mat_string = "m%s\n" % mat_num 

    # loop through isotopes and write mcnp-style mass fractions
    for idx, isotope in enumerate(sorted(comp.keys())):
        if abs(comp[isotope]) < 1e-10:
            continue
        mat_string += "     {0}.{1} -{2:.5e}".format(isotope, XS, 
                                                     comp[isotope])
        # no endline character for last isotope
        if idx < len(comp.keys()) - 1:
            mat_string += '\n'

    if sab:
        mat_string += "\nmt%s " % mat_num
        mat_string += sab
            
    return mat_string

class HomogeneousInput:
    """Write Homogeneous Input File.
    Class to write homogeneous MCNP burnup input files.
    """

    def __init__(self, **config):
        """Initialize geometric reactor parameters.

        Initialized Attributes:
        -----------------------
            z (float): reactor core height
            r (float): reactor core radius
            frac_fuel (float): fuel to coolant channel fraction
            refl_t (float): reflector thickness
        """
        
        config = config['config']
        
        # geometric quantities
        self.r = config.get('core_r')
        self.z = self.r*config.get('AR', 1)
        self.r_cool = config.get('r_cool', 0.5)
        self.c = config.get('c', 0.031)
        self.vfrac_fuel = config.get('fuel_frac')
        self.t_refl = self.r*config.get('ref_mult')

        # material choices
        self.fuel =     config.get('fuel', 'UN')
        self.enrich =   config.get('enrich', 0.93)
        self.matr =     config.get('matr', False)
        self.cool =     config.get('cool', 'CO2')
        self.clad =     config.get('clad', 'Inconel-718')
        self.refl =     config.get('refl', 'C')
        self.rho_cool = config.get('rho_cool', 0.079082)
        
        self.Q_therm = config.get('power')

        # calculate core volume
        self.core_vol = self.r * self.r * self.z * math.pi
        
        # get string for reflector material
        self.refl_mat_string = write_mat_string(2, md.mats[self.refl]['comp'],
                                                   md.mats[self.refl]['XS'],
                                                   md.mats[self.refl]['sab'])
        
        self.refl_mass = ((self.r + self.t_refl)**3 * math.pi - self.core_vol)\
                        * md.mats[self.refl]['rho']

    def homog_core(self):
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
        fuel_in_matr_frac = 1
        mfrac_matr = 0
        
        c_in_cool_frac = 1 - self.r_cool**2 / (self.r_cool + self.c)**2
        
        if self.matr:
            fuel_in_matr_frac = 0.6
            mfrac_matr = self.vfrac_fuel * md.mats[self.matr]['rho'] *\
                                          (1-fuel_in_matr_frac)
        # volume-weighted densities/masses
        mfrac_fuel = self.vfrac_fuel * md.mats[self.fuel]['rho']*\
                                       fuel_in_matr_frac
        mfrac_cool = (1- self.vfrac_fuel) * self.rho_cool * (1 - c_in_cool_frac)
        mfrac_clad = (1- self.vfrac_fuel) * c_in_cool_frac *\
                                            md.mats[self.clad]['rho']

        # smeared density
        self.rho = round(mfrac_fuel + mfrac_cool + mfrac_matr + mfrac_clad, 6)
        self.core_mass = self.core_vol * self.rho
        
        # get UN composition
        fuel_comp, self.MM_fuel = md.enrich_fuel(self.enrich, self.fuel)
        components = {
        # get isotopic mass fractions for fuel, matrix, coolant, cladding
            'normed_fuel_mfrac' : {iso : comp*mfrac_fuel / self.rho 
                                for iso, comp in fuel_comp.items()},
            'normed_cool_mfrac' : {iso : comp*mfrac_cool / self.rho 
                                for iso, comp in md.mats[self.cool]['comp'].items()},
            'normed_clad_mfrac' : {iso : comp*mfrac_clad / self.rho
                                for iso, comp in md.mats[self.clad]['comp'].items()}
                     }

        if self.matr:
            components.update({
            'normed_matr_mfrac' : {iso : comp*mfrac_matr / self.rho 
                                 for iso, comp in md.mats[self.matr]['comp'].items()}})
        
        # homogenize the material by merging components
        homog_comp = {}
        for frac in components:
            homog_comp = merge_comps(homog_comp, components[frac])
        
        self.core_mat_string = write_mat_string(1, homog_comp, 
                                                md.mats[self.fuel]['XS'])
        
    def write_input(self, name, header="Fuel Fraction Experiment"):
        """ Write MCNP6 input files.
        This function writes the MCNP6 input files for the leakage experiment using
        the template input string. It writes a bare and reflected core input file
        for each core radius.
        """
        # load template, substitute parameters and write input file
        input_tmpl = open('base_input.txt')
        templ = Template(input_tmpl.read())
        file_string = templ.substitute(model_information = header,
                                       r_core = self.r,
                                       core_z = self.z,
                                       refl_t = self.t_refl,
                                       refl_z = self.z + 2*self.t_refl,
                                       r_refl = self.r + self.t_refl,
                                       fuel_string = self.core_mat_string,
                                       refl_mat = self.refl_mat_string,
                                       volume = self.core_vol,
                                       fuel_rho = self.rho,
                                       shift_core = -(self.z + self.t_refl + 100),
                                       refl_rho = md.mats[self.refl]['rho'],
                                       ksrc_z = -(self.z + 99),
                                       thermal_power = self.Q_therm)
        # write the file
        ifile = open(name, 'w')
        ifile.write(file_string)
        ifile.close()


if __name__=='__main__':
    
    config = {'fuel' : 'UN',
              'matr' : 'W',
              'cool' : 'CO2',
              'clad' : 'Inconel-718',
              'rho_cool' : 233.89e-3,
              'fuel_frac' : 0.33,
              'core_r' : 16,
              'ref_mult' : 0.2
             }

    input = HomogeneousInput(config=config) 
    input.homog_core()
    input.write_input('test')
