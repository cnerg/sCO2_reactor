import math
from pyne.material import Material, MaterialLibrary
from string import Template
import material_data as md

def build_pyne_matlib(nucdata_file=None):
    """Fetch pyne material from compendium.

    This function builds PyNE material library
    for UWNR non-fuel components defined in reactor_data.py

    Arguments:
        nucdata_file (str)[-]: Filename 'nuc_data.h5' with its full path.
    """
    h5path = "/home/alex/.local/lib/python3.5/site-packages/pyne/nuc_data.h5"
    
    if nucdata_file:
        h5path = nucdata_file
    
    # Initialize material libraries.
    raw_matlib = MaterialLibrary()

    # Write entire PyNE material library.
    raw_matlib.from_hdf5(h5path,
                         datapath="/material_library/materials",
                         nucpath="/material_library/nucid")
    
    return raw_matlib

class HomogeneousInput:
    """Write Homogeneous Input File.
    Class to write homogeneous MCNP burnup input files.
    """
    
    def __init__(self, radius, length, pnnl_mats, thick_refl=15):
        """Initialize geometric reactor parameters.

        Initialized Attributes:
        -----------------------
            z (float): reactor core height
            r (float): reactor core radius
            frac_fuel (float): fuel to coolant channel fraction
            matlib (MaterialLibrary): PyNE material library
        """
        self.z = length
        self.r = radius
        self.refl_t = thick_refl
        self.matlib = pnnl_mats

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

        self.cell_vol = v_cool + v_clad + v_cermet
        # calculate normalized vfracs from total cell volume
        self.vfrac_cool = v_cool / self.cell_vol
        self.vfrac_clad = v_clad / self.cell_vol
        self.vfrac_cermet = v_cermet / self.cell_vol

        
    def homog_core(self, enrich=0.9, r_cool=0.5, 
                         PD=1.48, rho_cool=0.087, c=0.031):
        """Homogenize the fuel, clad, and coolant.
        
        Arguments:
        ----------
            enrich (float) (opt): uranium enrichment
            r_cool (float) (opt): coolant channel radius
            PD (float) (opt): pitch to diameter ratio
            rho_cool (float) (opt): coolant density
            c (float) (opt): claddinng thickness
        
        Modified Attributes:
        --------------------
            rho (float): fuel density
            fuel_string (str): MCNP material string 
        """
        # get volumes, volume fractions
        self.calc_vol_vfrac(r_cool, PD, c)
        
        # volume-weighted densities/masses
        fracs = {'fuel' : {
                      'volfrac' : self.vfrac_cermet*md.vfrac_UN,
                      'rho' : md.rho_UN},
                 'matr' : {
                      'volfrac' : self.vfrac_cermet*(1 - md.vfrac_UN),
                      'rho' : md.rho_W},
                 'cool' : {
                      'volfrac' : self.vfrac_cool,
                      'rho' : rho_cool}, 
                 'clad' : {
                      'volfrac' : self.vfrac_clad,
                      'rho' : md.rho_In}
                 }

        core_mass = 0
        # get component masses
        for comp in fracs:
            comp_mass = fracs[comp]['volfrac']*fracs[comp]['rho']
            fracs[comp].update({'mass' : comp_mass})
            core_mass += comp_mass
            
        self.homog_mat = Material()
        # mix normalized mass fractions
        for mat in fracs:
            # get UN material from custom composition
            if mat == 'fuel':
                pyne_mat = Material(md.enrich_fuel(enrich))
            else:
                pyne_mat = self.matlib[md.mats[mat]]
            pyne_mat.mass = fracs[mat]['mass'] / core_mass
            self.homog_mat += pyne_mat
        
        # total density [g/cc]
        self.rho = core_mass / self.core_vol

    def write_mat_string(self):
        """Prep the homogenized fuel material and write it in MCNP material card
        format.

        Modified Attributes:
        --------------------
            fuel_string (str): MCNP-style material card
        """
        # delete bad apple isotopes
        missing_nuclides = ['8018', '8017']
        del self.homog_mat[missing_nuclides]
        # set material number
        self.homog_mat.metadata['mat_number'] = 1
        # write mcnp-form string
        self.fuel_string = self.homog_mat.mcnp().strip('\n')

    def write_input(self):
        """ Write MCNP6 input files.
        This function writes the MCNP6 input files for the leakage experiment using
        the template input string. It writes a bare and reflected core input file
        for each core radius.

        Returns:
        --------
            filename (str): name of written MCNP6 input file
        """
        # load template, substitute parameters and write input file
        self.homog_core()
        self.write_mat_string()
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

if __name__=='__main__':
    matlib = build_pyne_matlib()
    test = HomogeneousInput(10, 5, matlib)
    test.write_input()
