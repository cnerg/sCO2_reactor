"""Mass fraction compositions for reactor materials.
"""
# fuel fraction in cermet 
vfrac_UN = 0.6 # (Webb and Charit 2012)
# densities [g/cc]
rho_W = 19.3   # PNNL mat compend #331 
rho_UN = 14.31 # PNNL mat compend #337
rho_In = 8.19  # PNNL mat compend #156

mats = {'clad' : 'Inconel-718',
        'cool' : 'Carbon Dioxide',
        'matr' : 'Tungsten' }

def enrich_fuel(enrich):
    """Enrich Uranium fuel and mix with specified compound.

    Arguments:
    ----------
        enrich (float): U-235 mass conc.
    Returns:
    --------
        fuel_comp (dict): isotopic mass vector of Uranium fuel compound.
    """
    UN = {92000 : 0.94441, 7015 : 0.05559}
    
    fuel = {92235 : UN[92000]*enrich,
            92238 : UN[92000]*(1-enrich),
            7015 : UN[7015]
           }

    return fuel 
