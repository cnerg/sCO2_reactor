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

def enrich_fuel(enrich, fuel=UN):
    """Enrich Uranium fuel and mix with specified compound.

    Arguments:
    ----------
        enrich (float): U-235 mass conc.
        fuel (dict) (optional): mass composition of fuel
    Returns:
    --------
        fuel_comp (dict): isotopic mass vector of Uranium fuel compound.
    """
    fuel_comp = {}
    # build fuel comp. starting with bounded compound
    for isotope in fuel:
        if isotope == 92000:
            # add enriched Uranium to fuel composition
            fuel_comp.update({92235 : fuel[92000]*enrich,
                              92238 : fuel[92000]*(1-enrich) })
        else:
            fuel_comp.update({isotope : fuel[isotope]})
    
    return fuel_comp 
