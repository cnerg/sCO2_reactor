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

def enrich_fuel(enrich, fuel):
    """Enrich Uranium fuel and mix with specified compound.

    Arguments:
    ----------
        enrich (float): U-235 mass conc.
        fuel (dict) (optional): mass composition of fuel
    Returns:
    --------
        fuel_comp (dict): isotopic mass vector of Uranium fuel compound.
    """
    # build fuel comp. starting with bounded compound
    del fuel[92236, 92234]
    
    comp = dict(fuel)
    
    mfrac_U = comp[922350000] + comp[922380000]

    fuel[922350000] = mfrac_U * enrich
    fuel[922380000] = mfrac_U * (1-enrich)
    
    return fuel 
