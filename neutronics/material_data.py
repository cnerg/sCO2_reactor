"""Mass fraction compositions for reactor materials.
"""
# fuel fraction in cermet 
vfrac_UN = 0.6 # (Webb and Charit 2012)
# densities [g/cc]
rho_W = 19.3   # PNNL mat compend #331 
rho_UN = 14.31 # PNNL mat compend #337
rho_In = 8.19  # PNNL mat compend #156

# Uranium Nitride
UN = {92000 : 0.94441, 7015 : 0.05559}

# Tungsten
W = {74180 : 1.1746e-03,
     74182 : 2.6227e-01,
     74183 : 1.4241e-01,
     74184 : 3.0658e-01,
     74186 : 2.8757e-01}

# Inconel-718
In = {6000   : 0.00073,
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
      42000  : 0.0305}

# Carbon Dioxide
CO2 = {6000 : 0.272912,
       8016 : 0.727088}

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
