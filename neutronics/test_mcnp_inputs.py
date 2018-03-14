from mcnp_inputs import HomogeneousInput

def test_homog_comp():
    """Test fuel, clad, coolant homogenization process.
    """
    exp_comp = {
        6000	: -0.0009387626,
        7015	: -0.0247170148,
        8016	: -0.0024120083,
        13027	: -0.0002288825,
        14000	: -0.0001455693,
        15031	: -6.40871086139384E-06,
        16000	: -6.40871086139384E-06,
        22000	: -0.0004119886,
        24000	: -0.0086975362,
        25055	: -0.0001455693,
        26000	: -0.007782006,
        27059	: -0.0004165662,
        28000	: -0.0240326657,
        29000	: -0.0001249699,
        41093	: -0.0023460459,
        42000	: -0.0013961834,
        74180	: -0.0005946713,
        74182	: -0.1327808871,
        74183	: -0.0720987003,
        74184	: -0.1552139565,
        74186	: -0.1455896584,
        92235	: -0.3779222222,
        92238	: -0.041991358
               }

    # build example input file
    obs = HomogeneousInput(15, 0.6, 15)
    obs.homog_core()

    obs_str = obs.fuel_string.split('\n')

    # extract compositions from string
    for line in obs_str:
        if 'm' not in line:
            data = list(filter(None, line.split(' ')))
            print(data)
            assert abs(exp_comp[float(data[0])] - float(data[1])) < 1e-6
