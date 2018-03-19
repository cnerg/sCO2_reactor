from mcnp_inputs import HomogeneousInput

def test_homog_comp():
    """Test fuel, clad, coolant homogenization process.
    """
    exp_comp = {
            6000 : 0.0011076084,
            7015 : 0.0277709767,
            8016 : 0.0028588435,
            13027 : 0.0002365987,
            14000 : 0.0001504768,
            15031 : 6.62476374069704E-06,
            16000 : 6.62476374069704E-06,
            22000 : 0.0004258777,
            24000 : 0.0089907508,
            25055 : 0.0001504768,
            26000 : 0.008044356,
            28000 : 0.024842864,
            27059 : 0.0004306096,
            29000 : 0.0001291829,
            41093 : 0.0024251367,
            42000 : 0.0014432521,
            74180 : 0.0005276074,
            74182 : 0.1178065703,
            74183 : 0.0639677953,
            74184 : 0.1377097583,
            74186 : 0.129170837,
            92235 : 0.4246171846,
            92238 : 0.0471796872
               }


    # build example input file
    obs = HomogeneousInput(15, 0.6, 15)
    comp = obs.homog_core()
    obs.write_mat_string(comp)
    obs_str = obs.fuel_string.split('\n')

    # extract compositions from string
    for line in obs_str:
        if 'm' not in line:
            data = list(filter(None, line.split(' ')))
            print(data)
            assert abs(exp_comp[float(data[0])] - -float(data[1])) < 1e-6
