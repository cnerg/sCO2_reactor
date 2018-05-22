import pytest
import os
import numpy as np
import parse_outputs as po


def test_parse():
    """Test the file parser
    """
    po.save_store_data('./ref/*o')
    
    ref_data = np.genfromtxt('./ref/ref_dep.csv', delimiter=',')
    obs_data = np.genfromtxt('depl_results.csv', delimiter=',')
    os.remove('depl_results.csv')
    
    print(ref_data)
    print(obs_data)
    
    for obs, exp in zip(obs_data, ref_data):
        assert obs == exp
