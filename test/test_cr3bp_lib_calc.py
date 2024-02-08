from src.cr3bp_lib_calc import lib_pt_loc
from src.cr3bp_char_quant import sys_chars
import numpy as np

def test_lib_pt_loc():
    sys_p1p2 = sys_chars("Earth","Moon")
    lib_loc = lib_pt_loc(sys_p1p2)

    expected_lib_loc_EM = np.array([[ 0.83691513,  0.        ,  0.        ],
       [ 1.15568216,  0.        ,  0.        ],
       [-1.00506265,  0.        ,  0.        ],
       [ 0.48784941,  0.8660254 ,  0.        ],
       [ 0.48784941, -0.8660254 ,  0.        ]])

    np.testing.assert_allclose(lib_loc, expected_lib_loc_EM, rtol=1e-8)  # Allow element-wise difference up to 0.00001