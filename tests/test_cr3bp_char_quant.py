from src.cr3bp_char_quant import sys_chars
import numpy as np

def test_bodies_char_compute():
    # Create class instance and call method
    sys_p1p2 = sys_chars('Earth','Moon')

    np.testing.assert_allclose(sys_p1p2.mu, 0.012150585350562453, rtol=1e-6)  # Allow difference up to 0.00001
    assert sys_p1p2.lstar == 384400
    np.testing.assert_allclose(sys_p1p2.tstar, 375190.25889262, rtol = 1.e-6)


