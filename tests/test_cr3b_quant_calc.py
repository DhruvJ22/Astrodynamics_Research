from src.cr3bp_quant_calc import mu_calc, tstar_calc, bodies_mu, bodies_dist

def test_mu_calc():
    assert mu_calc(1.0,1.0) == 0.5

def test_tstar_calc():
    assert tstar_calc(2.0,6.0,2.0) == 1.0

def test_bodies_mu():
    assert bodies_mu('Earth') == 398600.4415

def test_bodies_dist():
    assert bodies_dist('EarthMoon') == 384400
    