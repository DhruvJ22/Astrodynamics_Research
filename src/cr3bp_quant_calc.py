def mu_calc(mu_1, mu_2):
    """Calculate mu of CR3BP
    Parameters
    ----------
    mu_1 = mu of P1
    mu_2 = mu of P2

    Returns
    -------
    mu: float, M2/(M1+M2)
        M1 and M2 are mass of Primary Bodies and M2<M1
    """
    return mu_2 / (mu_1 + mu_2)


def tstar_calc(mu_1, mu_2, dist):
    """Calculate t* of CR3BP
    Parameters
    ----------
    mu_1: float, mu of P1
    mu_2: float, mu of P2
     dist: float, [nd]
        Non-dimensional distance between P1 and P2
   
    Returns
    -------
    t*: float, [nd]
        sqrt(dist^3/m*)
        Non-dimensional time of P1-P2 system
    """
    return (dist**3 / (mu_1 + mu_2)) ** 0.5

def bodies_mu(body_name):
    """Returns mu value of various celestial bodies and distance between P1-P2 systems

    Returns
    -------
    mu_body: dict, float
    distances: dict, float
    """
    # Body values, G*M_body
    mu_body = {}  # km^3 kg^-1 s^-2
    mu_body['Sun'] = 132712440017.99
    mu_body['Moon'] = 4902.8005821478
    mu_body['Earth'] = 398600.4415

    mu_body['Mars'] = 42828.314258067 # Mars, GM
    mu_body['Jupiter'] = 126712767.8578 # Jupiter, GM
    mu_body['Saturn'] = 37940626.061137 # Saturn, GM
    mu_body['Uranus'] = 5794549.0070719 # Uranus, GM
    mu_body['Neptune'] = 6836534.0638793 # Neptune, GM
    mu_body['Pluto'] = 981.600887707 # Pluto, GM
    
    mu_body["Phobos"] = 0.0007112  # Phobos, GM
    mu_body["Titan"] = 8978.1382  # Titan, GM
    mu_body["Ganymede"] = 9887.834  # Ganymede, GM
    mu_body["Titania"] = 228.2  # Titania, GM
    mu_body["Triton"] = 1427.598  # Triton, GM
    mu_body["Charon"] = 102.30  # Charon, GM
 
    try:
        return mu_body[body_name]
    except KeyError:
        print("Body not in database or misspelled.")

def bodies_dist(bodies_pair_name):
    """Returns mu value of various celestial bodies and distance between P1-P2 systems

    Returns
    -------
    mu_body: dict, float
    distances: dict, float
    """
    distances = {}  # km, diistance between the two primaries
    distances["EarthMoon"] = 384400
    distances["SunEarth"] = 149600000

    distances["SunMars"] = 227944135
    distances["SunJupiter"] = 778279959
    distances["SunSaturn"] = 1427387908
    distances["SunUranus"] = 2870480873
    distances["SunNeptune"] = 4498337290
    distances["SunPluto"] = 5907150229

    distances["MarsPhobos"] = 9376
    distances["JupiterGanymede"] = 1070400
    distances["SaturnTitan"] = 1221865
    distances["UranusTitania"] = 436300
    distances["NeptuneTriton"] = 354759
    distances["PlutoCharon"] = 17536
    
    return distances[bodies_pair_name]