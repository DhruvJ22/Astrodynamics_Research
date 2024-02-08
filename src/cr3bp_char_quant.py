"""
Created on 10 Nov 2021 22:43:38 2022
Updaed on 20 Mar 2022

@author: Dhruv Jain, Multi-Body Dynamics Research Group, Purdue University
        dhruvj9922@gmail.com

Objective: This file contains a class with methods to compute various characteristic quantities
           for a Cirular Restircted Three Body Problem (CR3BP) model setup

           Features:
               1. Compute CR3BP 'mu', mass ratio of P1-P2 primary bodies
               2. Compute CR3BP 'l*', characterisitic lenght of P1-P2 primary bodies
               3. Comptue CR3BP 't*', charcterisitic time of P1-P2 primary bodies

Useful ideas:
1. To non-dimensioanlize a physical position vector: divide by l*
2. To dimensioanlize a [nd] position vector: multiply by l*
3. To non-dimensionalize a physical veloctiy vecotr, multiply by t*/l*
4. To dimensionalize a [nd] veloctiy vecotr, multiply by l*/t*
5. Position of P1 = [-mu, 0, 0] [nd] {x,y,z}
5. Position of P2 = [1-mu, 0, 0] [nd] {x,y,z}

Physical Qunatities obtained from: JPL’s ephemerides file de405.spk and https://ssd.jpl.nasa.gov/?planet_pos
"""
from src.cr3bp_quant_calc import mu_calc, tstar_calc, bodies_mu, bodies_dist

class sys_chars:
    """ Computes and stores the properties(mu, l*, t*) of a CR3BP system
        Can be extended to compute properties for BCR4BP model, HR4BP model etc 
    """

    def __init__(self, p1='Earth', p2='Moon'):
        """
        Constructor

        Parameters
        ----------
        p1 : string
           'body_i'. The deafult is Earth
        p2 : string
           'body_j'. The deafult is Moon
        """
        
        mu, lstar, tstar = self.bodies_char_compute(p1,p2)
        self._mu = mu
        self._lstar = lstar
        self._tstar = tstar
        self._p1 = p1
        self._p2 = p2
    
    # All the attributes are made private to make the constant and not be mistakenly changed
    @property
    def mu(self):
        """Mass ration of P1-P2 primary bodies in CR3BP"""
        return self._mu
    @property
    def lstar(self):
        """ Characterisitc Length """
        return self._lstar
    @property
    def tstar(self):
        """ Characterisitc Time """
        return self._tstar
    @property
    def p1(self):
        """ Primary body P1 """
        return self._p1
    @property
    def p2(self):
        """ Primary body P2 """
        return self._p2    
    
    def bodies_char_compute(self, p1, p2):
        """Calculates mu, dist[km], t* [nd] of the 'body_i - body_j system'
        Parameters
        ----------
        p1 : string
           'body_i'
        p2 : string
           'body_j'
    
        Returns
        -------
        mu: float
            mu2/(mu1+mu2), mu2<mu1
        dist: float, [nd]
            Non-dimensional distance between P1 and P2
        tstar: float, [nd]
            sqrt(dist^3/m*)
            Non-dimensional time of P1-P2 system
        """
        
        mu_pi = []
        mu_body = [bodies_mu(p1), bodies_mu(p2)]
        
        # Sort P1 and P2 based on their mu values
        if mu_body[0] >= mu_body[1]:  # p1 is the bigger primary
            mu_pi.append(mu_body[0])
            mu_pi.append(mu_body[1])
            bodies = p1 + p2
        else:  # p2 is the bigger primary
            mu_pi.append(mu_body[1])
            mu_pi.append(mu_body[0])
    
            # To create string to get distance
            bodies = p2 + p1
    
        mu = mu_calc(mu_pi[0], mu_pi[1])
    
        try:
            dist = bodies_dist(bodies)
        except KeyError:
            print(
                "KeyError-> Error in combination of bodies P1-P2, typo/DNE/combo not created"
            )
            return 0, 0, 0
        else:
            tstar = tstar_calc(mu_pi[0], mu_pi[1], dist)
            return mu, dist, tstar