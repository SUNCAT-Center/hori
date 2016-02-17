#!/usr/bin/env python
"""Contains functions to calculate the change in free energy along the CO2
reduction pathway."""

import re
import numpy as np

from hori.common import countatoms


class Pathway:
    """Class to hold pathway and connectivity information, and to calculate
    the free energies along that pathway given binding (free) energies and
    chemical potentials. After initializing, add steps to the pathway with 
    the add_step function.

    Args:
        startingstep (str, optional):
        startingenergy (float, optional):
    """

    def __init__(self, startingstep='1', startingenergy=0.):
        self.startingstep = startingstep
        self.startingG = startingenergy
        self.steps = [['', self.startingstep, []]]

    def addstep(self, first, second, list):
        """Method to add a step to the pathway between the first and second
        states. Call this method with the state names and the list of
        changes between them. E.g., if state '1' to state '5' is::

            * + CO2 + (H+ + e-) -> *COOH

        then use this method with first='1',second='5', and list = 
        ['-s', '-CO2', '-pe', '+a_COOH']. In the list, reactants are
        preceded with a - and products with a +. Following the +/- is::

            s : bare surface
            pe : proton-electron pair
            a_<chemical formula> : adsorbate
            <chemical formula> : non-adsorbed species

        where the <>'s are left off the chemical formula. The info is

        sanity checked and added to an internal list of the form:
        
        .. code-block:: python

            [ [firststate,secondstate,reaction], ... ]

        Args:
            first (string)
            second (string)
            list (list)

        Raises:
            RuntimeError
        """
        # Convert the list of strings to a list of dictionaries
        reaction = []
        for list_string in list:
            reaction.append(self._parse_list_string(list_string))
        # Sanity check to make sure mass is conserved.
        surfaces = 0.
        elements = {'H': 0.}  # (for proton-electron pairs)
        for species in reaction:
            surfaces += species['count'] * species['adsorbed']
            elements['H'] += species['count'] * species['pe']
            fd = countatoms(species['formula'])  # formula dict
            for atom, subscript in fd.items():
                if atom not in elements.keys():
                    elements[atom] = 0.
                elements[atom] += species['count'] * subscript

        if sum(elements.values()) != 0 or surfaces != 0:
            raise RuntimeError('Inconsist mass balance in reaction.')
        # Sanity check to make sure the starting step exists.
        if first not in [step[1] for step in self.steps]:
            raise RuntimeError("Prior step does not exist.")
        # Sanity check to make sure exact same step doesn't already exist.
        for step in self.steps:
            if step[0] == first and step[1] == second:
                raise RuntimeError('Step exists -- cannot overwrite.')
        # Add to the internal list
        self.steps.append([first, second, reaction])

    def calculate_Gs(self, BG, mu, mu_pe):
        """Calculates the free energy along every step in the pathway, and
        returns this as a dictionary. Also store it in self.G.
        Takes as input BG, which is a dictionary of adsorbate free
        energies, mu, which is a dictionary of gas free energies, and
        mu_pe which is the free energy of a proton-electron pair (the last
        is a float, not a dictionary). These are usually taken from

        Args:
            BG (hori.thermo.AdsorbateThermodynamics.G)
            mu (hori.thermo.GasThermodynamics.G)
            mu_pe (hori.thermo.ProtonElectronThermodynamics.G)

        Raises:
            RuntimeError
        """
        G = {self.startingstep: self.startingG}
        for step in self.steps[1:]:
            priorstate, nextstate, rxn = step
            if priorstate not in G.keys():
                raise RuntimeError('Prior step does not exist. This method '
                                   'is set up assuming that all steps were '
                                   'added in such an order that the prior '
                                   'step occurs earlier in the list.')
            if nextstate not in G.keys():
                G[nextstate] = (G[priorstate] +
                                calculate_rxn_deltaG(rxn, BG, mu, mu_pe))
            else:
                testG = (G[priorstate] +
                         calculate_rxn_deltaG(rxn, BG, mu, mu_pe))
                if testG != G[nextstate]:
                    raise RuntimeError('Two different energies calculated '
                                       'for state %s' % nextstate)
        self.G = G
        return G

    def find_limiting_potential(self, path=None):
        """Uses the current values of the the free energy (stored in
        self.G; run self.calculate_Gs() first with mu_pe set to 0 V to
        establish this) to find the limiting potential and the step that
        limits the potential. This is defined as the largest uphill step
        -- after this is overcome, the pathway will be downhill.

        If the pathway is branched, then the path must be specified,

        similar to:
            path = ['1','28','4']

        Args:
            path (list, optional):

        Raises:
            RuntimeError
        """
        if path:
            steps = []
            for index in range(len(path) - 1):
                steps.append([path[index], path[index + 1]])
        else:
            # Check to be sure pathway is not branched.
            startsteps = []
            for item in self.steps:
                if item[0] in startsteps:
                    raise RuntimeError('path must be specified for branched'
                                       ' pathways.')
                startsteps.append(item[0])
            steps = self.steps[1:]

        limiting_potential = None
        for step in steps:
            dG = self.G[step[1]] - self.G[step[0]]
            if dG > limiting_potential:
                limiting_potential = dG
                limiting_step = [step[0], step[1]]
        limiting_potential *= -1.
        return limiting_potential, limiting_step

    def _parse_list_string(self, string):
        """Returns a dictionary with the count of atoms (and surfaces)
        in the strings in the list fed to Pathway.addstep().

        Args:
            string (string)

        Raises:
            RuntimeError
        """
        # Get sign.
        if string[0] == '-':
            count = -1
        elif string[0] == '+':
            count = +1
        else:
            raise RuntimeError('Misformatted string: %s' % string)
        string = string[1:]
        # Pull off any coefficient.
        pattern = re.compile('[0-9.]+')
        match = pattern.match(string)
        if match:
            count *= eval(match.group())
            string = string[match.end():]
        # Get if surface, adsorbate, or desorbed.
        if string == 's':  # Pure surface
            adsorbed = True
            chemicalformula = ''
            pe = 0.
        elif string == 'pe':  # proton-electron pair
            adsorbed = False
            chemicalformula = ''
            pe = 1.
        elif string[0] == 'a':  # Adsorbed chemical
            adsorbed = True
            chemicalformula = string[2:]
            pe = 0.
        else:  # Desorbed chemical
            adsorbed = False
            chemicalformula = string
            pe = 0.
        d = {'count': count,
             'adsorbed': adsorbed,
             'formula': chemicalformula,
             'pe': pe}
        return d


def calculate_rxn_deltaG(rxn, BG, mu, mu_pe):
    """Calculates delta-G of reaction for the reaction listed in rxn, given
    dictionaries of binding free energy (BG) and non-adsorbed species
    chemical potential (mu), as well as the chemical potential of a proton-
    electron pair at the current voltage (and pH, if applicable). rxn
    is in the format created by Pathway.addstep(). Note that in general,
    BG[''] will need to be defined; this is the clean slab's binding
    energy, which may be 0. or may be a large number, depending on the
    reference state chosen.

    Args:
        rxn (list)
        BG (hori.thermo.AdsorbateThermodynamics.G)
        mu (hori.thermo.GasThermodynamics.G)
        mu_pe (hori.thermo.ProtonElectronThermodynamics.G)
    """
    dG = 0.
    for species in rxn:
        if species['adsorbed'] == True:
            dG += species['count'] * BG[species['formula']]
        elif species['pe'] == 1.:
            dG += species['count'] * mu_pe
        else:
            dG += species['count'] * mu[species['formula']]
    return dG


def step2string(rxn):
    """Convert reaction step (as used in pathway) to a string.                                                                                                       

    Args:
        rxn (list): Description
    """
    spcs = rxn[2]
    reactants = r''
    products = r''
    # define mappings for species keywords to string representation
    adsorbed_str = {True: '*', False: ''}
    pe_str = {0.0: '', 1.0: ' H+ + e-'}
    for spc in spcs:
        s = spc['formula'] + adsorbed_str[spc['adsorbed']] + pe_str[spc['pe']] + ' + '
        if spc['count'] == -1:
            reactants += s
        elif spc['count'] == 1:
            products += s
    return reactants[:-2] + ' -> ' + products[:-2]
