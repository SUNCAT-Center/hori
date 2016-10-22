#!/usr/bin/env python
"""Contains common elements, such as reference energies."""


def countatoms(formula):
    """Given an atomic formula, e.g. CH2O, returns a dictionary of
    atomic counts, e.g, d = {'C':1,'H':2,'O':1}.

    Args:
        formula (string)
    """
    d = {}
    starts = []  # indices of string formula that start a new element
    for index, letter in enumerate(formula):
        if letter.isupper():
            starts.append(index)
    starts.append(len(formula))  # (gives searchable range for last el.)
    for index in range(len(starts) - 1):
        name = ''
        number = ''
        for letter in formula[starts[index]:starts[index + 1]]:
            if letter.isalpha():
                name = name + letter
            elif letter.isdigit():
                number = number + letter
        if len(number) == 0:
            count = 1
        else:
            count = eval(number)
        if name in d.keys():
            d[name] += count
        else:
            d[name] = count
    return d


def calculate_reference_energy(formula):
    """For a formula, such as CH2OH, counts the atoms and makes a
    reference energy (eV) based on the atomic energies in
    reference_energies.

    Args:
        formula (string)
    """
    referenceenergy = 0.
    counts = countatoms(formula)
    for key in counts.keys():
        referenceenergy += counts[key] * reference_energies[key]
    return referenceenergy

# The below energies are just used for reporting of data, not
# in internal calculations.
reference_energies = {
    'C': -155.456,  # graphene
    'H': 0.5 * (-32.0295),  # 1/2 H2
    'O': -469.7509 + 32.0295,  # H2O-H2
}


"""For convenience, all of the surfaces which have calculations
are listed in all_surfaces."""
all_surfaces = [
    'Rh-fcc211',
    'Al-fcc211',
    'Pd-fcc211',
    'Ag-fcc211',
    'Ir-fcc211',
    'Pt-fcc211',
    'Ni-fcc211',
    'Cu-fcc211',
    'Pb-fcc211',
    'Au-fcc211',
    'Tl-fcc211',
    'Cd-fcc211', ]


def create_full_pathway():
    """Returns the full CO2 reduction pathway that I have commonly looked
    at."""
    from hori.pathway import Pathway
    fp = Pathway(startingstep='1', startingenergy=0.)
    fp.addstep('1', '5', ['-s', '-CO2', '-pe', '+a_OCHO'])
    fp.addstep('5', '6', ['-a_OCHO', '-pe', '+s', '+HCOOH'])
    fp.addstep('5', '11', ['-a_OCHO', '-pe', '+a_OCH2O'])
    fp.addstep('11', '12', ['-a_OCH2O', '-s', '-pe', '+a_OCH3', '+a_O'])
    fp.addstep('12', '14', ['-a_OCH3', '-a_O', '-pe', '+a_OCH3', '+a_OH'])
    fp.addstep('14', '16', ['-a_OCH3', '-a_OH', '-pe', '+s', '+a_OCH3', '+H2O'])
    fp.addstep('16', '18', ['-a_OCH3', '-pe', '+a_O', '+CH4'])
    fp.addstep('18', '20', ['-a_O', '-pe', '+a_OH'])
    fp.addstep('20', '21', ['-a_OH', '-pe', '+s', '+H2O'])
    fp.addstep('1', '28', ['-s', '-CO2', '-pe', '+a_COOH'])
    fp.addstep('28', '3', ['-a_COOH', '-pe', '+a_CO', '+H2O'])
    fp.addstep('28', '4', ['-a_COOH', '-pe', '+CO', '+H2O', '+s'])
    fp.addstep('28', '6', ['-a_COOH', '-pe', '+HCOOH', '+s'])
    fp.addstep('28', '2', ['-a_COOH', '-s', '+a_CO', '+a_OH'])
    fp.addstep('3', '22', ['-a_CO', '-pe', '+a_COH'])
    fp.addstep('3', '23', ['-a_CO', '-pe', '+a_CHO'])
    fp.addstep('23', '25', ['-a_CHO', '-pe', '+a_CH2O'])
    fp.addstep('23', '24', ['-a_CHO', '-pe', '+a_CHOH'])
    fp.addstep('22', '24', ['-a_COH', '-pe', '+a_CHOH'])
    fp.addstep('25', '26', ['-a_CH2O', '-pe', '+a_CH2OH'])
    fp.addstep('24', '26', ['-a_CHOH', '-pe', '+a_CH2OH'])
    fp.addstep('25', '16', ['-a_CH2O', '-pe', '+a_OCH3'])
    fp.addstep('26', '27', ['-a_CH2OH', '-pe', '+CH3OH', '+s'])
    fp.addstep('22', '29', ['-a_COH', '-pe', '+a_C', '+H2O'])
    fp.addstep('29', '30', ['-a_C', '-pe', '+a_CH'])
    fp.addstep('24', '30', ['-a_CHOH', '-pe', '+H2O', '+a_CH'])
    fp.addstep('30', '31', ['-a_CH', '-pe', '+a_CH2'])
    fp.addstep('26', '31', ['-a_CH2OH', '-pe', '+H2O', '+a_CH2'])
    fp.addstep('31', '32', ['-a_CH2', '-pe', '+a_CH3'])
    fp.addstep('32', '21', ['-a_CH3', '-pe', '+CH4', '+s'])
    fp.addstep('1', '33', ['-s', '-pe', '+a_H'])
    fp.addstep('33', '34', ['-a_H', '-pe', '+H2', '+s'])
    fp.addstep('16', '27', ['-a_OCH3', '-pe', '+s', '+CH3OH'])
    fp.addstep('28', '43', ['-a_COOH', '-pe', '+a_COHOH'])
    fp.addstep('43', '22', ['-a_COHOH', '-pe', '+a_COH', '+H2O'])
    fp.addstep('5', '44', ['-a_OCHO', '-s', '+a_CHO', '+a_O'])
    fp.addstep('44', '45', ['-a_CHO', '-a_O', '-pe', '+a_CHO', '+a_OH'])
    fp.addstep('45', '23', ['-a_CHO', '-a_OH', '-pe', '+a_CHO', '+s', '+H2O'])
    fp.addstep('3', '46', ['-a_CO', '-s', '+a_C', '+a_O'])
    fp.addstep('46', '37', ['-a_C', '-a_O', '-pe', '+a_C', '+a_OH'])
    fp.addstep('2', '3', ['-a_CO', '-a_OH', '-pe', '+H2O', '+s', '+a_CO'])
    fp.addstep('22', '37', ['-a_COH', '-s', '+a_C', '+a_OH'])
    fp.addstep('37', '29', ['-a_C', '-a_OH', '-pe', '+a_C', '+s', '+H2O'])
    fp.addstep('23', '38', ['-a_CHO', '-s', '+a_CH', '+a_O'])
    fp.addstep('38', '47', ['-a_CH', '-a_O', '-pe', '+a_CH', '+a_OH'])
    fp.addstep('25', '48', ['-a_CH2O', '-s', '+a_CH2', '+a_O'])
    fp.addstep('48', '39', ['-a_CH2', '-a_O', '-pe', '+a_CH2', '+a_OH'])
    fp.addstep('47', '30', ['-a_CH', '-a_OH', '-pe', '+a_CH', '+s', '+H2O'])
    fp.addstep('16', '49', ['-a_OCH3', '-s', '+a_O', '+a_CH3'])
    fp.addstep('49', '50', ['-a_CH3', '-a_O', '-pe', '+a_CH3', '+a_OH'])
    fp.addstep('39', '31', ['-a_CH2', '-a_OH', '-pe', '+a_CH2', '+s', '+H2O'])
    fp.addstep('50', '32', ['-a_CH3', '-a_OH', '-pe', '+a_CH3', '+s', '+H2O'])
    fp.addstep('25', '40', ['-a_CH2O', '+0.5C2H4', '+a_O'])
    fp.addstep('40', '41', ['-a_O', '-pe', '+a_OH'])
    fp.addstep('41', '42', ['-a_OH', '-pe', '+s', '+H2O'])
    return fp


def hbond_dict():
    """Returns a dictionary of hydrogen bond stabilizations, in eV."""
    OH = -0.50  # OH directly on the surface
    ROH = -0.25  # a 'floppy' OH group
    CO = -0.1  # carbon monoxide
    d = {'COOH': ROH,
         '': 0.,
         'OCHO': 0.,
         'CO': CO,
         'CO2': 0.,
         'CH2O': 0.,
         'CH3O': 0.,
         'CHO': CO,
         'COCOH': CO + ROH,
         'OCCOH': CO + ROH,
         'CH2O': 0.,
         'CH3O': 0.,
         'OCH3': 0.,
         'O': 0.,
         'OH': OH,
         'S': 0.,
         'SH': 0.,
         'H': 0.,
         'COH': ROH,
         'H2O': ROH,
         'C': 0.,
         'CC': 0.,
         'CH': 0.,
         'CH2': 0.,
         'CH2_CH2': 0.,
         'CH2_CH2_TS': 0.,
         'C2H4': 0.,
         'CH3': 0.,
         'CHOH': ROH,
         'COHOH': ROH,
         'OCH2O': 0.,
         'CH2OH': ROH,
         'OCHCH2': 0.,
         'OCH2CH': 0.,
         'OHCHCH': ROH,
         'OCHCHO': 0.,
         'OCHCOH': CO + ROH,
         'OCH2CHO': 0.,
         'OCH2CH2O': 0.,
         'CHO_CHO': 2 * CO,
         'CH2O_CHO': CO,
         'CH2O_CH2O': 0,
         'OCH2CHOH': ROH,
         'OCHCHOH': ROH,
         'OCHCH': 0.,
         'OCHCH2OH': ROH,
         'OCH2CHOH': ROH,
         'OHCHCHOH': 2 * ROH,
         'OCHCH': 0.,
         'CHOCHO_OCHCHO_TS': 2 * CO,
         'CH2OCHO_OCH2CHO_TS': CO,
         'CH2OCH2O_OCH2CH2O_TS': 0,
         'CHOCO_TS': 2 * CO,
         'COCO_OCCO_TS': 2 * CO,
         'O_O': 0.,
         'CO_CO': 2 * CO,
         'COOH_COOH': 2 * ROH,
         'COHCOH': 2 * ROH,
         'OCH2CH2OH': ROH,
         'OHCH2CH2': ROH,
         'OCH2CH2': 0.0,
         'OH_O': OH,
         'OH_OH': 2 * OH,
         'CHOCO': 2 * CO,
         'CHO_CO': 2 * CO,
         'OCCO': 2 * CO,
         'O2H2HCCH2OH': 3 * ROH,
         'O2H2HCCHO2H2': 4 * ROH,
         'O2H2HCCH': 2 * ROH,
         'O2H2HCCH2': 2 * ROH,
         'OHHCCH2OH': 2 * ROH,
         'OHHCCH2': 1 * ROH,
         'HCCH2': 0.,
         'OHH2CCH2': 1 * ROH,
         'O2H2HCCHOH': 3 * ROH,
         'OHHCCHOH': 2 * ROH,
         'OHHCCH': ROH,
         'HCCH': 0.,
         'OHHCCH3': ROH,
         'O2H2HCCH3': 2 * ROH,
         'HCCH2OH': ROH,
         'HCCH3': 0.,
         'OCCHO': 2 * CO,
         'OCCH2O': CO,
         'OCHCH3': 0.,
         'OHCHCH2': ROH,
         'OHCHCH3': ROH,
         'CHCH2': 0.,
         'OHCH2CH': ROH,
         'NNH': 0.,
         'NNH2': 0.,
         'NHNH': 0.,
         'NHNH2': 0.,
         'NH': 0.,
         'N': 0.,
         'N2': 0.,
         'NH2': 0.,
         'NH3': 0.,
         'NH_H': 0.,
         'NH-H': 0.,
         'NH2_H': 0.,
         'NH2-H': 0.,
         'NH3+H+H+H': 0.,
         'NH3+H+H': 0.,
         'NH3+H': 0.,
         }
    return d


def hbond_dict_none():
    """Returns a dictionary of hydrogen bond stabilizations, in eV."""
    OH = 0.0  # OH directly on the surface
    ROH = 0.0  # a 'floppy' OH group
    CO = 0.0  # carbon monoxide
    d = {'COOH': ROH,
         'CO2': 0.,
         '': 0.,
         'OCHO': 0.,
         'CO': CO,
         'CHO': CO,
         'OCCOH': CO + ROH,
         'COCOH': CO + ROH,
         'CH2O': 0.,
         'CH3O': 0.,
         'OCH3': 0.,
         'O': 0.,
         'OH': OH,
         'S': 0.,
         'SH': 0.,
         'H': 0.,
         'COH': ROH,
         'C': 0.,
         'CC': 0.,
         'CH': 0.,
         'CH2': 0.,
         'CH2_CH2': 0.,
         'CH2_CH2_TS': 0.,
         'C2H4': 0.,
         'CH3': 0.,
         'CHOH': ROH,
         'COHOH': ROH,
         'OCH2O': 0.,
         'CH2OH': ROH,
         'OCHCH2': 0.,
         'OCHCHO': 0.,
         'OCH2CHO': 0.,
         'OCH2CH2O': 0.,
         'CHO_CHO': 2 * CO,
         'CH2O_CHO': CO,
         'CH2O_CH2O': 0,
         'CHOCHO_OCHCHO_TS': 2 * CO,
         'CH2OCHO_OCH2CHO_TS': CO,
         'CH2OCH2O_OCH2CH2O_TS': 0,
         'CHOCO_TS': 2 * CO,
         'COCO_OCCO_TS': 2 * CO,
         'O_O': 0.,
         'CO_CO': 2 * CO,
         'COOH_COOH': 2 * ROH,
         'COHCOH': 2 * ROH,
         'OCH2CH2OH': ROH,
         'OCH2CH2': 0.,
         'OHCH2CH2': ROH,
         'OH_O': OH,
         'OH_OH': 2 * OH,
         'CHOCO': 2 * CO,
         'CHO_CO': 2 * CO,
         'OCCO': 2 * CO,
         'O2H2HCCH2OH': 3 * ROH,
         'O2H2HCCHO2H2': 4 * ROH,
         'O2H2HCCH': 2 * ROH,
         'O2H2HCCH2': 2 * ROH,
         'OHHCCH2OH': 2 * ROH,
         'OHHCCH2': 1 * ROH,
         'HCCH2': 0.,
         'OHH2CCH2': 1 * ROH,
         'O2H2HCCHOH': 3 * ROH,
         'OHHCCHOH': 2 * ROH,
         'OHHCCH': ROH,
         'HCCH': 0.,
         'OHHCCH3': ROH,
         'O2H2HCCH3': 2 * ROH,
         'HCCH2OH': ROH,
         'NH3': 0.,
         'HCCH3': 0.
         }
    return d
