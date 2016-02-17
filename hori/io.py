#!/usr/bin/env python
"""Functions to read in or write out data."""

import os
import pickle
import numpy as np

from ase import Atom, Atoms

from hori import __path__ as horipath
datapath = os.path.join(horipath[0],'data')

def set_directory(path):
    """Use to set a directory of pickle files other than the current. The
    directory specified must have subdirectories 'electronic-energies',
    and/or 'generic-vibrations', depending on the data to be accessed.
    """
    global datapath
    datapath = path

def get_directory():
    """Returns the internal value of datapath."""
    return datapath

def electronicenergy(surface, adsorbate=None, verbose=True):
    """Reads in the electronic energy from the calculated data set. If the
    adsorbate is None or '', then returns the energy of the clean slab,
    if the energy doesn't exist in the data set, then returns nan. Alerts
    the user to this fact if verbose=True."""
    if adsorbate == '': adsorbate = None
    if adsorbate:
        filename = '_'.join([adsorbate,surface])
    else:
        filename = surface
    filename = os.path.join(datapath,'electronic-energies',filename)
    if os.path.isfile(filename):
        f = open(filename)
        d = pickle.load(f)
        f.close()
        if d.has_key('remark'):
            if d['remark']:
                print('Message from pickle %s:\n %s' % 
                      (filename, d['remark']))
        energy = d['electronic energy']
    else:
        if verbose:
            print('No file found at %s' % filename)
            print os.listdir(os.path.split(filename)[0])
        energy = np.nan
    return energy

def gasdata(species):
    """Reads in the electronic energy and vibrations from the calculated
    data set."""
    f = open(os.path.join(datapath,'electronic-energies',species))
    d = pickle.load(f)
    if d.has_key('remark'):
        if d['remark']:
            print('Message from pickle %s:\n %s' % (filename, d['remark']))
    f.close()
    # If using a new (version 2) gasenergypickle, then convert the
    # atomslist to an atoms object
    if d.has_key('version'):
        if d['version'] >= 2:
            d['atoms'] = atomsfromlist(d['atomslist'])
    return d

def atomsfromlist(atomslist):
    """Takes in a list of atomic symbols and coordinates, as in 
    [atom1, atom2, ...] where atomX = (symbol, (x,y,z)), and symbol is the
    atomic symbol (e.g. "Na") and x,y,z is the position, in Angstroms, of
    the atom. Returns an ASE atoms object."""
    atoms = Atoms()
    for atom in atomslist:
        atoms.append(Atom(atom[0], atom[1]))
    return atoms

def genericvibrations(adsorbate):
    """Reads in the generic vibrations used for adsorbates in general."""
    filename = os.path.join(datapath,'generic-vibrations',adsorbate)
    f = open(filename)
    d = pickle.load(f)
    if d.has_key('remark'):
        if d['remark']:
            print('Message from pickle %s:\n %s' % (filename, d['remark']))
    f.close()
    vib_energies = d['vibrations']
    # Check to see if any frequencies are imaginary, and report to user.
    if sum(np.iscomplex(vib_energies)) > 0:
        print('WARNING: Imaginary frequencies encountered for %s.' %
              adsorbate)
    return vib_energies

def specificvibrations(adsorbate, surface):
    """Reads in the vibrations for an adsorbate on a specific surface."""
    if surface is None:
        raise RuntimeError('A surface must be specified for specific '
                           'vibrations.')
    filename = '_'.join([adsorbate, surface])
    filename = os.path.join(datapath, 'specific-vibrations',
                            filename)
    f = open(filename)
    d = pickle.load(f)
    f.close()
    if d.has_key('remark'):
        if d['remark']:
            print('Message from pickle %s:\n %s' % (filename, d['remark']))
    if d.has_key('vibrations'):
	realvibs=[]
    	for vib in d['vibrations']:
        	if np.real(vib)>0:
            		realvibs.append(float(np.real(vib)))
    	vib_energies = np.array(realvibs)
        #vib_energies = d['vibrations']
    else:
        raise RuntimeError('No vibrations for %s on %s.' % (adsorbate,
                                                            surface))
    # Check to see if any frequencies are imaginary, and report to user.
    if sum(np.iscomplex(vib_energies)) > 0:
        print('WARNING: Imaginary frequencies encountered for %s on %s.' %
              (adsorbate, surface))
    return vib_energies

def ignoreimaginaryvibrations(adsorbate):
    """Reads in the generic vibrations used for adsorbates in general.
    Calculates only real part."""
    filename = os.path.join(datapath,'generic-vibrations',adsorbate)
    f = open(filename)
    d = pickle.load(f)
    if d.has_key('remark'):
        if d['remark']:
            print('Message from pickle %s:\n %s' % (filename, d['remark']))
    f.close()
    realvibs=[]
    for vib in d['vibrations']:
        if np.real(vib)>0:
            realvibs.append(float(np.real(vib)))
    vib_energies = np.array(realvibs)
    return vib_energies
    # Check to see if any frequencies are imaginary, and report to user.
    if sum(np.iscomplex(vib_energies)) > 0:
        print('WARNING: Imaginary frequencies encountered for %s on %s.' %
              (adsorbate, surface))

    
