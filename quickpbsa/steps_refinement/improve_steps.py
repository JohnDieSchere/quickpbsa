#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:59:08 2020

@author: johan
"""
import numpy as np

from .refinement_lowlevel import find_multstep
from .refinement_lowlevel import add_step

def improve_steps(data, jpos, means, variances, posbysic, maxmult = 5, lamb = 0.1, gamma0 = 0.5, \
                  multstep_fraction = 0.5, nonegatives = False, threshold = 1, maxadded = 5, \
                  splitcomb = 30000, combcutoff = 500000):
    ''' Refine the steps found by Kalafut Vischer. The algorithm will first find multiple steps, 
    then remove the steps in the inverse order in which they were found by KV and again find multiple steps.
    The first 2 steps are always left in place,  since they are crucial for the function of the algorithm.
    
    If no steps can be removed to yield an improved posterior,  steps will be added to improve the posterior.
    '''
    jpos_out = jpos
    success, sicmin, step_out = find_multstep(data, jpos, means, variances, maxmult, lamb, gamma0,
                                              multstep_fraction, nonegatives, threshold, splitcomb, combcutoff)
    nsteps_in = len(jpos)
    nsteps = nsteps_in
    counter = 0
    for I in range(1, len(jpos) - 2):
        ind_keep = np.in1d(jpos, posbysic[: - I], assume_unique = True)
        ind_keep[:2] = True # not allowed to remove the first 2 steps
        if np.sum(ind_keep) < nsteps:
            jpos = jpos[ind_keep]
            nsteps = len(jpos)
            means = np.array([np.mean(ar) for ar in np.split(data, jpos)])
            variances = np.array([np.var(ar) for ar in np.split(data, jpos)])
            success, sic, step = find_multstep(data, jpos, means, variances, maxmult, lamb, gamma0, multstep_fraction,\
                                               nonegatives, threshold, splitcomb, combcutoff)
            if sic < sicmin:
                sicmin = sic
                step_out = step
                jpos_out = jpos
    means = np.array([np.mean(ar) for ar in np.split(data, jpos_out)])
    variances = np.array([np.var(ar) for ar in np.split(data, jpos_out)])
    nsteps = len(jpos_out)
    if nsteps==nsteps_in:
        counter=0
        # add steps if no steps removed
        for I in range(maxadded):
            sicmin, jpos_out, step_out, means, variances = add_step(data, jpos_out, means, variances,
                                                                    step_out, lamb, gamma0)
            if len(step_out) == nsteps:
                # break if adding the last 2 steps did not help
                counter += 1
                if counter > 2:
                    break
            else:
                # added a step - >  continue
                counter = 0
                nsteps = len(step_out)
    jpos_out = np.array(jpos_out)
    step_out = np.array(step_out)
    return success, sicmin, jpos_out, step_out