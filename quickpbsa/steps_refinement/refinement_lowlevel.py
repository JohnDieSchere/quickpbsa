#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:59:08 2020

@author: johan
"""
import numpy as np
from sympy.utilities.iterables import multiset_permutations
import itertools as it
from scipy.special import factorial

def posterior_single(data, jpos, means, variances, df, i_in=0, lamb=0.1, gamma0=0.5):
    ''' single calculation of the posterior criterion from the Presse paper.
    '''
    # fluorophore and bg parameters from KV result
    mb = means[0]
    vb = variances[0]
    mf = means[1] - means[0]
    vf = variances[1] - variances[0]
    if vf < 0:
        vf = variances[1]
    i = np.cumsum([0] + list(df))
    if np.all(i >= 0):
        nphi = np.diff([0] + list(jpos) + [len(data)])
        varphi = i*vf + vb
        K = len(jpos)
        m = sum(map(abs, df))
        stp, ct = np.unique(np.abs(df), return_counts = True)
        sic = sum(nphi*np.log(varphi) + np.array(list(map(sum, (np.split(data, jpos) - i*mf - mb)**2)))/varphi)
        sic += 2*( - K*np.log(lamb) - np.log(factorial(m - K) * factorial(K) / factorial(m - 1)) + np.sum(np.log(factorial(ct))))
        sic += 2*gamma0*(m - K + 1)/K + np.log(m - K + 2) + np.log(m - K + 1) - np.log(m - K + 2 - (m - K + 1) * np.exp(-gamma0 / K))
    else:
        sic = np.float('Inf')
    return sic



def find_multstep(data, jpos, means, variances, maxmult=5, lamb=0.1, gamma0=0.5, multstep_fraction=0.5, 
                  nonegatives=False, threshold=1, splitcomb=30000, combcutoff=1000000):
    ''' find multiple steps according to the posterior criterion. The number of combinations get get extremely high very quickly,  
        which increases runtime. 
        There are 3 options to reduce the number of combinations:
        multstep_fraction (default 0.5): only ~ len(jpos)*multstep_fraction double steps are considered
        nonegatives (default False): no negative (i.e. fluorophores turning back on) double steps are considered
        threshold: steps below laststep*threshold in the mean difference are not considered as double steps
    '''
    success = True
    # last fluorophore and bg stats
    mb = means[0]
    vb = variances[0]
    mf = means[1] - means[0]
    vf = variances[1] - variances[0]
    if vf < 0:
        vf = variances[1]
    threshold = mf*threshold/2
    # initial steps
    stepsigns = np.sign(np.diff(means)).astype(int)
    step_out = stepsigns
    sicmin = posterior_single(data, jpos, means, variances, stepsigns)  
    # multstep = 1 means double steps and so on
    for multstep in range(1, maxmult):
        threshold = threshold*(multstep + 1)/multstep
        if multstep > 1:
            ind_mult = np.where(np.abs(step_out) == multstep)[0]
        else:
            # for no additional filtering (no need to consider the first step)
            multsel = np.arange(len(jpos)) > 0
            if nonegatives:
                # no negative double steps
                multsel = multsel*(stepsigns > 0)
            if threshold is not None:
                # only steps where mean difference is at least threshold
                multsel = multsel*(np.abs(np.diff(means)) >= threshold)
            ind_mult = np.where(multsel)[0]
        N_comb = (multstep + 1)**len(ind_mult)
        if N_comb > combcutoff:
            success = False
            break
        else:
            if len(ind_mult) > 0:
                # there are multiple steps to try out,  loop through combinations
                stepcomb = list(it.combinations_with_replacement(list(range(1, multstep + 2)), len(ind_mult)))
                # filter out old combinations and ones with too few fluorophores
                stepcomb = [comb for comb in stepcomb if np.sum(np.array(comb) == multstep + 1) > 0]
                stepcomb = [comb for comb in stepcomb if np.sum(np.array(comb)) > np.sum(np.abs(step_out[ind_mult])) - 3]
                if multstep == 1 and multstep_fraction < 1:
                    # only allow multstep_fraction of steps to be double steps
                    stepmax = int(np.ceil(len(jpos)*multstep_fraction)) + 1
                    stepcomb = stepcomb[:stepmax]
                for comb in stepcomb:
                    improved = False
                    # sic components independent of step placement
                    stp = np.hstack((np.delete(np.abs(step_out), ind_mult), comb))
                    m = np.sum(stp)
                    K = len(jpos)
                    stp, ct = np.unique(stp, return_counts = True)
                    sic_in = 2*( - K*np.log(lamb) - np.log(factorial(m - K)*factorial(K)/factorial(m - 1)) + np.sum(np.log(factorial(ct))))
                    sic_in += 2*gamma0*(m - K + 1)/K + np.log(m - K + 2) + np.log(m - K + 1) - np.log(m - K + 2 - (m - K + 1)*np.exp( - gamma0/K))
                    #combinations array
                    cc = np.array(list(multiset_permutations(comb)))
                    if len(jpos) != len(ind_mult):
                        # need to add rows of ones for not considered step positions
                        combarray = np.ones([np.size(cc, 0), len(jpos)], dtype = int)
                        combarray[:, ind_mult] = cc
                    # restrict combarray by threshold (for steps beyond 2)
                    if multstep > 1:
                        exclude = (np.diff(means) <= (threshold*(multstep + 1)))
                        sel = np.sum(combarray[:, exclude] == multstep + 1, 1).astype(bool)
                        combarray = combarray[np.invert(sel), :]
                        if np.sum(np.invert(sel)) == 0:
                            break
                    i_complete = np.cumsum(np.hstack((np.zeros([np.size(combarray, 0), 1], dtype = int), combarray*stepsigns)), 1)
                    # go through combinations in packages of splitcomb
                    for i in np.split(i_complete, np.arange(splitcomb, np.size(i_complete, 0), splitcomb)):
                        # filter out negative fluorophores
                        sel = np.invert(np.sum(i < 0, 1).astype(bool))
                        i = i[sel, :]
                        if np.size(i, 0) > 0:
                            # sic calculation
                            steploc = np.hstack((0, jpos, len(data)))
                            nphi = np.diff(steploc)
                            varphi = i*vf + vb
                            diffar = np.zeros(np.shape(i))
                            for I in range(len(steploc) - 1):
                                # difference array
                                dar = np.tile(data[steploc[I]:steploc[I + 1]], [np.size(i, 0), 1])
                                diffar[:, I] = np.sum((dar.T - mf*i[:, I] - mb)**2, 0)
                            sic = sic_in + np.sum(nphi*np.log(varphi) + diffar/varphi, 1)
                            if np.min(sic) < sicmin:
                                improved = True
                                sicmin = np.min(sic)
                                step_out = np.diff(i[np.argmin(sic), :])
                    if multstep == 1 and not improved:
                        # no need to try more double steps
                        break
    return success, sicmin, step_out



def add_step(data, jpos, means, variances, steps, lamb = 0.1, gamma0 = 0.5):
    ''' Try adding a negative or positive step to the steps already found. Only the part of the trace
    before the penultimate step will be considered.
    '''
    mb = means[0]
    vb = variances[0]
    mf = means[1] - means[0]
    vf = variances[1] - variances[0]
    if vf < 0:
        vf = variances[1]
    sz = len(data)
    
    # initial sic
    sicmin = posterior_single(data, jpos, means, variances, steps)
    stp, ct = np.unique(np.abs(steps), return_counts = True)
    # add one to steps with zero occupancy
    ct[0] += 1
    step_out = steps
    jpos_out = jpos
    
    # build arrays for squared diff,  variance and length
    i_in = np.cumsum(np.hstack((0, steps)))
    steploc = np.hstack((0, jpos, sz))
    diffar = np.tile(np.array(list(map(sum, (np.split(data, jpos) - i_in*mf - mb)**2))), [sz - jpos[1], 1])
    diffar = np.hstack((diffar, np.zeros([sz - jpos[1], 1])))
    varphi = np.tile(np.hstack((i_in*vf + vb, 0)), [sz - jpos[1], 1])
    nphi = np.tile(np.hstack((np.diff(steploc), 0)), [sz - jpos[1], 1])
    
    sics = []
    steppos = []
    for stp in (-1, 1):
        for I in range(2, len(steploc) - 1):
            # length
            stepl = int(steploc[I + 1] - steploc[I])
            if I < len(steploc) - 2:
                # shift values after current range by one and recalculate with i  + / -  1
                try:
                    dd = np.array(list(map(sum, (np.split(data[steploc[I + 1]:], jpos[I + 1:] - steploc[I + 1]) - (i_in[I + 1:] + stp)*mf - mb)**2)))
                except ValueError:
                    dd = np.array(list(map(sum, (np.split(data[steploc[I + 1]:], jpos[I + 1:] - steploc[I + 1]) - np.expand_dims((i_in[I + 1:] + stp)*mf - mb, axis = 1))**2)))
                diffar[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I + 2:] = np.tile(dd, [stepl, 1])
                varphi[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I + 2:] = np.tile((i_in[I + 1:] + stp)*vf + vb, [stepl, 1])
                nphi[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I + 2:] = np.tile(np.diff(steploc[I + 1:]), [stepl, 1])
            
            # diffs under the diagonal (before added step)
            dar = np.tile(data[steploc[I]:steploc[I + 1]], [stepl, 1])
            dar[np.triu(np.ones([stepl, stepl]), 1).astype(bool)] = i_in[I]*mf + mb # diff values will be zeros
            diffar[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I] = np.sum((dar - i_in[I]*mf - mb)**2, 1)
            
            # diffs above diagonal (after added step)
            dar = np.tile(data[steploc[I]:steploc[I + 1]], [stepl, 1])
            dar[np.tril(np.ones([stepl, stepl]), 0).astype(bool)] = (i_in[I] + stp)*mf + mb
            diffar[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I + 1] = np.sum((dar - (i_in[I] + stp)*mf - mb)**2, 1)
            
            # varphi for the current range
            varphi[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I] = i_in[I]*vf + vb
            varphi[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I + 1] = (i_in[I] + stp)*vf + vb
            
            # length for current range
            nphi[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I] = np.arange(stepl) + 1
            nphi[steploc[I] - jpos[1]:steploc[I + 1] - jpos[1], I + 1] = stepl - np.arange(stepl) - 1
        
        # SIC calculation
        m = np.sum(np.abs(steps)) + 1
        K = len(jpos) + 1
        # This sometimes throws warnings invalid value in log
        sic = np.sum(nphi*np.log(varphi) + diffar/varphi, 1)
        sic += 2*( - K*np.log(lamb) - np.log(factorial(m - K)*factorial(K)/factorial(m - 1)) + np.sum(np.log(factorial(ct))))
        sic += 2*gamma0*(m - K + 1)/K + np.log(m - K + 2) + np.log(m - K + 1) - np.log(m - K + 2 - (m - K + 1)*np.exp( - gamma0/K))
        
        # exclude positions around existing steps and edges
        jignore = np.ravel(np.tile(jpos, [5, 1]).T + np.arange( - 2, 3)) - jpos[1]
        jignore = jignore[(jignore > 0)*(jignore < len(sic) - 2)]
        ind = np.arange(len(sic) - 2)
        ind = np.delete(ind, jignore)
        
        sic = sic[ind]
        sics.append(np.min(sic))
        steppos.append(ind[np.argmin(sic)] + jpos[1] + 1)
    
    # nanmin in case of negative fluorophores,  maybe can be done more elegantly
    # sometimes warns about all nan axis
    if np.nanmin(sics) < sicmin:
        # found new step,  add it
        sicmin = np.nanmin(sics)
        newstep = steppos[np.nanargmin(sics)]
        stepsign = np.nanargmin(sics)*2 - 1
        insertpos = np.max(np.where(jpos < newstep)[0]) + 1
        step_out = np.hstack((steps[:insertpos], stepsign, steps[insertpos:]))
        jpos_out = sorted(np.hstack((jpos, newstep)))
        
        # calculate variances and means
        variances = [np.var(ar) for ar in np.split(data, jpos_out)]
        means = [np.mean(ar) for ar in np.split(data, jpos_out)]
        
    return sicmin, jpos_out, step_out, means, variances