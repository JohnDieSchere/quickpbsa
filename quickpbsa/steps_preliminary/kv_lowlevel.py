#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import numpy as np
import math


def sic_single(data, steps):
    '''
    calculates the Schwartz Inforamtion Criterion for the data and step locations given in jtest
    '''
    N = len(data)
    splar = np.split(data, steps)
    SIC = (len(steps) + 2) * math.log(N) + N * math.log(math.fsum([np.var(ar) * len(ar)/N for ar in splar]))
    return SIC


def kv_single_fast(data, threshold, maxiter):
    ''' Search steps with the Kalafut Vischer Algorithm based on the Schwartz Information
    Criterion. steps below threshold will be ignored,  the maximum number of detected steps is
    maxiter. '''
    # typecast to save on memory usage (won't need the precision)
    data = data.astype('float32')
    sz = len(data)
    jpos = np.array([]).astype(int)
    posbysic = []
    sicmin = sic_single(data, [])
    # do not consider steps in last and first 3 datapoints
    indices = np.arange(3, sz - 3)
    
    counter = 0
    for k in range(maxiter):
        
        # build variance array for SIC calculation
        variancesances = np.hstack(([np.var(ar) for ar in np.split(data, jpos)], 0))
        varar = np.tile(variancesances, [sz, 1])
        # build length array for SIC calculation
        steploc = np.hstack((0, jpos, sz))
        lar = np.tile(np.hstack((np.diff(steploc), 0)), [sz, 1])
        
        # run over intervals between previously detected steps to place steps inbetween
        for I in range(len(steploc) - 1):
            # format lengths into lar (length array)
            stepl = int(steploc[I + 1] - steploc[I])
            lar[steploc[I]:steploc[I + 1], I] = np.arange(stepl) + 1
            lar[steploc[I]:steploc[I + 1], -1] = stepl - np.arange(stepl) - 1
            
            # tiled array of data
            dar = np.tile(data[steploc[I]:steploc[I + 1]], [stepl, 1])
            # tiled data values below diagonal (before new step)
            dar_subdiag = dar.copy()
            dar_subdiag[np.triu(np.ones([stepl, stepl]), 1).astype(bool)] = np.nan
            varar[steploc[I]:steploc[I + 1], I] = np.nanvar(dar_subdiag, 1)
            # variances of above diagonal into varar[:,  - 1]
            dar_abovediag = dar.copy()
            dar_abovediag[np.tril(np.ones([stepl, stepl]), 0).astype(bool)] = np.nan
            varar[steploc[I]:steploc[I + 1] - 1, -1] = np.nanvar(dar_abovediag[:-1, :], 1)
            
        # calculate SIC array (for entire trace)
        sic = len(steploc)*np.log(sz) + sz*np.log(np.sum(varar*lar/sz, 1))
        # ignore steps not in indices
        sic = sic[indices]
        
        if np.min(sic) < sicmin:
            # step found
            stp = indices[np.argmin(sic)] + 1
            jtest = sorted(np.hstack((jpos, stp)))
            means = [np.mean(ar) for ar in np.split(data, jtest)]
            stpind = np.where(jtest == stp)[0]
            if np.abs(np.diff(means))[stpind] > threshold:
                # step above threshold  =>  add
                jpos = np.array(jtest).astype(int)
                sicmin = np.min(sic)
                posbysic.append(stp)
                counter = 0
            else:
                counter += 1
                # stop if last maxiter/5 steps were below threshold
                if counter > maxiter/5:
                    break
            # ignore area around step
            jignore = np.argmin(sic) + np.arange( - 2, 3)
            jignore = jignore[(jignore > 0)*(jignore < len(indices) - 1)]
            indices = np.delete(indices, jignore)
        else:
            break
    
    # calculate final variances and means
    variances = np.array([np.var(ar) for ar in np.split(data, jpos)])
    means = np.array([np.mean(ar) for ar in np.split(data, jpos)])
    posbysic = np.array(posbysic)
    
    return jpos, means, variances, posbysic, k


def kv_single(data, threshold, maxiter):
    ''' Search steps with the Kalafut Vischer Algorithm based on the Schwartz Information
    Criterion. steps below threshold will be ignored,  the maximum number of detected steps is
    maxiter. '''
    posbysic = []
    jignore = np.array([]).astype(int)
    sicmin = sic_single(data, [])
    jpos = np.array([]).astype(int)
    
    counter = 0
    for I in range(maxiter):
        # exclude steps around jpos and in jignore
        ind_jpos = np.ravel([jpos - i for i in range(5)])
        ind_jpos = ind_jpos[(ind_jpos >= 0)*(ind_jpos < (len(data) - 5))]
        ind = np.delete(np.arange(2, len(data) - 2), np.hstack((ind_jpos, jignore)))
        # SIC calculation for all step locations
        SIC = np.array([sic_single(data, sorted(np.hstack((jpos, j2)))) for j2 in ind])
        if np.min(SIC) < sicmin:
            #Step improves SIC
            jtest = sorted(np.hstack((jpos, ind[np.argmin(SIC)])))
            means = [np.mean(ar) for ar in np.split(data, jtest)]
            if np.sum(np.abs(np.diff(means)) < threshold):
                # Step ads a mean step below threshold  =>  ignore
                jignore = np.hstack((jignore, ind[np.argmin(SIC)] - np.arange(5)))
                jignore = jignore[(jignore >= 0) & (jignore < (len(data) - 5))]
                counter += 1
                # stop if last maxiter/5 steps were below threshold
                if counter > maxiter/5:
                    break                
            else:
                # Step is valid,  add it to list
                jpos = np.array(jtest)
                sicmin = np.min(SIC)
                posbysic.append(ind[np.argmin(SIC)])
                counter = 0
        else:
            break
        
    means = np.array([np.mean(ar) for ar in np.split(data, jpos)])
    variances = np.array([np.var(ar) for ar in np.split(data, jpos)])
    posbysic = np.array(posbysic)
    
    return jpos, means, variances, posbysic, I 