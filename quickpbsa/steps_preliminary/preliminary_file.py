#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import numpy as np
import os
import pandas as pd
import logging
import time

# relative imports
from .kv_lowlevel import kv_single_fast
from .kv_lowlevel import kv_single
from .helpers import read_tracedf
from .helpers import crop_traces
from .helpers import kv_from_json, kv_to_json
from ..helpers import export_csv

def kv_file(infile, threshold, maxiter, outfolder=None, norm=1, max_memory=2.0, crop=True, bgframes = 500):
    # file names and outfolder
    folder, filename = os.path.split(infile)
    filename, ext = os.path.splitext(filename)
    if outfolder is None:
        outfolder = os.path.normpath(folder) + '/'
    else:
        outfolder = os.path.normpath(outfolder) + '/'
    os.makedirs(outfolder, exist_ok=True)
    # outfiles
    jsonfile = outfolder + filename + '_kvout.json'
    outfile = outfolder + filename + '_kvout.csv'
    # logging
    logfile = outfolder + filename + '_kv.log'
    logger = logging.getLogger('pbsalog')
    if not logger.isEnabledFor(20):
        # set up logger if it is not enabled for INFO level
        logging.basicConfig(filename = logfile,  filemode = 'w',  level = logging.INFO,
                            format = '%(asctime)s %(levelname)s %(message)s')
        logger = logging.getLogger('pbsalog')
    logger.info('Starting KV on ' + str(infile))
    
    # read in Traces
    Traces, basedf, comment, parameters, N_traces, N_frames = read_tracedf(infile)
    # add KV parameter to parameter dictionary
    parameters.update({'KV_threshold': threshold,
                       'maxiter': maxiter,
                       'norm': norm,
                       'crop': crop,
                       'bgframes': bgframes})
    # Norm and flip traces
    Traces = np.fliplr(Traces)/norm
    
    # crop traces if specified (does not actually cut the data, only ignores part of it for analysis purposes)
    if crop:
        crop_index = crop_traces(Traces, threshold/2, bgframes)

    # Prepare output dataframe
    outputs = ['trace', 'kv_mean', 'kv_sdev', 'fluors_intensity', 'fluors_kv']
    if np.size(basedf) == 0:
        result_out = pd.DataFrame()
    else:
        result_out = basedf.iloc[np.repeat(np.arange(N_traces), len(outputs))]
        result_out = result_out.reset_index(drop=True)
    result_out['crop_index'] = 0
    result_out['kv_time [s]'] = 0
    result_out['kv_iter'] = 0
    result_out['laststep'] = 0
    result_out['sdev_laststep'] = 0
    result_out['bg'] = 0
    result_out['sdev_bg'] = 0
    result_out['flag'] = 0
    result_out['type'] = outputs*N_traces
    
    # result array
    result_array = np.zeros([N_traces*len(outputs), N_frames])
    result_array[np.arange(0, N_traces*len(outputs), len(outputs)), :] = Traces
    
    # import json
    if os.path.isfile(jsonfile):
        starttrace = kv_from_json(jsonfile, parameters, result_out, result_array)
    else:
        starttrace = 0
    
    tracetime = []        
    # Main Loop over traces
    for K in range(starttrace, N_traces):
        tic = time.time()
        # cropped trace
        trace = Traces[K, crop_index[K]:]
        if N_frames**2*4/1e9 < max_memory:
            # faster but memory consuming for long traces ( > 10000)
            steppos, means, variances, posbysic, numiter = kv_single_fast(trace, threshold, maxiter)
        else:
            # slower but less memory consuming
            steppos, means, vari, posbysic, numiter = kv_single(trace, threshold, maxiter)
        tracetime.append(time.time() - tic)
        steppos = steppos + crop_index[K]
        posbysic = posbysic + crop_index[K]
        if len(steppos) > 0:
            # calculate fluors by intensity
            means_sub = means - means[0]
            fluors_intensity = np.round(means_sub/means_sub[1]).astype(int)
            # calculate fluors by steps
            stepsigns = np.sign(np.diff(means)) #sign of steps
            fluors_kv = np.cumsum(np.hstack((0, stepsigns)))
            # write into result array
            diffs = np.diff(np.hstack([0, steppos, N_frames]))
            result_array[5*K + 1, :] = np.repeat(means, diffs)
            result_array[5*K + 2, :] = np.repeat(variances, diffs)
            result_array[5*K + 3, :] = np.repeat(fluors_intensity, diffs)
            result_array[5*K + 4, :] = np.repeat(fluors_kv, diffs)
            result_out.loc[5*K:5*(K + 1), 'KV_cutind'] = N_frames - crop_index[K]
            result_out.loc[5*K:5*(K + 1), 'KV_time [s]'] = tracetime[-1]
            result_out.loc[5*K:5*(K + 1), 'KV_iter'] = numiter
            result_out.loc[5*K:5*(K + 1), 'laststep'] = means[1]
            result_out.loc[5*K:5*(K + 1), 'sdev_laststep'] = np.sqrt(variances[1])
            result_out.loc[5*K:5*(K + 1), 'bg'] = means[0]
            result_out.loc[5*K:5*(K + 1), 'sdev_bg'] = np.sqrt(variances[0])
            result_out.loc[5*K:5*(K + 1), 'flag'] = 1
        else:
            # no steps in this trace
            result_out.loc[5*K:5*(K + 1), 'KV_cutind'] = N_frames - crop_index[K]
            result_out.loc[5*K:5*(K + 1), 'KV_time [s]'] = tracetime[-1]
            result_out.loc[5*K:5*(K + 1), 'KV_iter'] = numiter
            result_out.loc[5*K:5*(K + 1), 'bg'] = means[0]
            result_out.loc[5*K:5*(K + 1), 'sdev_bg'] = np.sqrt(variances[0])
            result_out.loc[5*K:5*(K + 1), 'flag'] = -1
        
        # update json file
        kv_to_json(jsonfile, parameters, steppos, means, variances, posbysic,
                   fluors_intensity, fluors_kv, crop_index[K], tracetime[-1], numiter)
        # log
        logmessage = 'trace {:d}, {:d} iterations, runtime {:.2f}s'.format(K, numiter, tracetime[-1])
        logger.info('Finished KV ' + logmessage)
    
    #flip result back
    result_array = np.fliplr(result_array)
    # write result
    result_out = export_csv(result_out, result_array, outfile, parameters, comment)
    # log
    logmessage = '{:d} Traces, avg runtime {:.2f}s'.format(len(tracetime), np.mean(tracetime))
    logger.info('Finished KV for ' + logmessage)
    
    return result_out, outfile, jsonfile



