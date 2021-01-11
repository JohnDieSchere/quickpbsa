#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import numpy as np
import os
import logging
import time
import multiprocessing as mp

# relative imports
from .kv_lowlevel import kv_single_fast
from .kv_lowlevel import kv_single
from .helpers import read_tracedf
from .helpers import crop_traces
from .helpers import kv_from_json, kv_to_json, kvresult_from_json
from ..helpers import export_csv


def worker(trace_index, trace, threshold, maxiter, max_memory, Q):
    tic = time.time()
    result = [trace_index]
    if len(trace)**2*4/1e9 < max_memory:
        # faster but memory consuming for long traces ( > 10000)
        result.append(kv_single_fast(trace, threshold, maxiter))
    else:
        # slower but less memory consuming
        result.append(kv_single(trace, threshold, maxiter))
    result.append(time.time() - tic)
    Q.put(result)
    return result


def listener(Q, crop_index, jsonfile, parameters, logger):
    while True:
        result = Q.get()
        if result == 'kill':
            break
        trace_index = result[0]
        # kv result
        steppos, means, variances, posbysic, kv_iter = result[1]
        # correct for cropping
        steppos = steppos + crop_index[trace_index]
        posbysic = posbysic + crop_index[trace_index]
        if len(steppos) > 0:
            # calculate fluors by intensity
            means_sub = means - means[0]
            fluors_intensity = np.round(means_sub/means_sub[1]).astype(int)
            # calculate fluors by steps
            stepsigns = np.sign(np.diff(means)) #sign of steps
            fluors_kv = np.cumsum(np.hstack((0, stepsigns)))
        else:
            # fluors (zero)
            fluors_intensity = np.array([0])
            fluors_kv = np.array([0])
        kv_time = result[2]
        # write to json file
        kv_to_json(jsonfile, parameters, trace_index, steppos, means, variances, posbysic,
                    fluors_intensity, fluors_kv, 0, kv_time, kv_iter)
        # write to log
        logmessage = 'trace {:d}, {:d} iterations, runtime {:.2f}s'.format(trace_index, kv_iter, kv_time)
        logger.info('Finished KV ' + logmessage)
    return
        
    


def kv_file(infile, threshold, maxiter, outfolder=None, norm=1, max_memory=2.0, crop=True, bgframes = 500, num_cores=2):
    '''
    Run preliminary step detection on a .csv file with rows as photobleaching traces
    '''
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
                       'bg_frames': bgframes})
    # Norm and flip traces
    Traces = np.fliplr(Traces)/norm
    
    # crop traces if specified (does not actually cut the data, only ignores part of it for analysis purposes)
    if crop:
        crop_index = crop_traces(Traces, threshold/2, bgframes)
    
    # import json
    if os.path.isfile(jsonfile):
        trace_indices = kv_from_json(jsonfile, parameters, N_traces)
    else:
        trace_indices = np.arange(N_traces)
    
    # need at least 2 cores for listener and worker
    if num_cores < 2:
        num_cores = 2
    
    # Set up pool of workers
    manager = mp.Manager()
    Q = manager.Queue()    
    pool = mp.Pool(num_cores)
    
    # start listener
    pool.apply_async(listener,
                     (Q, crop_index, jsonfile, parameters, logger))
    
    jobs = []      
    # Start workers on traces
    for K in trace_indices:
        # cropped trace
        trace = Traces[K, crop_index[K]:]
        job = pool.apply_async(worker,
                               (K, trace, threshold, maxiter, max_memory, Q))
        jobs.append(job)
        
    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #now we are done, kill the listener
    Q.put('kill')
    pool.close()
    pool.join()
        
    # write result
    result_out, result_array = kvresult_from_json(jsonfile, Traces, basedf)
    result_out = export_csv(result_out, result_array, outfile, parameters, comment)
    # log
    logmessage = '{:d} Traces'.format(len(trace_indices))
    logger.info('Finished KV for ' + logmessage)
    
    return result_out, outfile, jsonfile



