#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:59:08 2020

@author: johan
"""

import numpy as np
import json
import logging
import os
import time
import multiprocessing as mp

from .improve_steps import improve_steps_single
from ..helpers import export_csv
from .helpers import read_kvresult
from .helpers import step2_from_json
from .helpers import step2_to_json
from .helpers import step2result_from_json



def generate_flags(kvjson, KVthreshold, subtracted=True, percentile_step=90, length_laststep=20):
    '''
    Flag out traces based on the last steps from the preliminary step detection
    '''
        
    # number of steps
    nsteps = np.array([len(steppos) for steppos in kvjson['steppos']])
    # initialize flags
    flags = np.zeros(len(kvjson['steppos']), dtype = int)
    flags[nsteps == 0] =  -1 # zero steps
    flagsel1 = np.where(nsteps > 0)[0]
    flags[nsteps == 1] = 1  # 1 step,  no need for further analysis
    flagsel2 = np.where(nsteps > 1)[0]
    
    # background
    bg = np.array([means[0] for means in kvjson['means']])
    if subtracted:
        lower_bg =  -1*KVthreshold
        upper_bg = KVthreshold
    else:
        lower_bg = 0
        upper_bg = np.min(bg) + KVthreshold
    sel = np.invert((bg > lower_bg) & (bg < upper_bg))        
    flags[sel] =  -2
        
    # last step
    laststep = np.array([means[1] - means[0] for means in kvjson['means'] if len(means) > 1])
    if isinstance(percentile_step, int) or isinstance(percentile_step, float):
        # only one border provided,  use threshold as lower border
        if subtracted:
            lower = KVthreshold
        else:
            lower = np.min(bg) + KVthreshold
        upper = np.percentile(laststep, percentile_step)
    else:
        lower = np.percentile(laststep, percentile_step[0])
        upper = np.percentile(laststep, percentile_step[1])
    sel = np.invert((laststep >= lower) & (laststep <= upper))
    # to not override old flags
    ind_new = np.intersect1d(flagsel1[sel], np.where(flags >= 0)[0])
    flags[ind_new] =  -3    
    
    # trace means goes into negative values
    sel = np.array([np.any(np.array(means) < lower_bg) for means in kvjson['means']])
    ind_new = np.intersect1d(np.where(sel)[0], np.where(flags == 0)[0])
    flags[ind_new] =  -4
    
    # KV found negative fluorophore numbers:
    sel = [np.any(np.cumsum(np.sign(np.diff(means))) < 0) for means in kvjson['means'] if len(means) > 2]
    ind_new = np.intersect1d(flagsel2[sel], np.where(flags == 0)[0])
    flags[ind_new] =  -5   
    
    # length of last step
    lastint = np.array([steppos[1] - steppos[0] for steppos in kvjson['steppos'] if len(steppos) > 1])
    sel = (lastint < length_laststep)
    ind_new = np.intersect1d(flagsel2[sel], np.where(flags == 0)[0])
    flags[ind_new] =  -6
    
    # calculate average last step
    avg_laststep = np.mean(laststep[(flags[flagsel1] == 0) | (flags[flagsel1] == 1)])
    
    # reorder according to trace index
    index = np.argsort(kvjson['trace_index'])
    flags = flags[index]
    
    return flags, avg_laststep



def worker(trace_index, trace, refinement_args, Q):
    tic = time.time()
    result = [trace_index]
    res = improve_steps_single(trace, *refinement_args)
    result.append(res)
    if res[0]:
        flag=1
    else:
        flag=-7
    result.append(flag)
    result.append(time.time() - tic)
    Q.put(result)
    return result


def listener(Q, jsonfile, parameters, logger):
    while True:
        result = Q.get()
        if result == 'kill':
            break
        trace_index = result[0]
        success, sic_final, steppos, steps = result[1]
        flag = result[2]
        step2_time = result[3]
        # write to json file
        step2_to_json(jsonfile, parameters, trace_index, steppos, steps, step2_time, sic_final, flag)
        # write to log
        if flag == 1:
            logmessage = 'successful for {:d}, runtime {:.2f}s'.format(trace_index, step2_time)
        else:
            logmessage = 'unsuccessful for {:d}, runtime {:.2f}s'.format(trace_index, step2_time)            
        logger.info('Step refinement ' + logmessage)
    return
    
    



def improve_steps_file(kvout, kvjson, subtracted=True, percentile_step=90, length_laststep=20,
                       maxmult=5, lamb=0.1, gamma0=0.5, multstep_fraction=0.5,
                       nonegatives=False, mult_threshold=1, maxadded=5, splitcomb=30000,
                       combcutoff=1000000, num_cores=2):
    '''
    Run the step refinement on the output files from the preliminary step detection
    '''
    
    # file name (cleave off _kvout.csv)
    filename = '_'.join(kvout.split('_')[:-1])
    jsonfile = filename + '_result.json'
    outfile = filename + '_result.csv'
    logfile = filename + '_step2.log'
    
    # Start logging
    logger = logging.getLogger('pbsalog')
    if not logger.isEnabledFor(20):
        logging.basicConfig(filename = logfile,  filemode = 'w',  level = logging.INFO, format = '%(asctime)s %(levelname)s %(message)s')
        logger = logging.getLogger('pbsalog')
    logger.info('Starting Step 2 on ' + str(kvout))

    # Read in Traces
    kvresult, Traces, basedf, comment, parameters, N_traces, N_frames = read_kvresult(kvout)
    # update parameter dictionary
    parameters.update({'subtracted': subtracted,
                       'percentile_step': percentile_step,
                       'length_laststep': length_laststep,
                       'maxmult': maxmult,
                       'lamb': lamb,
                       'gamma0': gamma0,
                       'multstep_fraction': multstep_fraction,
                       'nonegatives': nonegatives,
                       'mult_threshold': mult_threshold,
                       'maxadded': maxadded,
                       'splitcomb': splitcomb,
                       'combcutoff': combcutoff})
    # flip traces
    Traces = np.fliplr(Traces)
    
    # read in kvjson
    fp = open(kvjson)
    kvjsondata = json.load(fp)
    fp.close()
    # generate flages
    flags, avg_laststep = generate_flags(kvjsondata, parameters['KV_threshold'],
                                         subtracted, percentile_step, length_laststep)

    # import json
    if os.path.isfile(jsonfile):
        trace_indices = step2_from_json(jsonfile, parameters, N_traces)
    else:
        trace_indices = np.arange(N_traces)
        
    # remove flagged out traces
    flags_sel = flags[trace_indices]
    trace_indices = np.delete(trace_indices, np.where(flags_sel != 0)[0])
    
    # need at least 2 cores for listener and worker
    if num_cores < 2:
        num_cores = 2
    
    # Set up pool of workers
    manager = mp.Manager()
    Q = manager.Queue()    
    pool = mp.Pool(num_cores)
    
    # start listener
    pool.apply_async(listener,
                      (Q, jsonfile, parameters, logger))
    
    jobs = []           
    # Main Loop over traces
    for K in trace_indices:
        # inputs from kvjson
        kvindex = int(np.where(kvjsondata['trace_index']==K)[0])
        steppos_kv = np.array(kvjsondata['steppos'][kvindex])
        means = np.array(kvjsondata['means'][kvindex])
        variances = np.array(kvjsondata['variances'][kvindex])
        posbysic = np.array(kvjsondata['posbysic'][kvindex])
        # trace
        trace = Traces[K, :]
        # refinement args
        refinement_args = (steppos_kv, means, variances, posbysic, maxmult, lamb, gamma0,
                           multstep_fraction, nonegatives, mult_threshold, maxadded,
                           splitcomb, combcutoff)
        job = pool.apply_async(worker,
                               (K, trace, refinement_args, Q))
        jobs.append(job)
        
    # collect results from the workers through the pool result queue
    for job in jobs:
        job.get()

    #now we are done, kill the listener
    Q.put('kill')
    pool.close()
    pool.join()
        
    # write result
    result_out, result_array = step2result_from_json(jsonfile, kvout, flags, avg_laststep)
    result_out = export_csv(result_out, result_array, outfile, parameters, comment)
    # log
    logmessage = '{:d} Traces '.format(N_traces)
    logger.info('Finished step 2 for ' + logmessage)
    
    return result_out
    