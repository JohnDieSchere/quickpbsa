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

from .improve_steps import improve_steps_single
from ..helpers import export_csv
from .helpers import read_kvresult
from .helpers import step2_from_json
from .helpers import step2_to_json



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
    sel = np.invert((laststep > lower) & (laststep < upper))
    # to not override old flags
    ind_new = np.intersect1d(flagsel1[sel], np.where(flags == 0)[0])
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
    
    return flags, avg_laststep



def improve_steps_file(kvout, kvjson, subtracted=True, percentile_step=90, length_laststep=20,
                       maxmult=5, lamb=0.1, gamma0=0.5, multstep_fraction=0.5,
                       nonegatives=False, mult_threshold=1, maxadded=5, splitcomb=30000, combcutoff=1000000):
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
    
    # Prepare output dataframe
    outputs = ['trace', 'kv_mean', 'kv_variance', 'fluors_intensity',
               'fluors_kv', 'fluors_avgintensity', 'fluors_full']
    result_out = basedf.iloc[np.repeat(np.arange(N_traces), len(outputs))]
    result_out = result_out.reset_index(drop=True)
    result_out['step2_time [s]'] = 0
    result_out['sic_final'] = 0
    result_out['type'] = outputs*N_traces

    # read in kvjson
    fp = open(kvjson)
    kvjsondata = json.load(fp)
    fp.close()
    # generate flages
    flags, avg_laststep = generate_flags(kvjsondata, parameters['KV_threshold'],
                                         subtracted, percentile_step, length_laststep)
    result_out['flag'] = np.repeat(flags, 7)
    # result array
    result_array = np.zeros([N_traces*len(outputs), N_frames])
    # kvresult into result array
    for i in range(5):
        indices = np.arange(i, N_traces*len(outputs), len(outputs))
        result_array[indices, :] = np.fliplr(kvresult.loc[kvresult['type'] == outputs[i], '0':])
    # fluorophores by average intensity
    fluors_avgintensity = np.round(Traces/avg_laststep).astype(int)
    result_array[np.arange(5, N_traces*len(outputs), len(outputs)), :] = fluors_avgintensity
    
    # import json
    if os.path.isfile(jsonfile):
        starttrace = step2_from_json(jsonfile, parameters, result_out, result_array)
    else:
        starttrace = 0
        
    tracetime = []        
    # Main Loop over traces
    for K in range(starttrace, N_traces):
        # inputs from kvjson
        steppos_kv = np.array(kvjsondata['steppos'][K])
        means = np.array(kvjsondata['means'][K])
        variances = np.array(kvjsondata['variances'][K])
        posbysic = np.array(kvjsondata['posbysic'][K])
        if flags[K] == 0:
            tic = time.time()
            # call to improve steps
            success, sicmin, steppos_out, steps_out = improve_steps_single(Traces[K, :], steppos_kv,
                                                                           means, variances, posbysic,
                                                                           maxmult, lamb, gamma0, multstep_fraction,
                                                                           nonegatives, mult_threshold, maxadded,
                                                                           splitcomb, combcutoff)
            tracetime.append(time.time()-tic)
            if success:
                result_out.loc[7*K:7*(K + 1), 'flag'] = 1
                # update json
                step2_to_json(jsonfile, parameters, steppos_out, steps_out,
                              tracetime[-1], sicmin, 1)
                # log
                logmessage = '{:d} analysed succesfully, runtime {:.2f}s'.format(K, tracetime[-1])
                logger.info('Trace ' + logmessage)
            else:
                result_out.loc[7*K:7*(K + 1), 'flag'] = -7
                step2_to_json(jsonfile, parameters, steppos_out, steps_out,
                              tracetime[-1], sicmin, -7)
                # log
                logmessage = '{:d} not analysed succesfully, runtime {:.2f}s'.format(K, tracetime[-1])
                logger.info('Trace ' + logmessage)
            result_out.loc[7*K:7*(K + 1), 'step2_time'] = tracetime[-1]
            result_out.loc[7*K:7*(K + 1), 'sic_final'] = sicmin
            # diffs (frames between steps)
            diffs = np.diff(np.hstack([0, steppos_out, N_frames]))
            # fluors from steps
            fluors = np.cumsum(np.hstack([0, steps_out]))
            # write result into array
            result_array[7*K + 6, :] = np.repeat(fluors, diffs)
        elif flags[K] == 1:
            # only one step in trace
            result_out.loc[7*K:7*(K + 1), 'flag'] = 1
            # steps
            steps_out = np.array([1])
            # update json
            step2_to_json(jsonfile, parameters, steppos_kv, steps_out, 0, 0, 1)
            # diffs (frames between steps)
            diffs = np.diff(np.hstack([0, steppos_kv, N_frames]))
            # fluors from steps
            fluors = np.cumsum(np.hstack([0, steps_out]))
            # write result into array
            result_array[7*K + 6, :] = np.repeat(fluors, diffs)
            logmessage = '{:d} skipped, only one step'.format(K)
            logger.info('Trace ' + logmessage)
        else:
            # trace flagged out, only update json
            step2_to_json(jsonfile, parameters, np.array([]), np.array([]), 0, 0, flags[K])
            result_out.loc[7*K:7*(K + 1), 'flag'] = flags[K]
            # log
            logmessage = '{:d} flagged out: flag {:d}'.format(K, flags[K])
            logger.info('Trace ' + logmessage)
    
    #flip result back
    result_array = np.fliplr(result_array)
    # write result
    result_out = export_csv(result_out, result_array, outfile, parameters, comment)
    # log
    logmessage = '{:d} Traces, avg runtime {:.2f}s'.format(len(tracetime), np.mean(tracetime))
    logger.info('Finished step 2 for ' + logmessage)
    
    return result_out
    