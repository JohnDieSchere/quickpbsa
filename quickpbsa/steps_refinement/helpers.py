#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:59:08 2020

@author: johan
"""

import numpy as np
import pandas as pd
import os
import json


def read_kvresult(filename):
    fid = open(filename)
    # read in comment and parameters from it
    comment = fid.readline()[:-1]
    parameters = eval(comment[comment.find('{'):])
    comment = comment[:comment.find('{')]
    # read in KVresult and separate traces
    kvresult = pd.read_csv(filename, header=1)
    Traces = np.array(kvresult.loc[kvresult['type']=='trace', '0':])
    N_traces, N_frames = np.shape(Traces)
    # separate base dataframe
    basedf = kvresult.loc[kvresult['type']=='trace', :'0']
    basedf = basedf.drop(['0', 'type'], axis=1)
    return kvresult, Traces, basedf, comment, parameters, N_traces, N_frames


def step2_from_json(jsonfile, parameters, N_traces):
    fp = open(jsonfile)
    jsondata = json.load(fp)
    fp.close()
    # traces yet to be analysed
    trace_indices = np.arange(N_traces)
    if jsondata['parameters'] == list(parameters.values()):
        analysed_traces = jsondata['trace_index']
        # delete trace indices already in json file
        trace_indices = np.delete(trace_indices, analysed_traces)
    return trace_indices


def step2_to_json(jsonfile, parameters, trace_index, steppos, steps, step2_time, sic_final, flag):
    # convert to list to make arrays serializable
    jsondata_out = {'trace_index': [int(trace_index)],
                    'steppos': [steppos.tolist()],
                    'steps': [steps.tolist()],
                    'step2_time': [step2_time],
                    'sic_final': [sic_final],
                    'flag': [int(flag)]}
    jsondata_out['parameters'] = list(parameters.values())
    if os.path.isfile(jsonfile):
        fp = open(jsonfile)
        jsondata_file = json.load(fp)
        fp.close()
        if jsondata_file['parameters'] == list(parameters.values()):
            jsondata_out = jsondata_file
            # convert to list to make arrays serializable
            jsondata_out['trace_index'].append(int(trace_index))
            jsondata_out['steppos'].append(steppos.tolist())
            jsondata_out['steps'].append(steps.tolist())
            jsondata_out['step2_time'].append(step2_time)
            jsondata_out['sic_final'].append(sic_final)
            jsondata_out['flag'].append(int(flag))
    fp = open(jsonfile, 'w')
    json.dump(jsondata_out, fp)
    fp.close()
    return


def step2result_from_json(jsonfile, kvout, flags, avg_laststep):
    # Read in Traces
    kvresult, Traces, basedf, comment, parameters, N_traces, N_frames = read_kvresult(kvout)
    types_kv = int(np.size(kvresult, 0)/N_traces)
    # Prepare output dataframe
    outputs = ['trace', 'kv_mean', 'kv_variance', 'fluors_intensity',
               'fluors_kv', 'fluors_avgintensity', 'fluors_full']
    No = len(outputs)
    result_out = basedf.iloc[np.repeat(np.arange(N_traces), No)]
    result_out = result_out.reset_index(drop=True)
    result_out['step2_time [s]'] = 0
    result_out['sic_final'] = 0
    result_out['type'] = outputs*N_traces
    result_out['flag'] = np.repeat(flags, No)
    # result array
    result_array = np.zeros([N_traces*No, N_frames])
    # kvresult into result array
    for i in range(types_kv):
        indices = np.arange(i, N_traces*No, No)
        result_array[indices, :] = np.fliplr(kvresult.loc[kvresult['type'] == outputs[i], '0':])
    # fluorophores by average intensity
    fluors_avgintensity = np.round(Traces/avg_laststep).astype(int)
    indices = np.arange(types_kv, N_traces*No, No)
    result_array[indices, :] = fluors_avgintensity
    # valid traces with only one fluorophore
    indices = (np.where(flags==1)[0] + No - 1) * No
    result_array[indices, :] = result_array[indices-No+types_kv, :]
    # read in json
    fp = open(jsonfile)
    jsondata = json.load(fp)
    fp.close()
    # write to result
    num_entries = len(jsondata['trace_index'])
    for I in range(num_entries):
        K = jsondata['trace_index'][I]
        result_out.loc[No*K:No*(K + 1), 'step2_time'] = jsondata['step2_time'][I]
        result_out.loc[No*K:No*(K + 1), 'sic_final'] = jsondata['sic_final'][I]
        result_out.loc[No*K:No*(K + 1), 'flag'] = jsondata['flag'][I]
        steppos = jsondata['steppos'][I]
        steps = jsondata['steps'][I]
        # diffs (frames between steps)
        diffs = np.diff(np.hstack([0, steppos, N_frames]))
        # fluors from steps
        fluors = np.cumsum(np.hstack([0, steps]))
        # write result into array
        result_array[No*K + No -1, :] = np.repeat(fluors, diffs)
        
    result_array = np.fliplr(result_array)
    
    return result_out, result_array
        