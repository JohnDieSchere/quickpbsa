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


def step2_to_json(jsonfile, parameters, steppos, steps, step2_time, sic_final, flag):
    # convert to list to make arrays serializable
    jsondata_out = {'steppos': [steppos.tolist()],
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
            jsondata_out['steppos'].append(steppos.tolist())
            jsondata_out['steps'].append(steps.tolist())
            jsondata_out['step2_time'].append(step2_time)
            jsondata_out['sic_final'].append(sic_final)
            jsondata_out['flag'].append(int(flag))
    fp = open(jsonfile, 'w')
    json.dump(jsondata_out, fp)
    fp.close()
    return


def step2_from_json(jsonfile, parameters, resultdf, result_array):
    fp = open(jsonfile)
    jsondata = json.load(fp)
    fp.close()
    if jsondata['parameters'] == list(parameters.values()):
        starttrace = len(jsondata['flag'])
        for K in range(starttrace):
            if jsondata['flag'][K] == 1 or jsondata['flag'][K] == -7:
                # differences, i.e. frames between steps
                diffs = np.diff([0] + jsondata['steppos'][K] + [np.size(result_array, 1)])
                # fluors from steps
                fluors = np.cumsum([0] + jsondata['steps'][K])
                # write results into array
                result_array[7*K + 6, :] = np.repeat(fluors, diffs)
                # write metadata into results dataframe
                resultdf.loc[7*K:7*(K + 1), 'step2_time [s]'] = jsondata['step2_time'][K]
                resultdf.loc[7*K:7*(K + 1), 'sic_final'] = jsondata['sic_final'][K]
                resultdf.loc[7*K:7*(K + 1), 'flag'] = jsondata['flag'][K]
            else:
                # Trace flagged out
                resultdf.loc[7*K:7*(K + 1), 'flag'] = jsondata['flag'][K]                
    else:
        starttrace = 0
    return starttrace