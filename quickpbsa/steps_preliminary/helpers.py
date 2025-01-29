#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import pandas as pd
import numpy as np
import json
import os

def read_tracedf(filename):
    fid = open(filename)
    # read in comment
    firstline = fid.readline()[:-1]
    if firstline[0] == '#':
        # first line is comment
        try:
            parameters = eval(firstline[firstline.find('{'):])
            comment = firstline[:firstline.find('{')]
            skiprows = 1
        except:
            comment = ''
            parameters = {}
            skiprows = 0
    fid.close()
    try:
        tracedf = pd.read_csv(filename, header=0, skiprows=skiprows)
        basedf = tracedf.loc[:,:'0']
        basedf = basedf.iloc[:,:-1]
    except:
        tracedf = pd.read_csv(filename, header=None)
        basedf = tracedf.loc[:,:'0']
        basedf = basedf.iloc[:,:-1]
    Traces = np.array(tracedf.loc[:, '0':])
    N_traces, N_frames = np.shape(Traces)
    if basedf.shape[1] == 0:
        basedf['id'] = np.arange(N_traces)
    return Traces, basedf, comment, parameters, N_traces, N_frames


def crop_traces(Traces, threshold, bgframes, crop):
    # clip of nans
    indar = (np.invert(np.isnan(Traces)) * np.arange(np.size(Traces, 1)))
    indar[indar == 0] = np.size(Traces, 1)
    crop_index_nan = np.min(indar, 1)
    crop_index_nan[crop_index_nan == np.size(Traces, 1)] = 0
    if crop:
        # check where differences go over threshold
        indar = (np.diff(Traces) > threshold) * np.arange(np.size(Traces, 1) - 1)
        indar[indar == 0] = np.size(Traces, 1)
        # index where to crop
        crop_index = np.min(indar, 1)
        crop_index[crop_index == np.size(Traces, 1)] = 0
        # go back bg frames
        crop_index = crop_index - bgframes
        crop_index[crop_index < 0] = 0
        crop_index_nan = np.maximum(crop_index, crop_index_nan)
    return crop_index_nan


def kv_to_json(jsonfile, parameters, trace_index, steppos, means, variances, posbysic,
               fluors_intensity, fluors_kv, crop_index, kv_time, kv_iter):
    # convert to list to make arrays serializable
    jsondata_out = {'trace_index': [int(trace_index)],
                    'steppos': [steppos.tolist()],
                    'means': [means.tolist()],
                    'variances': [variances.tolist()], 
                    'posbysic': [posbysic.tolist()],
                    'fluors_intensity': [fluors_intensity.tolist()],
                    'fluors_kv': [fluors_kv.tolist()],
                    'crop_index': [int(crop_index)],
                    'kv_time': [kv_time],
                    'kv_iter': [kv_iter]}
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
            jsondata_out['means'].append(means.tolist())
            jsondata_out['variances'].append(variances.tolist())
            jsondata_out['posbysic'].append(posbysic.tolist())
            jsondata_out['fluors_intensity'].append(fluors_intensity.tolist())
            jsondata_out['fluors_kv'].append(fluors_kv.tolist())
            jsondata_out['crop_index'].append(int(crop_index))
            jsondata_out['kv_time'].append(kv_time)
            jsondata_out['kv_iter'].append(kv_iter)
    fp = open(jsonfile, 'w')
    json.dump(jsondata_out, fp)
    fp.close()
    return


def kv_from_json(jsonfile, parameters, N_traces):
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


def kvresult_from_json(jsonfile, Traces, basedf):   
    N_traces, N_frames = np.shape(Traces)
    # Prepare output dataframe
    outputs = ['trace', 'kv_mean', 'kv_variance', 'fluors_intensity', 'fluors_kv']
    No = len(outputs)
    if np.size(basedf) == 0:
        result_out = pd.DataFrame()
    else:
        result_out = basedf.iloc[np.repeat(np.arange(N_traces), No)]
        result_out = result_out.reset_index(drop=True)
    result_out['crop_index'] = 0
    result_out['kv_time [s]'] = 0.0
    result_out['kv_iter'] = 0
    result_out['laststep'] = 0.0
    result_out['sdev_laststep'] = 0.0
    result_out['bg'] = 0.0
    result_out['sdev_bg'] = 0.0
    result_out['flag'] = 0
    result_out['type'] = outputs*N_traces
    
    # result array
    result_array = np.zeros([N_traces*No, N_frames])
    result_array[np.arange(0, N_traces*No, No), :] = Traces
    
    # read from json file
    fp = open(jsonfile)
    jsondata = json.load(fp)
    fp.close()
    
    num_entries = len(jsondata['trace_index'])
    for I in range(num_entries):
        K = jsondata['trace_index'][I]
        if len(jsondata['steppos'][I]) > 0:
            # differences, i.e. frames between steps
            diffs = np.diff([0] + jsondata['steppos'][I] + [np.size(result_array, 1)])
            # write results into array
            result_array[No*K + 1, :] = np.repeat(jsondata['means'][I], diffs)
            result_array[No*K + 2, :] = np.repeat(jsondata['variances'][I], diffs)
            result_array[No*K + 3, :] = np.repeat(jsondata['fluors_intensity'][I], diffs)
            result_array[No*K + 4, :] = np.repeat(jsondata['fluors_kv'][I], diffs)
            # write metadata into results dataframe
            result_out.loc[No*K:No*(K + 1) - 1, 'crop_index'] = jsondata['crop_index'][I]
            result_out.loc[No*K:No*(K + 1) - 1, 'kv_time [s]'] = jsondata['kv_time'][I]
            result_out.loc[No*K:No*(K + 1) - 1, 'kv_iter'] = jsondata['kv_iter'][I]
            result_out.loc[No*K:No*(K + 1) - 1, 'laststep'] = jsondata['means'][I][1]
            result_out.loc[No*K:No*(K + 1) - 1, 'sdev_laststep'] = np.sqrt(jsondata['variances'][I][1])
            result_out.loc[No*K:No*(K + 1) - 1, 'bg'] = jsondata['means'][I][0]
            result_out.loc[No*K:No*(K + 1) - 1, 'sdev_bg'] = np.sqrt(jsondata['variances'][I][0])
            result_out.loc[No*K:No*(K + 1) - 1, 'flag'] = 1
        else:
            # no steps
            result_array[K*No + 1, :] = jsondata['means'][I][0]
            result_array[K*No + 2, :] = jsondata['variances'][I][0]
            result_out.loc[No*K:No*(K + 1) - 1, 'crop_index'] = jsondata['crop_index'][I]
            result_out.loc[No*K:No*(K + 1) - 1, 'kv_time [s]'] = jsondata['kv_time'][I]
            result_out.loc[No*K:No*(K + 1) - 1, 'kv_iter'] = jsondata['kv_iter'][I]
            result_out.loc[No*K:No*(K + 1) - 1, 'bg'] = jsondata['means'][I][0]
            result_out.loc[No*K:No*(K + 1) - 1, 'sdev_bg'] = np.sqrt(jsondata['variances'][I][0])
            result_out.loc[No*K:No*(K + 1) - 1, 'flag'] =  -1
            
    result_array = np.fliplr(result_array)
        
    return result_out, result_array



