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
    comment = fid.readline()[:-1]
    if comment[0] == '#':
        # if there is a comment line try to read parameters from it
        try:
            parameters = eval(comment[comment.find('{'):])
            comment = comment[:comment.find('{')]
        except:
            parameters = {}
        fid.close()
        tracedf = pd.read_csv(filename, header=1)
    else:
        # if there is not comment line only read in data
        comment = ''
        parameters = {}
        tracedf = pd.read_csv(filename)
    Traces = np.array(tracedf.loc[:, '0':])
    N_traces, N_frames = np.shape(Traces)
    basedf = tracedf.loc[:, :'0']
    basedf = basedf.drop('0', axis=1)
    return Traces, basedf, comment, parameters, N_traces, N_frames


def crop_traces(Traces, threshold, bgframes):
    # check where differences go over threshold
    indar = (np.diff(Traces) > threshold) * np.arange(np.size(Traces, 1) - 1)
    indar[indar == 0] = np.size(Traces, 1)
    # index where to crop
    crop_index = np.min(indar, 1)
    crop_index[crop_index == np.size(Traces, 1)] = 0
    # go back bg frames
    crop_index = crop_index - bgframes
    crop_index[crop_index < 0] = 0
    return crop_index


def kv_to_json(jsonfile, parameters, steppos, means, variances, posbysic,
               fluors_intensity, fluors_kv, crop_index, kv_time, kv_iter):
    # convert to list to make arrays serializable
    jsondata_out = {'steppos': [steppos.tolist()],
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


def kv_from_json(jsonfile, parameters, resultdf, result_array):
    fp = open(jsonfile)
    jsondata = json.load(fp)
    fp.close()
    if jsondata['parameters'] == list(parameters.values()):
        starttrace = len(jsondata['steppos'])
        for K in range(starttrace):
            if len(jsondata['steppos'][K]) > 0:
                # differences, i.e. frames between steps
                diffs = np.diff([0] + jsondata['steppos'][K] + [np.size(result_array, 1)])
                # write results into array
                result_array[5*K + 1, :] = np.repeat(jsondata['means'][K], diffs)
                result_array[5*K + 2, :] = np.repeat(jsondata['variances'][K], diffs)
                result_array[5*K + 3, :] = np.repeat(jsondata['fluors_intensity'][K], diffs)
                result_array[5*K + 4, :] = np.repeat(jsondata['fluors_kv'][K], diffs)
                # write metadata into results dataframe
                resultdf.loc[5*K:5*(K + 1), 'crop_index'] = jsondata['crop_index'][K]
                resultdf.loc[5*K:5*(K + 1), 'kv_time [s]'] = jsondata['kv_time'][K]
                resultdf.loc[5*K:5*(K + 1), 'kv_iter'] = jsondata['kv_iter'][K]
                resultdf.loc[5*K:5*(K + 1), 'laststep'] = jsondata['means'][K][1]
                resultdf.loc[5*K:5*(K + 1), 'sdev_laststep'] = np.sqrt(jsondata['variances'][K][1])
                resultdf.loc[5*K:5*(K + 1), 'bg'] = jsondata['means'][K][0]
                resultdf.loc[5*K:5*(K + 1), 'sdev_bg'] = np.sqrt(jsondata['variances'][K][0])
                resultdf.loc[5*K:5*(K + 1), 'flag'] = 1
            else:
                # no steps
                result_array[K*5 + 1, :] = jsondata['means'][K][0]
                result_array[K*5 + 2, :] = jsondata['variances'][K][0]
                resultdf.loc[5*K:5*(K + 1), 'crop_index'] = jsondata['crop_index'][K]
                resultdf.loc[5*K:5*(K + 1), 'kv_time [s]'] = jsondata['kv_time'][K]
                resultdf.loc[5*K:5*(K + 1), 'kv_iter'] = jsondata['kv_iter'][K]
                resultdf.loc[5*K:5*(K + 1), 'bg'] = jsondata['means'][K][0]
                resultdf.loc[5*K:5*(K + 1), 'sdev_bg'] = np.sqrt(jsondata['variances'][K][0])
                resultdf.loc[5*K:5*(K + 1), 'flag'] =  -1
    else:
        starttrace = 0
    return starttrace



