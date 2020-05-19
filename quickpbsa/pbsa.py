#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import numpy as np
import pandas as pd
import os
import logging
from glob import glob

from .steps_preliminary import kv_file
from .steps_refinement import improve_steps_file

def pbsa_file(infile, threshold, maxiter, outfolder=None,
              preliminary_optional={}, filter_optional={}, refinement_optional={}):
    ''' Run the full photobleaching step analysis on a single file.
    
    Parameters
    ----------
        infile : str
            Input file should be a .csv file with one photobleaching trace
            per line.
        threshold : float
            Kalafut Vischer threshold for preliminary step detection. Try to
            set this at approximately half the height of a typical step.
        maxiter : int
            Maximum iterations in the preliminary step detection. This should
            be significantly higher than the expected number of fluorophores.
        outfolder : str, optional
            output folder for the result files. If not provided output will be saved
            to the directory of the input file.
        preliminary_optional : dict, optional
            Optional parameters for preliminary step detection. For a list of
            possible parameters see documentation.
        filter_optional : dict, optional
            Optional parameters for the exclusion of traces from the analysis. For a list of
            possible parameters see documentation.
        refinement_optional : dict, optional
            Optional parameters for the improvement of the preliminary result.  For a list of
            possible parameters see documentation.
            
    Returns
    _______
        result : DataFrame
            The result DataFrame, with the full analysis result designated with the
            type 'fluors_full'
    '''
    # default Kalafut Vischer Parameters
    kv_param = {'norm': 1,
                'max_memory': 2.0,
                'crop': True,
                'bg_frames': 500}
    kv_param.update(preliminary_optional)
    filter_param = {'subtracted': True,
                    'percentile_step': 90,
                    'length_laststep': 20}
    filter_param.update(filter_optional)
    # default step2 parameters
    step2_param = {'multstep_fraction': 0.5,
                   'nonegatives': False,
                   'mult_threshold': 1.0,
                   'combcutoff': 1000000,
                   'splitcomb': 30000,
                   'maxmult': 5,
                   'maxadded': 5,
                   'lambda': 0.1,
                   'gamma0': 0.5,}
    step2_param.update(refinement_optional)
    # outfolder
    folder, filename = os.path.split(infile)
    if outfolder is None:
        outfolder = os.path.normpath(folder) + '/'
    else:
        outfolder = os.path.normpath(outfolder) + '/'
    logfile = outfolder + filename[:-4] + '.log'
    #remove old logging handlers
    logger = logging.getLogger()
    logger.handlers = []
    logging.basicConfig(filename = logfile,  filemode = 'w',  level = logging.INFO,
                        format = '%(asctime)s %(levelname)s %(message)s')
    # preliminary step detection
    kvout, kvfile, jsonfile = kv_file(infile, threshold, maxiter, outfolder,
                                      kv_param['norm'],
                                      kv_param['max_memory'],
                                      kv_param['crop'],
                                      kv_param['bgframes'])
    # step refinement
    result_out = improve_steps_file(kvfile, jsonfile,
                                    filter_param['subtracted'],
                                    filter_param['percentile_step'],
                                    filter_param['length_laststep'],
                                    step2_param['maxmult'],
                                    step2_param['lambda'],
                                    step2_param['gamma0'],
                                    step2_param['multstep_fraction'],
                                    step2_param['nonegatives'],
                                    step2_param['mult_threshold'],
                                    step2_param['maxadded'],
                                    step2_param['splitcomb'],
                                    step2_param['combcutoff'],
                                    )
    return result_out

def run_pbsa(folder, threshold, maxiter, file_id='_difference.csv',  outfolder=None,
             preliminary_optional={}, filter_optional={}, refinement_optional={}):
    ''' Run the full photobleaching step analysis on a all files in folder
    
    Parameters
    ----------
        folder : str
            Target directory
        threshold : float
            Kalafut Vischer threshold for preliminary step detection. Try to
            set this at approximately half the height of a typical step.
        maxiter : int
            Maximum iterations in the preliminary step detection. This should
            be significantly higher than the expected number of fluorophores.
        file_id : str, optional
            Default '_difference.csv'. Should be a identifier of the files to process.
        kv_optional : dict, optional
            Optional parameters for preliminary step detection. For a list of
            possible parameters see documentation.
        step2_optional : dict, optional
            Optional parameters for the improvement of the preliminary steps. Full
            list of possible parameters can be found in the documentation.
        outfolder : str, optional
            output folder for the result files. If not provided output will be saved
            to input folder.
    '''
    folder = os.path.normpath(folder)
    files = sorted(glob(folder + '/*' + file_id))
    for file in files:
        pbsa_file(file, threshold, maxiter, outfolder,
                  preliminary_optional, filter_optional, refinement_optional)
    return

def summarize_results(folder, outfile=None, avg_over=5, file_id='_result.csv'):
    ''' Summarize result files in a single DataFrame / .csv file
    
    Parameters
    ----------
        folder : str
            Target directory, subfolders will be considered
        outfile : str, optional
            Output file, csv summary will be written into this file. Will overwrite!
        avg_over : int, optional
            First avg_over frames will be averaged to provide initial Intensity
            for each trace. Default is 5.
        file_id : str, optional
            identifier for result files. Default is _result.csv.
            
    Returns
    -------
        df_complete : DataFrame
            Results summary DataFrame. The fluorophore counts, e.g. in 'fluors_full',
            are taken from the first timepoint in the result traces.
    '''
    # list of subfolders
    subfolders = [x[0] for x in os.walk(folder)]
    # output dataframe
    df_complete = pd.DataFrame()
    # run through folders
    for subfolder in subfolders:
        # result files
        files = sorted(glob(subfolder + '/*' + file_id))
        for file in files:
            # Load result dataframe
            data = pd.read_csv(file, header = 1)
            # crop of trace data
            df_file = data.loc[data['type'] == 'fluors_full', :'0']
            df_file = df_file.drop(['type', '0'], axis=1)
            # first avg_over frames
            first_frames = np.array(data.loc[data['type']=='trace', '0':str(avg_over-1)])
            # Intensity from first avg_over frames
            df_file['Intensity'] = np.mean(first_frames, 1)
            # fluors at frame 0
            for typ in ['fluors_intensity', 'fluors_kv', 'fluors_avgintensity', 'fluors_full']:
                df_file[typ] = np.array(data.loc[data['type'] == typ, '0']).astype(int)
            # folder and file into dataframe
            df_file.insert(0, 'folder', subfolder[len(folder):])
            df_file.insert(1, 'file', file.split('/')[-1].split('.')[0])
            # append to output dataframe
            df_complete = df_complete.append(df_file, ignore_index = True)
    if outfile is not None:
        fid = open(outfile, 'w')
        fid.write('# PBSA result summary ' + folder + '\n')
        df_complete.to_csv(fid)
        fid.close()
    return df_complete
    
    

