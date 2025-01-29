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

def pbsa_file(infile, threshold, maxiter, outfolder=None, num_cores=1,
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
            Optional parameters for preliminary step detection. Possible Parameters
            detailed below.
        filter_optional : dict, optional
            Optional parameters for the exclusion of traces from the analysis. Possible parameters
            detailed below.
        refinement_optional : dict, optional
            Optional parameters for the improvement of the preliminary result. Possible parameters
            detailed below.
            
    Returns
    _______
        result : DataFrame
            The result DataFrame, with the full analysis result designated with the
            type 'fluors_full'
            
    Notes
    _____
    
    Possible Parameters in preliminary_optional:
        - 'norm': divide all traces by this before analysis (default 1)
        - 'max_memory': Maximum Memory to use in GB. If this is exceeded for long traces,\
        a slower but less memory concuming implementation is used. (default 4.0)
        - 'crop': Automatically crop traces for the analysis if True, based on the threshold (default True)
        - 'bg_frames': How many timepoints to include after the crop point (default 500)
        
    Possible Parameters in filter_optional:
        - 'subtracted': If `True` it is assumed that traces are background corrected.\
        This sets the bounds on the background intensity to [-threshold, threshold]\
        and the default lower bound on the single fluorophore intensity to threshold.\
        If False the bounds on the background intensity are set based on the minimum\
        background intensity in the dataset min_bg: [min_bg, min_bg + threshold].\
        If False the default lower bound on the single fluorophore intensity is also\
        min_bg + threshold. (Default True)
        - 'percentile_step': Sets the bounds on the single fluorophore intensity.\
        If one value is provided, the upper bound on the single fluorophore intensity\
        is set at this percentile. If two values are provided, \
        lower and upper bounds are set at the percentiles respectively. (default 90)  
        - 'length_laststep': Minimum number of frames between the last two steps. (default 20)
    
    Possible Parameters in refinement_optional:      
        - 'multstep_fraction': Maximum fraction of steps with an occupancy higher than 1 (default 0.5) 
        - 'nonegatives': If True, no negative double steps are considered.\
        This means that arrangements where 2 or more fluorophore turn back\
        on at the same time are not considered. (default False)
        - 'mult_threshold' Only steps where the difference in the mean is above\
        mult_threshold multiplied by the last fluorophore intensity are considered\
        as steps with occupancy higher than 1. (default 1.0) 
        - 'combcutoff': Maximum number of arrangements to test. If this is exceeded,\
        the trace is flagged out with flag `-7`. (default 2E6)
        - 'splitcomb': How many arrangements to test simultaneously (vectorized). (default 30000)  
        - 'maxmult': Maximum considered occupancy, i.e. how many fluorophores can bleach simultaneously. (default 5)  
        - 'maxadded': Maximum number of added single steps if no steps are removed\
        to yield an improved posterior. (default 10) 
        - 'lambda': Hyperparameter lambda. (default 0.1)  
        - 'gamma0': Hyperparameter gamma_0. (default 0.5)
    
    '''
    # default Kalafut Vischer Parameters
    kv_param = {'norm': 1,
                'max_memory': 4.0,
                'crop': False,
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
                   'combcutoff': 2000000,
                   'splitcomb': 30000,
                   'maxmult': 5,
                   'maxadded': 10,
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
                                      kv_param['bg_frames'],
                                      num_cores)
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
                                    num_cores)
    return result_out

def run_pbsa(folder, threshold, maxiter, file_id='_difference.csv',  outfolder=None, num_cores=1,
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
        outfolder : str, optional
            output folder for the result files. If not provided output will be saved
            to input folder.
        file_id : str, optional
            Default '_difference.csv'. Should be a identifier of the files to process.
        preliminary_optional : dict, optional
            Optional parameters for preliminary step detection, details in the doc of pbsa_file.
        filter_optional : dict, optional
            Optional parameters for the exclusion of traces from the analysis, details in the doc of pbsa_file.
        refinement_optional : dict, optional
            Optional parameters for the improvement of the preliminary result, details in the doc of pbsa_file.
    '''
    folder = os.path.normpath(folder)
    if outfolder is None:
        outfolder = os.path.normpath(folder) + '/'
    else:
        outfolder = os.path.normpath(outfolder) + '/'
    os.makedirs(outfolder, exist_ok=True)
    files = sorted(glob(folder + '/*' + file_id))
    for file in files:
        pbsa_file(file, threshold, maxiter, outfolder, num_cores,
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
            # maximum fluors in the trace
            for typ in ['fluors_intensity', 'fluors_kv', 'fluors_avgintensity', 'fluors_full']:
                data_typ = np.array(data.loc[data['type'] == typ, '0':]).astype(int)
                df_file[typ] = np.max(data_typ, 1)
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
    
    

