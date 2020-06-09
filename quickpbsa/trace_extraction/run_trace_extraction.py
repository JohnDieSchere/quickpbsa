#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import os
from glob import glob

# relative imports
from .from_localization import extract_traces_localization
from .from_mask import extract_traces_mask


def run_extraction_localization(folder, r_peak, r_bg1, r_bg2, min_dist,
                                resultsfolder='results', locfile_id='tsoutput.csv',
                                filters={}, binning=1, pix_size=None):
    ''' Extract bleach and background traces from all tiff stacks in a folder based
    on localization output (e.g. ThunderSTORM) and write them into .csv files.
    
    Parameters
    __________
    
    folder : str
        Folder where tiff stacks are located
    r_peak : float
        peak radius [pix]
    r_bg1 : float
        inner bg radius [pix]
    r_bg2 : float
        outer bg radius [pix]
    min_dist: float
        minimum distance between peaks [pix]
    resultsfolder : str, optional
        Name of subfolder where localization files are located and where traces will be saved (default results)\
        ROI images in the resultsfolder named <filename>_roi.tif will be used.
    locfile_id : str, optional
        identifier for localization files (default 'tsoutput.csv')
    filters : dict, optional
        Filter localizations by any column in the localization file, syntax is e.g. filters={'sigma': [40,60]}
    binning : int, optional
        binning in time (frames) (default 1)
    pix_size : float, optional
        Define pixel size in nm. In most cases this is not necessary, since the
        pixel size is read out from the tiff tags. (default None)
    '''
    
    folder = os.path.normpath(folder)
    files = sorted(glob(folder + '/*.tif'))
    for file in files:
        fp, filename = os.path.split(file)
        filename, ext = os.path.splitext(filename)
        locfile = folder + '/' + resultsfolder + '/' + filename + '_' + locfile_id
        roifile = folder + '/' + resultsfolder + '/' + filename + '_roi.tif'
        if not os.path.isfile(roifile):
            roifile = None
        if os.path.isfile(locfile):
            extract_traces_localization(file, locfile, r_peak, r_bg1, r_bg2, min_dist,
                                        roifile, filters, binning)
    return


def run_extraction_mask(folder, dist, r_bg,
                        resultsfolder='results', maskfile_id='mask.tif',
                        range_area=None, range_size=None, binning=1):
    ''' Extract bleach and background traces from all tiff stacks in folder based on mask images
    
    Parameters
    __________
    
    folder : str
        Folder where tiff stacks are located
    dist : float
        Distance between the selextion and the background ring [pix]
    r_bg : float
        Width of the background ring [pix]
    resultsfolder : str, optional
        Name of subfolder where localization files are located and where traces will be saved (default results)\
        ROI images in the resultsfolder named <filename>_roi.tif will be used.
    maskfile_id : str, optional
        identifier for mask files (default 'mask.tif')
    range area : list, optional
        Limits for the area of selections [pix]
    range_size : list, optional
        Limits for the length of the major axis of selections [pix]
    binning : int, optional
        binning in time (frames) (default 1)
    '''
    
    folder = os.path.normpath(folder)
    files = sorted(glob(folder + '/*.tif'))
    for file in files:
        fp, filename = os.path.split(file)
        filename, ext = os.path.splitext(filename)
        maskfile = folder + '/' + resultsfolder + '/' + filename + '_' + maskfile_id
        roifile = folder + '/' + resultsfolder + '/' + filename + '_roi.tif'
        if not os.path.isfile(roifile):
            roifile = None
        if os.path.isfile(maskfile):
            extract_traces_mask(file, maskfile, dist, r_bg,
                                roifile, range_area, range_size, binning)            
    return