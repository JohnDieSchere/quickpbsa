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
    folder = os.path.normpath(folder)
    files = sorted(glob(folder + '/*.tif'))
    for file in files:
        fp, filename = os.path.split(file)
        filename, ext = os.path.splitext(filename)
        locfile = folder + '/' + resultsfolder + '/' + filename + '_' + locfile_id
        roifile = folder + '/' + resultsfolder + '/' + filename + '_' + 'roi.tif'
        if not os.path.isfile(roifile):
            roifile = None
        if os.path.isfile(locfile):
            extract_traces_localization(file, locfile, r_peak, r_bg1, r_bg2, min_dist,
                                        roifile=roifile, filters=filters,
                                        binning=binning)            
    return


def run_extraction_mask(folder, dist, r_bg,
                        resultsfolder='results', maskfile_id='mask.tif',
                        binning=1, range_area=None, range_size=None):
    folder = os.path.normpath(folder)
    files = sorted(glob(folder + '/*.tif'))
    for file in files:
        fp, filename = os.path.split(file)
        filename, ext = os.path.splitext(filename)
        maskfile = folder + '/' + resultsfolder + '/' + filename + '_' + maskfile_id
        roifile = folder + '/' + resultsfolder + '/' + filename + '_' + 'roi.tif'
        if not os.path.isfile(roifile):
            roifile = None
        if os.path.isfile(maskfile):
            extract_traces_mask(file, maskfile, dist, r_bg, roifile=roifile,
                                range_area=range_area, range_size=range_size)            
    return