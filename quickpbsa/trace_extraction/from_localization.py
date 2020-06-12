#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import numpy as np
import pandas as pd
import tifffile
import os

# relative imports
from .helpers import pardict_from_image
from .helpers import export_mask
from .helpers import traces_from_stack
from ..helpers import export_csv

def extract_traces_localization(tiffstack, locfile, r_peak, r_bg1, r_bg2, min_dist,
                                roifile=None, filters={}, binning=1, pix_size=None):
    ''' Extract bleach and background traces from tiff stack based on localization output (e.g. ThunderSTORM)
    and write them into .csv files
    
    Parameters
    __________
    
    tiffstack : str
        Path to input file (should be a tiff stack)
    locfile : str
        Path to localization file (should be .csv, containing the columns 'x [nm]' and 'y [nm]')
    r_peak : float
        peak radius [pix]
    r_bg1 : float
        inner bg radius [pix]
    r_bg2 : float
        outer bg radius [pix]
    min_dist: float
        minimum distance between peaks [pix]
    roifile : str, optional
        8 bit image to select region of interest (default None)
    filters : dict, optional
        Filter localizations by any column in the localization file, syntax is e.g. filters={'sigma': [40,60]}
    binning : int, optional
        binning in time (frames) (default 1)
    pix_size : float, optional
        Define pixel size in nm. In most cases this is not necessary, since the
        pixel size is read out from the tiff tags. (default None)
        
    Returns
    _______
    
    peak : array
        Peak traces
    bg : array
        Background traces
    difference : array
        Difference Traces
    '''
    
    # generate parameter dictionary
    pardict = pardict_from_image(tiffstack)
    if pix_size is not None:
        pardict['pix_size'] = pix_size
    pardict.update({'r_peak': r_peak,
                    'r_bg1': r_bg1,
                    'r_bg2': r_bg2,
                    'min_dist': min_dist,
                    'binning': binning})
    for k in filters:
        pardict[k + '_range'] = filters[k]

    # Read in coordinate file
    loc_df = pd.read_csv(locfile)
    # ERROR HERE IF NO x [nm] IN loc_df
    x = np.array(loc_df['x [nm]']) / pardict['pix_size']
    y = np.array(loc_df['y [nm]']) / pardict['pix_size']
    # exclude spots to close to the edge
    edge = int(np.ceil(r_bg2))
    coordsel = (x > edge) & (x < (pardict['pix_x'] - edge - 1))
    coordsel = coordsel & (y > edge) & (y < (pardict['pix_y'] - edge - 1))
    x = x[coordsel]
    y = y[coordsel]
    loc_df = loc_df.iloc[coordsel, :]
    
    #### build arrays of pixels surrounding centers ###########################
    # format to 1d
    center = pardict['pix_x'] * np.floor(y).astype(int) + np.floor(x).astype(int)
    loc_df['center_pix'] = center    
    # create a surrounding points array (respective to center)
    radind = np.tile(np.arange( -edge * pardict['pix_x'],
                               (edge + 1) * pardict['pix_x'],
                               pardict['pix_x']),
                     [2*edge + 1, 1])
    radind = radind.T + np.arange( -edge, edge + 1)    
    # surrounding indices for each point
    ind = (np.tile(np.ravel(radind), [len(center), 1]).T + center).astype(int)    
    # pixel coordinates from ind
    x_ar = ind % pardict['pix_y'] + 0.5
    y_ar = np.floor(ind / pardict['pix_x']) + 0.5    
    # distance to center array
    dist = (x_ar - x)**2 + (y_ar - y)**2
    
    #### Remove points with too close neighbours ##############################
    # repeated center coordinates
    x_rep, y_rep = np.meshgrid(x, y)
    # distances from center to center
    point_distance = (x_rep.T - x)**2 + (y_rep - y)**2
    # boolean array over distance
    selector = point_distance < min_dist**2
    # ignore diagonal
    selector = selector & np.invert(np.eye(np.size(center)).astype(bool))
    # cleaner for density filtering
    distfilter = np.where(np.sum(selector, 0) == 0)[0]
    if np.size(distfilter) < 2:
        pass
        # PUT ERROR HERE
    # peak selector
    peaksel = ind.copy()
    peaksel[dist > r_peak**2] = 0
    peaksel = peaksel[:, distfilter]
    peakpix = np.sum(peaksel > 0, 0)
    #background selector
    bgsel = ind.copy()
    bgsel[dist > r_bg2**2] = 0
    bgsel[dist < r_bg1**2] = 0
        
    #### remove intersection with peaks from background #######################    
    # remove intersection with rbg_1 from other points
    bgsel = np.ravel(bgsel)
    bg1sel = ind[dist < r_bg1**2]
    intersect, intersectind1, intersectind2 = np.intersect1d(bgsel, bg1sel,
                                                             return_indices=True)
    while len(intersect) > 0:
        bgsel[intersectind1] = 0
        intersect, intersectind1, intersectind2 = np.intersect1d(bgsel, bg1sel,
                                                                 return_indices=True)
    # reshape and remove distance filtered points
    bgsel = np.reshape(bgsel, np.shape(ind))
    bgsel = bgsel[:, distfilter]
    # filter for points without bg roi
    bgpix = np.sum(bgsel > 0, 0)
    bgfilter = (bgpix > 0)
    peaksel = peaksel[:, bgfilter]
    peakpix = peakpix[bgfilter]
    bgsel = bgsel[:, bgfilter]
    bgpix = bgpix[bgfilter]
    # filter loc_df and center
    loc_df = loc_df.iloc[distfilter, :]
    loc_df = loc_df.iloc[bgfilter, :]
    center = center[distfilter]
    center = center[bgfilter]
    
    #### filter according to ROI mask #########################################
    if roifile is not None:
        # assumes mask is 8bit image with white selection
        roi = np.ravel(tifffile.imread(roifile))
        roi_ind = np.where(roi)[0]
        intersect, ind1, ind2 = np.intersect1d(center, roi_ind,
                                               return_indices=True)
        center = center[ind1]
        loc_df = loc_df.iloc[ind1, :]
        peaksel = peaksel[:, ind1]
        peakpix = peakpix[ind1]
        bgsel = bgsel[:, ind1]
        bgpix = bgpix[ind1]
    
    #### optional filtering ###################################################
    for k in filters:
        data_k = np.array(loc_df[k])
        filter_selection = (data_k > filters[k][0]) & (data_k < filters[k][1])
        loc_df = loc_df.iloc[filter_selection, :]
        peaksel = peaksel[:, filter_selection]
        peakpix = peakpix[filter_selection]
        bgsel = bgsel[:, filter_selection]
        bgpix = bgpix[filter_selection]
        
    #### save peak and background selection masks #############################
    outpath, fn = os.path.split(locfile)
    p, fname = os.path.split(tiffstack)
    fname, ext = os.path.splitext(fname)
    loc_df = loc_df.reset_index(drop = True)
    # write peak selection to file
    outfile_peaksel = outpath + '/' + fname + '_peaksel.csv'
    export_csv(loc_df, peaksel.T, outfile_peaksel, pardict)
    # write background selection to file
    outfile_bgsel = outpath + '/' + fname + '_bgsel.csv'
    export_csv(loc_df, bgsel.T, outfile_bgsel, pardict)
    # save peak selection as 8bit image
    outfile_peakmask = outpath + '/' + fname + '_mask_peak.tif'
    export_mask(peaksel, outfile_peakmask, pardict)
    # save bg selection as 8bit image
    outfile_bgmask = outpath + '/' + fname + '_mask_bg.tif'
    export_mask(bgsel, outfile_bgmask, pardict)
    
    #### Extract traces from stack and save to csv ############################
    peak, bg = traces_from_stack(tiffstack, peaksel, bgsel, peakpix, bgpix, binning)
    outfile_peak = outpath + '/' + fname + '_peak.csv'
    export_csv(loc_df, peak.T, outfile_peak, pardict)
    outfile_bg = outpath + '/' + fname + '_bg.csv'
    export_csv(loc_df, bg.T, outfile_bg, pardict)
    difference = peak - bg
    outfile_difference = outpath + '/' + fname + '_difference.csv'
    export_csv(loc_df, difference.T, outfile_difference, pardict)
    
    return peak, bg, difference
