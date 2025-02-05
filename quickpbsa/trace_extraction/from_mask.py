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
from scipy.ndimage import label

# relative imports
from .helpers import pardict_from_image
from ..helpers import export_csv
from .helpers import export_mask
from .helpers import traces_from_stack

def extract_traces_mask(tiffstack, maskfile, dist, r_bg, roifile=None,
                        range_area=None, range_size=None, binning=1):
    ''' Extract bleach and background traces from tiff stack based on mask image
    
    Parameters
    __________
    
    tiffstack : str
        Path to input file (should be a tiff stack)
    maskfile : str
        Path to mask file (should be 8bit .tif, with white selection on black background)
    dist : float
        Distance between the selextion and the background ring [pix]
    r_bg : float
        Width of the background ring [pix]
    roifile : str, optional
        8 bit image to select region of interest (default None)
    range area : list, optional
        Limits for the area of selections [pix]
    range_size : list, optional
        Limits for the length of the major axis of selections [pix]
    binning : int, optional
        binning in time (frames) (default 1)
        
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
    pardict.update({'dist': dist,
                    'r_bg': r_bg, \
                    'range_area': range_area, \
                    'range_size': range_size, \
                    'binning': binning})

    # Read in mask(s)
    mask = tifffile.imread(maskfile)
    maskpix = np.where(np.ravel(mask))
    if roifile is not None:
        roi = tifffile.imread(roifile).astype(bool)
        mask = mask & roi
    # split segmentation
    labels, Nlabels = label(mask)
    labels = np.ravel(labels)
    # create constant array for surrounding
    xx, yy = np.meshgrid(np.arange(-dist, dist + 1),
                         np.arange(-dist*pardict['pix_x'], pardict['pix_x']*dist + 1, pardict['pix_x']))
    const_ex = np.ravel(xx + yy)
    dbg = dist + r_bg
    xx, yy = np.meshgrid(np.arange(-dbg, dbg + 1),
                         np.arange(-dbg*pardict['pix_x'], pardict['pix_x']*dbg + 1, pardict['pix_x']))
    const_bg = np.ravel(xx + yy)
    
    #### Create arrays for peak selection, exclusion zone, and bg selection ###
    peakpix = []
    bgpix = []
    peaksel = []
    exsel = []
    bgsel = []
    majax = []
    for I in range(1, Nlabels + 1):
        peak = np.where(labels == I)[0]
        # exclude hot pixels
        if len(peak) > 1:
            peaksel.append(peak)
            peakpix.append(len(peak))
            # exclusion zone around roi
            ex = np.unique(np.repeat(peak, len(const_ex)) + np.tile(const_ex, len(peak)))
            # bg zone (including roi)
            bg = np.unique(np.repeat(peak, len(const_bg)) + np.tile(const_bg, len(peak)))
            exsel.append(ex)
            bgsel.append(bg)
            bgpix.append(len(bg))
            # calculate major axis length of ROI (for size filtering)
            if range_size is not None:
                pcomb1 = np.repeat(peak, len(peak))
                pcomb2 = np.tile(peak, len(peak))
                dx = (pcomb1 % pardict['pix_x']) - (pcomb2 % pardict['pix_x'])
                dy = np.floor(pcomb1 / pardict['pix_x']) - np.floor(pcomb2 / pardict['pix_x'])
                dmax = np.max(np.sqrt(dx**2 + dy**2))
                majax.append(dmax)
    # convert to arrays
    peakpix = np.array(peakpix)
    majax = np.array(majax)
    maxpeakpix = np.max(peakpix)
    maxbgpix = np.max(bgpix)
    peaksel = np.vstack([np.pad(p, [0, maxpeakpix - len(p)]) for p in peaksel])
    exsel = np.hstack(exsel)
    bgsel = np.vstack([np.pad(p, [0, maxbgpix - len(p)]) for p in bgsel])
    
    #### clear background #####################################################
    # remove intersection with other exclusion zones from bg
    bgsel = np.ravel(bgsel)    
    inter, arr1, arr2 = np.intersect1d(bgsel, exsel, return_indices = True)
    bgsel[arr1] = 0
    while len(inter > 0):
        inter, arr1, arr2 = np.intersect1d(bgsel, exsel, return_indices = True)
        bgsel[arr1] = 0
    # remove remaining mask pixels from bg
    inter, arr1, arr2 = np.intersect1d(bgsel, maskpix, return_indices = True)
    bgsel[arr1] = 0
    while len(inter > 0):
        inter, arr1, arr2 = np.intersect1d(bgsel, maskpix, return_indices = True)
        bgsel[arr1] = 0
    bgsel = np.reshape(bgsel, [np.size(peaksel, 0),  - 1])
    
    # remove edges
    bgsel[np.floor(bgsel%pardict['pix_x']) < (dbg + 1)] = 0
    bgsel[np.floor(bgsel%pardict['pix_x']) > (pardict['pix_x'] - dbg - 1)] = 0
    bgsel[np.floor(bgsel/pardict['pix_x']) < (dbg + 1)] = 0
    bgsel[np.floor(bgsel/pardict['pix_x']) > (pardict['pix_y'] - dbg - 1)] = 0
    
    bgpix = np.sum(bgsel > 0, 1)
        
    #### filter ###############################################################    
    bgfilter = bgpix > 0
    #remove points without bg
    peaksel = peaksel[bgfilter, :]
    bgsel = bgsel[bgfilter, :]
    peakpix = peakpix[bgfilter]
    bgpix = bgpix[bgfilter]   
    if range_size is not None:
        majax = majax[bgfilter]      
    #filter over area
    if range_area is not None:
        sel = (peakpix >= range_area[0]) & (peakpix <= range_area[1])
        peaksel = peaksel[sel, :]
        bgsel = bgsel[sel, :]
        peakpix = peakpix[sel]
        bgpix = bgpix[sel]
        majax = majax[sel]        
    # filter over major axis length
    if range_size is not None:
        sel = (majax > range_size[0]) & (majax < range_size[1])
        peaksel = peaksel[sel, :]
        bgsel = bgsel[sel, :]
        peakpix = peakpix[sel]
        bgpix = bgpix[sel]
        majax = majax[sel]
        
    #### save peak and background selection masks #############################
    outpath, fn = os.path.split(maskfile)
    p, fname = os.path.split(tiffstack)
    fname = fname.split('.')[0]
    basedf = pd.DataFrame()
    basedf['id'] = np.arange(len(bgpix))
    basedf['peak_pix'] = peakpix
    basedf['bg_pix'] = bgpix
    # write peak selection to file
    outfile_peaksel = outpath + '/' + fname + '_peaksel.csv'
    export_csv(basedf, peaksel.T, outfile_peaksel, pardict)
    # write background selection to file
    outfile_bgsel = outpath + '/' + fname + '_bgsel.csv'
    export_csv(basedf, bgsel.T, outfile_bgsel, pardict)
    # save peak selection as 8bit image
    outfile_peakmask = outpath + '/' + fname + '_mask_peak.tif'
    export_mask(peaksel, outfile_peakmask, pardict)
    # save bg selection as 8bit image
    outfile_bgmask = outpath + '/' + fname + '_mask_bg.tif'
    export_mask(bgsel, outfile_bgmask, pardict)
    
    #### Extract traces from stack and save to csv ############################
    peak, bg = traces_from_stack(tiffstack, peaksel.T, bgsel.T, peakpix, bgpix, binning)
    outfile_peak = outpath + '/' + fname + '_peak.csv'
    export_csv(basedf, peak.T, outfile_peak, pardict)
    outfile_bg = outpath + '/' + fname + '_bg.csv'
    export_csv(basedf, bg.T, outfile_bg, pardict)
    difference = peak - bg
    outfile_difference = outpath + '/' + fname + '_difference.csv'
    export_csv(basedf, difference.T, outfile_difference, pardict)
    
    return peak, bg, difference
