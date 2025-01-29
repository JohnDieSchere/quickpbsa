#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import numpy as np
import tifffile



def pardict_from_image(imfile):
    pardict = {}
    with tifffile.TiffFile(imfile) as tif:
        # image size
        pardict['pix_x'] = tif.pages[0].asarray().shape[1]
        pardict['pix_y'] = tif.pages[0].asarray().shape[0]
        resolution_unit = int(tif.pages[0].tags['ResolutionUnit'].value)
        resolution = tif.pages[0].tags['XResolution'].value
        if resolution_unit == 3:
            # resolution in cm
            pardict['pix_size'] = 1e7 / resolution[0]
        elif resolution_unit == 1:
            # no unit defined, usually ratio of resolution tag values
            pardict['pix_size'] = 1e3 * resolution[1] / resolution[0]
    return pardict

def export_mask(indices, filename, pardict):
    mask = np.ones(pardict['pix_x']*pardict['pix_y'], dtype = 'uint8')*255
    mask[indices] = 0
    mask = np.reshape(mask, [pardict['pix_y'], pardict['pix_x']])
    tifffile.imwrite(filename, data=mask)
    return

def traces_from_stack(stack, peaksel, bgsel, peakpix, bgpix, binning):
    # read in stack
    data = tifffile.imread(stack)
    # ravel x and y
    data = data.reshape(data.shape[0], -1)
    # set 0,0 pixel to 0, since that is the coordinate to be sorted out
    data[:, 0] = 0
    # peak traces
    peak = data[:, peaksel]
    peak = np.sum(peak, 1) / peakpix
    # bg traces
    bg = data[:, bgsel]
    bg = np.sum(bg, 1) / bgpix
    # binning
    if binning > 1:
        # number of bins
        nbins = peak.shape[0]//binning
        # peak
        peak = peak[:binning*nbins, :]
        peak = peak.reshape(nbins, -1, peak.shape[1])
        peak = np.sum(peak, 1) / binning
        # bg
        bg = bg[:binning*nbins, :]
        bg = bg.reshape(nbins, -1, bg.shape[1])
        bg = np.sum(bg, 1) / binning
    return peak, bg