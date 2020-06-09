#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import numpy as np
from glob import glob
import os
import tifffile


def save_with_metadata(file, data, TiffFile):
    tags = TiffFile.pages[0].tags
    taglist = []
    for k in tags:
        if k not in ('IJMetadata', 'MicroManagerMetadata'):
            if isinstance(tags[k].value, str):
                tags[k].value = tags[k].value.encode('ascii', 'ignore').decode('ascii')
            tagtuple = (tags[k].code, tags[k].dtype.strip('1'),
                        tags[k].count, tags[k].value, True)
            taglist.append(tagtuple)
    tifffile.imwrite(file, data,
                     extratags=taglist)
    return


def generate_avgimages(folder, avgover=5):
    ''' run through a folder with tiffstacks and create average image of
    1st images and last image in subfolder.
    folder: folder with tiffstacks
    avgover: Number of images to average in beginning of stack (default 5)
    '''
    folder = os.path.normpath(folder) + '/'
    outfolder = folder + 'avgimages/'
    try:
        os.mkdir(outfolder)
    except:
        pass
    filelist_out = []
    for stack in glob(folder + '*.tif'):
        fdir, fname = os.path.split(stack)
        fname = fname.split('.')[0]
        image = tifffile.TiffFile(stack)
        # average image
        image_avg = image.asarray(range(avgover))
        image_avg = np.mean(image_avg, 0)
        image_avg = np.round(image_avg).astype('uint16')
        save_with_metadata(outfolder + fname + '_avg.tif',
                           image_avg, image)
        # last image
        nframes = len(image.pages)
        image_last = image.asarray(nframes-1)
        save_with_metadata(outfolder + fname + '_last.tif',
                           image_last, image)      
    return filelist_out



