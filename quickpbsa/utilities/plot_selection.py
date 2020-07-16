#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tifffile


def outline_boolean(data):
    '''
    outline the True values
    '''
    data = np.pad(data, 1, mode = 'constant').astype(bool)
    pix = np.shape(data)
    ## edge detect horizontal
    horizontal = (data[1:, :]^data[: - 1, :])
    edges = np.where(np.diff(np.ravel(horizontal)))[0]
    ##convert to plot with lines
    lends = (edges%pix[1]) + 1
    y = np.repeat(np.arange(pix[0] - 1), pix[1])[edges] + 1
    y = y[::2] - 1
    hstarts = lends[::2] - 1
    hends = lends[1::2] - 1
    ## edge detect vertical
    data = np.transpose(data)
    vertical = (data[1:, :]^data[: - 1, :])
    edges = np.where(np.diff(np.ravel(vertical)))[0]
    ## convert to plot
    lends = (edges%pix[0]) + 1
    x = np.repeat(np.arange(pix[1] - 1), pix[0])[edges] + 1
    x = x[::2] - 1
    vstarts = lends[::2] - 1
    vends = lends[1::2] - 1
    return y, hstarts, hends, x, vstarts, vends


def selection_svg_fromtif(file, fig=None, ax=None, outfile=None, color='w'):
    data = tifffile.imread(file).astype(bool)
    data = np.invert(data)
    pix = np.shape(data)
    y, hstarts, hends, x, vstarts, vends = outline_boolean(data)
    # set edges to true to avoid edge issues
    if ax is None:
        fig = plt.figure()
        ax = fig.add_axes([0.0, 0.0, pix[1]/np.max(pix), pix[0]/np.max(pix)])
    ax.hlines(y - .5, hstarts - .5, hends - .5, color=color)
    ax.vlines(x - .5, vstarts - .5, vends - .5, color=color)
    # plot
    if ax is None:
        ax.set_facecolor([0, 0, 0, 0])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim([pix[0], 0])
        ax.set_xlim([0, pix[1]])
        fig.set_size_inches([5, 5])
    if outfile is not None:
        plt.savefig(outfile)
    return fig, ax



def selection_svg_fromcsv(file, ids, outfile = None, fig = None, ax = None, colors='w'):
    # read in csv data
    df = pd.read_csv(file, header = 1)
    comment = open(file).readline()
    pardict = eval(comment[comment.find('{'):])
    pix = [pardict['pix_x'], pardict['pix_y']]
    if type(ids) is int:
        ids = [ids, ]
    newfig = False
    if fig is None:
        fig = plt.figure()
        ax = fig.add_axes([0.0, 0.0, pix[1]/np.max(pix), pix[0]/np.max(pix)])
        newfig = True
    i = 0
    if isinstance(colors, str):
        colors = [colors,] * len(ids)
    for ID in ids:
        center = int(df.loc[df['id'] == ID, 'center_pix'])
        ar = np.squeeze(np.array(df.loc[df['id'] == ID, '0':]).astype(bool))
        size = int(np.sqrt(len(ar)))
        ar = ar.reshape([size, size])
        y, hstart, hend, x, vstart, vend = outline_boolean(ar)
        cy = int(center/pix[0])
        cx = center%pix[1]
        yoff = cy - int(size/2) - .5
        xoff = cx - int(size/2) - .5
        ax.hlines(y + yoff, hstart + xoff, hend + xoff, color=colors[i])
        ax.vlines(x + xoff, vstart + yoff, vend + yoff, color=colors[i])
        i += 1
    # plot
    if newfig:
        ax.set_facecolor([0, 0, 0, 0])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim([pix[0], 0])
        ax.set_xlim([0, pix[1]])
        fig.set_size_inches([5, 5])
    if outfile is not None:
        plt.savefig(outfile)
    return fig, ax



