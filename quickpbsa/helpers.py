#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:29:12 2019

@author: johan
"""

import pandas as pd



def export_csv(base_df, array, filename, pardict, comment='# '):
    fid = open(filename, 'w')
    fid.write(comment + str(pardict) + '\n')
    df = pd.concat([base_df, pd.DataFrame(array)], axis = 1)
    df.to_csv(fid, index = False)
    fid.close()
    return df