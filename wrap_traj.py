#!/usr/bin/env python -u
import math
import numpy as np
import matplotlib.pyplot as plt

def wrap(   
    Zxyzs,
    unit_cell,
    ):

    xhi = unit_cell[0,1] - unit_cell[1,0]
    yhi = unit_cell[1,1] - unit_cell[2,0]
    zhi = unit_cell[2,1] - unit_cell[3,0]

    new_Zxyzs = np.zeros_like(Zxyzs)
    for ind,Zxyz in enumerate(Zxyzs):
        new_Zxyzs[ind,:,1] = np.mod(Zxyz[:,1] / (xhi),1) * xhi 
        new_Zxyzs[ind,:,2] = np.mod(Zxyz[:,2] / (yhi),1) * yhi 
        new_Zxyzs[ind,:,3] = np.mod(Zxyz[:,3] / (zhi),1) * zhi 

    return new_Zxyzs


def make_whole(
    Zxyzs,
    ):

    



    
