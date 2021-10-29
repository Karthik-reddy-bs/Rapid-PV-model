# -*- coding: utf-8 -*-
"""
This code has been taken from 
"""

from __future__ import print_function
from builtins import range

import numpy as np
from shadowingfunction_wallheight_13 import shadowingfunction_wallheight_13
from shadowingfunction_wallheight_23 import shadowingfunction_wallheight_23
import linecache
import sys



def Worker(dsm, scale, building_slope, building_aspect, voxelheight, sizey, sizex, vegdsm, vegdsm2,
             wheight, waspect, albedo, psi, radmatI, radmatD, radmatR, usevegdem, calc_month):
    a = dsm
    scale = scale
    slope = building_slope
    aspect = building_aspect
    walls = wheight
    dirwalls = waspect
    vegdem = vegdsm
    vegdem2 = vegdsm2
    # Parameters
    deg2rad = np.pi/180
    Knight = np.zeros((sizex, sizey))
    Energyyearroof = np.copy(Knight)
    Direct = np.copy(Knight)
    Diffuse = np.copy(Knight)
    reflected = np.copy(Knight)



    if usevegdem == 1:
        # amaxvalue
        vegmax = vegdem.max()
        amaxvalue = a.max() - a.min()
        amaxvalue = np.maximum(amaxvalue, vegmax)

        # Elevation vegdsms if buildingDEM includes ground heights
        vegdem = vegdem+a
        vegdem[vegdem == a] = 0
        vegdem2 = vegdem2+a
        vegdem2[vegdem2 == a] = 0

        #% Bush separation
        bush = np.logical_not((vegdem2*vegdem))*vegdem
    else:
        psi = 1

    # Creating wallmatrix (1 meter interval)
    wallcol, wallrow = np.where(np.transpose(walls) > 0)    # row and col for each wall pixel
    # wallrow, wallcol = np.where(walls > 0.2)    #row and col for each wall pixel
    wallstot = np.floor(walls * (1 / voxelheight)) * voxelheight
    wallsections = np.floor(np.max(walls) * (1 / voxelheight))     # finding tallest wall
    wallmatrix = np.zeros((np.shape(wallrow)[0], int(wallsections)))
    Energyyearwall = np.copy(wallmatrix)
    # Energymonthwall = np.zeros(np.shape(wallmatrix[0]), np.shape(wallmatrix[1]), 12)

    # Main loop - Creating skyvault of patches of constant radians (Tregeneza and Sharples, 1993)
    skyvaultaltint = np.array([6, 18, 30, 42, 54, 66, 78, 90])
    aziinterval = np.array([30, 30, 24, 24, 18, 12, 6, 1])

    if usevegdem == 1:
        wallshve = np.zeros(np.shape(a))
        vegrow, vegcol = np.where(vegdem > 0)  # row and col for each veg pixel
        vegdata = np.zeros((np.shape(vegrow)[0], 3))
        for i in range(0, vegrow.shape[0] - 1):
            vegdata[i, 0] = vegrow[i] + 1
            vegdata[i, 1] = vegcol[i] + 1
            vegdata[i, 2] = vegdem[vegrow[i], vegcol[i]]
    else:
        vegdata = 0

    index = 0
    for i in range(skyvaultaltint.size):
        for j in range(aziinterval[i]):

            #################### SOLAR RADIATION POSITIONS ###################
            # Solar Incidence angle (Roofs)
            suniroof = np.sin(slope) * np.cos(radmatI[index, 0] * deg2rad) * \
                        np.cos((radmatI[index, 1]*deg2rad)-aspect) + \
                        np.cos(slope) * np.sin((radmatI[index, 0] * deg2rad))

            suniroof[suniroof < 0] = 0

            # Shadow image
            if usevegdem == 1:
                vegsh, sh, _, wallsh, wallsun, wallshve, _, facesun = shadowingfunction_wallheight_23(a,
                                    vegdem, vegdem2, radmatI[index, 1], radmatI[index, 0], scale, amaxvalue,
                                                                                bush, walls, dirwalls * deg2rad)
                shadow = np.copy(sh-(1.-vegsh)*(1.-psi))
            else:
                sh, wallsh, wallsun, facesh, facesun = shadowingfunction_wallheight_13(a, radmatI[index, 1],
                                                            radmatI[index, 0], scale, walls, dirwalls * deg2rad)
                shadow = np.copy(sh)

            # roof irradiance calculation
            # direct radiation
            if radmatI[index, 2] > 0:
                I = shadow * radmatI[index, 2] * suniroof
            else:
                I = np.copy(Knight)

            # roof diffuse and reflected radiation
            D = radmatD[index, 2] * shadow
            R = radmatR[index, 2] * (shadow*-1 + 1)
            Direct = np.copy(Direct+I)
            Diffuse = np.copy(Diffuse+D)
            reflected = np.copy(reflected+R)

            Energyyearroof = np.copy(Energyyearroof+D+R+I)

        
            index = index + 1


    return Energyyearroof, Direct, Diffuse, reflected

