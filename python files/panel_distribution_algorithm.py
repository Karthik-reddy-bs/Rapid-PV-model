# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 19:05:09 2021

@author: karthik reddy

Panel Distribution Algorithm: is used to geometrically place PV panels
on rooftop segments considering the slope, aspect and geometry
of rooftop segments.
"""
import os
import ogr
from math import ceil
import numpy as np
import pandas as pd
import geopandas as gpd
import math
import shapely
import operator
from itertools import chain
# import time
from shapely.geos import TopologicalError

"""
    Main Input parameters for the PDA:
    -----------------
    The primary input for this algorithm is an attribute table in the form
    of a shapefile, which consists of rooftop segment ID, geometry of
    rooftop segment, slope and aspect of the rooftop segment
    
    crs          : co-ordinate reference system of the shapefile
    row_distance : Module row spacing, can be optimized according to technical & Economical parameters
    flat_slope   : Optimal slope of PV panels on flat rooftops, can be optimized according to 
                   technical & Economical parameters
    panel_length : Length of PV panel, can be changed according to the requirement
    panel_width  : Width of PV panel, can be changed according to the requirement
               
    """

shapefile = gpd.read_file("segments.shp")
crs = "EPSG:25833"
row_distance = 1.5
flat_slope = 40
panel_length = 1.6
panel_width = 1

# Temporary and output locations of the shapefiles.
temp_loc = "temp.shp"
output_loc = "pv_panels.shp"

# split dataframe by row and get the crs
splits = [shapefile.loc[[i]] for i in shapefile.index]
gdf = gpd.GeoDataFrame(crs=crs)


def fishnet(poly, count, aspect_polygon, mode, slope, row_distance, num=0):
    
    """
    Input Parameters:
    -----------------
    
        poly = Rooftop segment polygon
        count = To count the PV panels
        aspect_polygon = Orientation of the rooftop segment (0 - 360 deg)
        mode = Landscape or portrait mode
        slope = Slope of rooftop segment
        row_distance = module row spacing between the rows of PV arrays.
        num = A temporary constant (0 - default, 1 - portrait, 2- landscape)
        
    Returns
    -------
        PV panels geometries on rooftop segment.
        
    """
# To get bounds of the polygon for creating the fishnet grid
    bounds2 = poly.bounds
    xmin, ymin, xmax, ymax = bounds2
    
    
    def main(outputGridfn, xmin, xmax, ymin, ymax, gridHeight, gridWidth, slope, row_distance):
        
        """
         Input Parameters:
         ----------------
        
            outputGridfn = Location of output file.
            xmin, xmax, ymin, ymax = Bounds of the rooftop segment polygon
            gridHeight = Length of PV panel
            gridWidth = Width of PV panel
            slope = Slope of rooftop segment
            row_distance = module row spacing between the rows of PV arrays.
            
         Returns
         -------
            Fishnet grid according to the slope and module row spacing of the rooftop segment.
            
        """
    
        # convert bounds to float
        diff1 = float(xmax) - float(xmin)
        diff2 = float(ymax) - float(ymin)
        
        # Creating the boundary rectangle for fishnet grid
        if diff1 > diff2:
            
            xmin = float(xmin) - 5
            xmax = float(xmax) + 5
            ymin = float(ymin) - (diff1)
            ymax = float(ymax) + (diff1)
            
        else:
            xmin = float(xmin) - (diff2)
            xmax = float(xmax) + (diff2)
            ymin = float(ymin) - 5
            ymax = float(ymax) + 5
    
        gridWidth = float(gridWidth)
        gridHeight = float(gridHeight)
        
        # get rows
        rows = ceil((ymax-ymin)/gridHeight)
        # get columns
        cols = ceil((xmax-xmin)/gridWidth)
        
        # start grid cell envelope
        ringXleftOrigin = xmin
        ringXrightOrigin = xmin + gridWidth
        ringYtopOrigin = ymax
        ringYbottomOrigin = ymax-gridHeight
    
        # create output file
        outDriver = ogr.GetDriverByName('ESRI Shapefile')
        if os.path.exists(outputGridfn):
            os.remove(outputGridfn)
        outDataSource = outDriver.CreateDataSource(outputGridfn)
        outLayer = outDataSource.CreateLayer(outputGridfn, geom_type=ogr.wkbPolygon)
        featureDefn = outLayer.GetLayerDefn()
    
        # create grid cells
        countrows = 0
        while countrows < rows:
            countrows += 1
    
            # reset envelope for rows
            ringXleft = ringXleftOrigin
            ringXright = ringXrightOrigin
            
            countcols = 0
    
            while countcols < cols:
                countcols += 1
                ring = ogr.Geometry(ogr.wkbLinearRing)
                ring.AddPoint(ringXleft, ringYtopOrigin)
                ring.AddPoint(ringXright, ringYtopOrigin)
                ring.AddPoint(ringXright, ringYbottomOrigin)
                ring.AddPoint(ringXleft, ringYbottomOrigin)
                ring.AddPoint(ringXleft, ringYtopOrigin)
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ring)
    
                # add new geom to layer
                outFeature = ogr.Feature(featureDefn)
                outFeature.SetGeometry(poly)
                outLayer.CreateFeature(outFeature)
                outFeature.Destroy
    
                # new envelope for next poly
                
                ringXleft = ringXleft + gridWidth
                ringXright = ringXright + gridWidth
    
            # new envelope for next poly
            if slope == flat_slope:
                ringYtopOrigin = ringYtopOrigin - gridHeight - row_distance
                ringYbottomOrigin = ringYbottomOrigin - gridHeight - row_distance
            else:
                ringYtopOrigin = ringYtopOrigin - gridHeight
                ringYbottomOrigin = ringYbottomOrigin - gridHeight
    
        # Close DataSources
        outDataSource.Destroy()
    
    
    
    if __name__ == "__main__":
            
        outfp = temp_loc
        area = mode[0] * mode[1]
        main(outfp, xmin, xmax, ymin, ymax, mode[1], mode[0], slope, row_distance)
    
   # stores the fishnet grid as a temporary shapefile "temp.shp" 
    temp = gpd.read_file("temp.shp")
    temp.crs = crs

    # Rotating the fishnet grid according to the aspect of the rooftop segment
    a = list(chain(*list(temp.loc[0, 'geometry'].centroid.coords)))
    b = list(chain(*list(temp.loc[temp.shape[0] - 1, 'geometry'].centroid.coords)))
    c = np.average(np.array([a, b]), axis=0)
    for index, row in temp.iterrows():    
        rotated = shapely.affinity.rotate(row['geometry'], 180 + aspect_polygon, (c[0], c[1]))
        temp.loc[index, 'geometry'] = rotated
    
    # Translates the fishnet grid horizontally to accomodate more PV panels
    opt_x = {}
    
    for x in np.arange(0,mode[1],0.1):
        
        shift_x = temp.translate(x * round(math.cos(math.radians(180 + aspect_polygon)), 2), \
                                 * round(math.sin(math.radians(180 + aspect_polygon)), 2))
        shift_x.to_file("shifttemp.shp")
        shifttemp = gpd.read_file("shifttemp.shp")
        shifttemp.crs = crs
        
        try:
            res_intersection = gpd.overlay(polygon_df, shifttemp, how='intersection')
            splits1 = [res_intersection.loc[[i]] for i in res_intersection.index]
            for i in res_intersection.index:
                polygon_df1 = splits1[i]
                poly3 = polygon_df1.loc[i, 'geometry']
                poly3.crs = crs
                res_intersection.loc[i, 'pv_area'] = poly3.area
            
            ful_polygon = res_intersection[res_intersection['pv_area'] > area-0.01]
            opt_x[x] = ful_polygon.shape[0]
            
        except TopologicalError:
            pass
    
    opt_x = max(opt_x.items(), key=operator.itemgetter(1))[0]
    shift = temp.translate(opt_x * round(math.cos(math.radians(180 + aspect_polygon)), 2),\
                           opt_x * round(math.sin(math.radians(180 + aspect_polygon)), 2))
    
    # Translates the fishnet grid vertically to accomodate more PV panels
    opt_y = {}
    for y in np.arange(0,mode[0],0.2):
        
        shift_y = shift.translate(y * round(math.cos(math.radians(90 + 180 + aspect_polygon)), 2),\
                                  y * round(math.sin(math.radians(90 + 180 + aspect_polygon)), 2))
        shift_y.to_file("shifttemp.shp")
        shifttemp = gpd.read_file("shifttemp.shp")
        shifttemp.crs = crs
        try:
            res_intersection = gpd.overlay(polygon_df, shifttemp, how='intersection')
            splits1 = [res_intersection.loc[[i]] for i in res_intersection.index]
            for i in res_intersection.index:
                polygon_df1 = splits1[i]
                poly3 = polygon_df1.loc[i, 'geometry']
                poly3.crs = crs
                res_intersection.loc[i, 'pv_area'] = poly3.area
            
            ful_polygon = res_intersection[res_intersection['pv_area'] > area-0.01]
            opt_y[y] = ful_polygon.shape[0]

        except TopologicalError:
            pass
    
    opt_y = max(opt_y.items(), key=operator.itemgetter(1))[0]
    
    # Stores the best possible location of fishnet grid in a temporary file for PV slicing
    shifted = shift.translate(opt_y * round(math.cos(math.radians(90 + 180 + aspect_polygon)), 2),\
                              opt_y * round(math.sin(math.radians(90 + 180 + aspect_polygon)), 2))
    shifted.to_file("shifttemp.shp")
    shifttemp = gpd.read_file("shifttemp.shp")
    shifttemp.crs = crs
    
    # PV panel slicing is done by intersecting the rooftop segment polygon with the fishnet grid.
    # Sometimes the vector datasets have topological errors due to intersection of multiple polygons and linestrings.
    try:
        res_intersection = gpd.overlay(polygon_df, shifttemp, how='intersection')
        splits1 = [res_intersection.loc[[i]] for i in res_intersection.index]
        for i in res_intersection.index:
            polygon_df1 = splits1[i]
            poly3 = polygon_df1.loc[i, 'geometry']
            poly3.crs = crs
            res_intersection.loc[i, 'pv_area'] = poly3.area
        
        ful_polygon = res_intersection[res_intersection['pv_area'] > area-0.01]
        count[num] = ful_polygon.shape[0]
        return ful_polygon

    except TopologicalError:
        pass


for i in range(shapefile.shape[0]):
    
    print(i)
    balance = 1
    polygon_df = splits[i]
    slope = polygon_df.loc[i, 'Slope']
    aspect_polygon = polygon_df.loc[i, 'Seg_Aspect']
    
    # Orthogonal projection of Length and width of PV panel accoding to the slope and arrangement (landscape or portrait)
    length1 = round(math.cos(math.radians(slope)), 2) * panel_length
    width1 = panel_width
    portrait = [length1, width1]

    width2 = round(math.cos(math.radians(slope)), 2) * panel_width
    length2 = panel_length
    landscape = [length2, width2]
    
    while balance != 0:
        
        poly2 = polygon_df.loc[i, 'geometry']
        poly2.crs = crs
        count = {}
        
        # Fishnet grid in portrait mode
        fishnet(poly2, count, aspect_polygon, portrait, slope, row_distance, 1)
        # Fishnet grid in landscape mode
        fishnet(poly2, count, aspect_polygon, landscape, slope, row_distance, 2)
        
        # Checks in which mode, more number of PV panels can be placed on rooftop segment
        try:
            highest_mode = max(count.items(), key=operator.itemgetter(1))[0]
            balance = max(count.items(), key=operator.itemgetter(1))[1]
            
        except:
            balance = 0
        
        # Checks if there is no space on the rooftop segment for a PV panel and breaks the iteration
        if balance == 0:
            break
            
        # Performs PV panel slicing according to the highest mode
        if highest_mode == 1:
            pvp_full = fishnet(poly2, count, aspect_polygon, portrait, slope, row_distance)
        else:
            pvp_full = fishnet(poly2, count, aspect_polygon, landscape, slope, row_distance)
            
        # Joins the PV panels together that belongs to a rooftop segment
        gdf = gpd.GeoDataFrame(pd.concat([gdf, pvp_full], ignore_index=True))
        
        # Checks for duplicate PV panel geometries and breaks the iteration
        if not gdf["geometry"].is_unique:
            break
    
        try:
            polygon_df = gpd.overlay(polygon_df, pvp_full, how='difference')
        except TopologicalError:
            break
            
        # Checks if rooftop polygon can accomodate anymore PV panels and breaks the iteration
        if polygon_df.empty or slope == flat_slope:
            break
