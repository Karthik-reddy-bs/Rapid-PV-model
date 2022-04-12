# Rapid forecasting of Solar Photovoltaic energy production across cities

This is a repository which contains source code and documentation.

## Contents

[Introduction](https://github.com/Karthik-reddy-bs/Rapid-PV-model/edit/main/README.md#introduction)

[PV Models](https://github.com/Karthik-reddy-bs/Rapid-PV-model/edit/main/README.md#pv-models)

[Berlin Datasets](https://github.com/Karthik-reddy-bs/Rapid-PV-model/edit/main/README.md#berlin-datasets)

[Master Branch](https://github.com/Karthik-reddy-bs/Rapid-PV-model/edit/main/README.md#master-branch)

[Workflow](https://github.com/Karthik-reddy-bs/Rapid-PV-model/edit/main/README.md#workflow)

[Input parameters](https://github.com/Karthik-reddy-bs/Rapid-PV-model/edit/main/README.md#input-parameters)

## Introduction

PV models for rapid and continuous forecasting of the potential and production of photovoltaic energy at the neighbourhood and city scales. A combination of LiDAR (Light Detection and Ranging) geospatial data and meteorological forecasting data (SUEWS-UMEP) will achieve this.

### PV Models
**Model BR** - A building resolving model with detailed, explicit 3-dimensional PV panel placement considering rooftop geometry and topography to optimize each PV panel in an urban setting is considered as reference model. This is combined with a time-resolving urban radiation scheme (UMEP-SEBE) to calculate solar irradiance at 1 m2 resolution rooftop area. As this is computationally intensive, this model is used to: (i) parametrise other simpler PV models, and (ii) create a dataset to evaluate the simpler models.

**Model 1P** - Single panel model assumes one flat PV panel per grid cell, not impacted by geometry (e.g., shading or slope), only meteorology (ambient air temperature, incoming radiation fluxes). The neighbourhood scale output uses the roof area that can be covered by panels and the net losses (e.g. inefficient geometric placement and shading losses from objects such as trees, and buildings).

**Model 4P** - Four panel model assumes there are four tilted panels per grid cell, oriented in the four cardinal directions (i.e.  North 0°, East 90°, South, 180° and West 270°) with a optimized slope. On flat roofs the panel is assumed to be South facing (in our northern hemisphere application) but tilted, whereas on pitched roofs the PV power output is the total from the four panels. The grid value is the sum of the weighted  projected roof area of flat and pitched roofs, after accounting for the net losses (e.g. placement geometry, shading  from objects, trees and buildings and self-shading).

**Model MP** - Multi-panel model requires information on each roof segment’s - (i) area, (ii) slope and (iii) aspect. For flat roofs the same assumption as 4P flat roofs are made. For each pitched roof, a single panel is assumed with the actual slope and aspect of the rooftop segment. The power output per segment is calculated and aggregated to the roof area, but without losses (placement and shading). The grid cell power output is the total roof area with losses (placement and shading from surrounding objects, trees, buildings and self-shading) accounted for. 

## Berlin Datasets

For model training and evaluation, building morphology data for the city of Berlin are used, which are based on [vector data](https://www.businesslocationcenter.de/berlin3d-downloadportal/#/export) and [2 m resolution LiDAR data](https://www.berlin.de/umweltatlas/en/land-use/building-and-vegetation-heights/). Survey data is available as 2D LOD2 data containing information on simplified roof and wall locations. Within the LOD2 dataset, every building rooftop is segmented into rooftop segments (uniformly sloped roof parts, but not resolving obstructions such as chimneys, windows etc.). Rooftop segment is a reasonably large part of a roof with the same slope and same aspect. Each building has at least one segment but can have multiple rooftop segments (e.g., pitched roof), but segments cannot span several buildings even if slope and aspect of neighbouring buildings are similar. Along with buildings, the [vegetation heights](https://www.berlin.de/umweltatlas/en/land-use/building-and-vegetation-heights/) also considered to include the effect of vegetation on shading and PV energy.

The input parameters for the reference BR are the building and vegetation DSMs, the segmented LOD2 vector dataset that delineate roof segments, and meteorological data. For Models 1P, 4P and MP, in each grid cell, integral urban form parameters are derived which are Coordinates of grid cell centre (Lon, Lat), number of rooftop segments per ha, average building height, standard deviation of building height, plan area fraction of buildings, fraction of flat segments, and average segment length.  Additionally, meteorological data (Incident solar direct and diffuse radiation, ambient temperature, solar elevation and zenith angles), which changes over time but is considered spatially homogeneous across the grid cell is used as input. When integrated into a mesoscale model, the weather can be different for different grid cells but is uniform for a single grid cell. For Model MP, additional information is needed, such as size, slope, and aspect of every roof segment in the grid cell, which needs to be preprocessed once at model set-up. 

## Master Branch

Here contain two folders in which under the folder Datasets contains sample simulation files () which can be used to run the PV models for a certain urban area in Berlin. And the other folder python files contain the python scripts for all PV models.

## Workflow

	1. Download the repository
	2. Select an PV model by opening the respective python script in the folder python files
	3. Go through the description of the model in the script which states about the model, how it is developed and what are its advantages and limitations.
	4. The Input parameters in the script varies according to the PV model and some of these input parameters need to be estimated beforehand.
	5. Once the input parameters are set, the script can be executed and the output of the file will be a data frame which contains hourly PV power output [MW km-2] for the simulated area.
  
## Input Parameters

	1. N: grid cell total number of rooftop segments
	2. A_G: Total plan area of grid cell [m2]
	3. A_R: Projected plan area of all rooftop segments in a grid cell [m2]
	4. λ_p: Projected plan area fraction of building segments (A_R/A_G). It is proposed as a scaling factor for the area of roofs present in each grid cell that can be equipped with PV panels. As no overlapping roofs are considered, λ_p ranges between 0 (no buildings in grid cell present) and 1. Values of λ_p are commonly present as inputs in mesoscale meteorological models, often identified as plan area of buildings for urban surface parameterization schemes. This information can be extracted from aerial imagery, LIDAR or LOD1 or 2 datasets.
	5. l_s: Average building segment length [m] is the typical distance of the length of the roof, which tells you how far the shadows get on the roof and is calculated by taking the square root of the average segment area ((l_s = [A_R / N] 0.5 ))
	6. f_f: Fraction of projected rooftop segment areas that are flat (0-1) is used to parameterize placement losses and consider orientation effects of panels on flat roofs.
	7. f_p: Scaling factor coefficient (f_p = 1 if all the rooftop segments are covered with panels)
	8. α_s: System albedo for an urban area(0-1)
	9. l: Length of PV panel in portrait mode [m]
	10. w: Width of PV panel in portrait mode [m]
	11. d: Panel row spacing, can be optimized according to technical & Economical parameters.
	12. β_o: Optimal slope of PV panels on flat rooftops, can be optimized according to location. In this script, β_0 = 40° has been used for the city of Berlin.
	13. η_p: Efficiency of PV panel [%]
	14. α_T: Temperature coefficient of PV panel [1/K]
	15. η_W: Wiring losses
	16. metadata: Input meteorological data specifically formatted using 
	      UMEP -> Pre-processing -> Meteorological data -> Prepare existing data
	      Here, a sample Random Meteorological Day (RMD) for the month of June has been given as an example. Annual weather data for the city of Berlin has also been attached for reference.
	17. sol_angles: Sun position angles for the urban area. Here, a sample set of solar angles for the RMD - June has been given as example


