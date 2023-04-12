# Lagrangian_Particle_Tracking_SuppAnimations
Code and datasets for: Lilly et al. (2022). Using a Lagrangian particle tracking model to evaluate impacts of El Niño-related advection on euphausiids in the southern California Current System. Deep-Sea Research I.

***NOTE: Must contact for two additional large files ('ETOPO1_Ice_g_gmt4.grd', 'Western.shp').



% README: Laura E. Lilly (lauralilly4@gmail.com) - 22 Mar 2022


% Summary: 
Supplemental animations of particle backtracks for the distributions of six euphausiid species off Central and Southern California, 2008-2017

% Please cite as:
Lilly, Laura E.; Cornuelle, Bruce D.; Ohman, Mark D. (2022). Data from: Using a Lagrangian particle tracking model to evaluate impacts of El Niño-related advection on euphausiids in the southern California Current System. UC San Diego Library Digital Collections. https://doi.org/10.6075/J0KK9BZP

% Primary associated publication:
Lilly, L. E., Cornuelle, B. D., and Ohman, M., D. Submitted. Using a Lagrangian particle tracking model to evaluate impacts of El Niño-related anomalous advection on euphausiids in the southern California Current System. 

% Inquiries about the manuscript, animation files, or datasets can be sent to: 
Laura E. Lilly (lauralilly4@gmail.com)


% Description of contents:
The above supplemental animations depict particle backtracks from the spring distributions of six euphausiid species off Central and Southern California to the prior winter (four month backtrack). Euphuasiid distributions were developed from data from the CalCOFI program, which consists of in situ measurements of euphausiid biomass on a set grid of sampling stations one time per spring, between mid-March and mid-May. We first objectively mapped the biomass values from the individual CalCOFI stations onto a smooth, interpolated surface to depict the distribution across the entire region for a single year. From that surface, we identified the threshold contours encompassing the regions of >80% biomass and 50-80% biomass for that year's distribution. We then seeded particles within those contours and ran them back for four months (Mar 31 to the prior Dec 1) to determine the winter origins of waters impacting each population in spring. The animations capture those particle backtracks. The associated datafiles and Matlab scripts allow the user to recreate the supplemental animations.


% Methods:
Distributions for particle backtracks were seeded based on actual species distributions measured by the California Cooperative Oceanic Fisheries Investigations (CalCOFI) program in Apr-May of each year. Dark blue particles originated from the region containing biomass values that were >80% of the maximum biomass value of that year (e.g., if max biomass for a year = 1 mgC m-3, dark blue particles were seeded for the region encompassing all biomasses > 0.8 mgC m-3). Turquoise particles were seeded for regions of biomass between 50-80% of a year's maximum biomass. Particles were backtracked using the Runge-Kutta interpolation method, with 100 interpolation steps per 24-hr period. Backtracks were forced by modeled u and v flow fields from the California State Estimate (CASE). 


% List of Datafiles: 

- 'Matlab codes and datafiles (ParticleTrackCodes_Matlab.zip)' - Folder containing Matlab scripts and associated files necessary to re-create the particle backtrack animations.
	- 'Lillyetal_EuphPrtTrcks_SuppAnms_Step1.m' - Matlab script, Part 1 (run first)
	- 'Lillyetal_EuphPrtTrcks_SuppAnms_Step2.m' - Matlab script, Part 2 (run second)
	- 'Lillyetal_EuphPrtTrcks_SuppAnms_Step3.m' - Matlab script, Part 3 (run third)
	- 'Epaci_CSB_nt_tot_unt_bm.csv' - Example file of station-specific biomass for Euphausia pacifica (necessary to run Matlab scripts)
	- 'FMT.mat' - Matlab structure necessary to run scripts
	- 'grid.mat' - Matlab structure necessary to run scripts
	- 'MIT_CCS_LATLON_500m.mat' - Matlab structure necessary to run scripts
	- 'read_mitstate_CCS.m' - Matlab sub-code necessary to run scripts
	- 'TIME_CCS_2007_2020.bin' - Binary file of time-steps for U and V velocities. Necessary to run Matlab scripts
	- 'Coastline_files' (folder)
		- 'bc_entidad.dbf'
		- 'bc_entidad.prg'
		- 'bc_entidad.shp'
		- 'bc_entidad.shx'
		- 'bcs_entidad.dbf'
		- 'bcs_entidad.prj'
		- 'bcs_entidad.shp'
		- 'bcs_entidad.shx'
		- 'ETOPO1_Ice_g_gmt4.grd'
		- 'son_entidad.dbf'
		- 'son_entidad.prj'
		- 'son_entidad.shp'
		- 'son_entidad.shx'
		- 'Western.dbf'
		- 'Western.prj'
		- 'Western.shp'
		- 'Western.shx'
		
- 'U velocity file (U_MIT_CCS_2007_2020_fulldomain.bin)' - Binary file of U velocities (daily step) from Jan 1, 2007 - Mar 21, 2017. NOTE file size: 70 GB. File necessary to run Matlab scripts. Place in same directory with above files and scripts

- 'V velocity file (V_MIT_CCS_2007_2020_fulldomain.bin)' - Binary file of V velocities (daily step) from Jan 1, 2007 - Mar 21, 2017. NOTE file size: 70 GB. File necessary to run Matlab scripts. Place in same directory with above files and scripts

- 'Particle Backtrack Animations - [species name] - Backtrack animation files partitioned by the six euphausiid species (Euphausia pacifica, Nematoscelis difficilis, Thysanoessa spinifera, Nyctiphanes simplex, Euphausia eximia, Euphausia gibboides) and run for each of 10 springs (2008-2017) for which a given species was present. Forty-eight .avi files (10-15 MB each) total. Each file shows the particle backtrack from Mar 31 to the prior Dec 1 for a given year and euphausiid species. Not every species was present in every year. Individual files contain backtrack dates ('MarYYYY_DecYYYY').



% CODING STEPS
NOTE: You must separately download and format the file of euphausiid species biomass. A test file is provided for Euphausia pacifica; all other species must be requested.

1. Download Matlab scripts and datafiles
2. Place all datafiles in the same folder as scripts.
3. Email mohman@ucsd.edu to request access to biomass data for euphausiid species from the Brinton-Townsend Euphausiid Database (https://oceaninformatics.ucsd.edu/euphausiid/). Note that you can access the site directly as a public user, but you will only be able to download abundance data, not biomass. Request *untransformed carbon biomass* values for all available years and stations for the species you are interested in. 
4. From the dataset you receive, input the following variables into a .csv file, in order: RowNumber, Cruise, Mon, Day, Year, Line, Station, Region, Latitude, Longitude, Biomass
5. Add the .csv file to your datafiles folder, and in the Matlab script for 'Step2', Line 15, change the folder path and/or file name to your version of the file
6. Run the three files in order by 'Step' number. In Step2, you will be prompted to input 'Species name' and 'Input year'. The species name is to select which .csv file to run; the input year subsets the requested year for objective mapping and particle tracking. Input the species name in single quotes (') and the year without quotes. 


% Data dictionary:
- 'Epaci_CSB_nt_tot_unt_bm.csv':
	- 'RowNumber' - Consecutive number for each row in the original file, arranged by line and station; rows have been rearranged to order by date. (numeric value) 
	- 'Cruise' - Cruise identification number using CalCOFI labeling protocol ('CalCOFI_YYMM') or other cruise-specific labeling. (text value)
	- 'Mon' - Month of year (1-12). (numeric value)
	- 'Day' - Day of month (1-31). (numeric value)
	- 'Year' - Year (1951-2018). (numeric value)
	- 'Line' - CalCOFI sampling line (extending from coast to offshore) from which the sample was taken. See https://calcofi.org/sampling-info/station-positions/ for map and geographic positions of lines. (numeric value)
	- 'Station' - CalCOFI sampling station, along a CalCOFI line, from which the sample was taken. See link from 'Line' for station geographic positions. (numeric value)
	- 'Region' - CalCOFI-delineated region off California: 'CC' - Central California (San Francisco, CA, to Pt. Conception, CA); 'SC' - Southern California (Pt. Conception, CA, to U.S.-Mexico border); 'NBC' - Northern Baja California; 'CBC' - Central Baja California'; 'SBC' - Southern Baja California. See the following link for the map of historical sampling stations that included the three Baja regions: https://calcofi.org/about/history/. (text value)
	- 'Latitude' - Exact latitude of sampling location, in decimal degrees. Note that a line-station combination always has approximately the same latitude and longitude. (numeric value)
	- 'Longitude' - As for latitude above, except longitude of sampling location, in decimal degrees. (numeric value)
	- 'Biomass' - Biomass of species sampled at sampling location, in mg carbon / m2. (numeric value)

- 'FMT.mat' - Matlab structure with two character components to specify the subset slide of 'TIME_CCS_2007_2020.bin' file. Components are: 'fmt' = 'real*4'; 'Ieee' = 'b'.

- 'MIT_CCS_LATLON_500m.mat':
	- 'ZC' - List of all depth levels in CASE model domain. Depth increments in meters. 34 entries, ranging from -5 to -485, where negative value indicates meters below sea level. (numeric value). 
	- NOTE: Although this structure has 'XC' and 'YC' variables, those are NOT used in this script because they depict a smaller subset of the CASE model. The versions of those variables used here come from 'grid.mat' (below). 
	- NOTE: Although 'grid.mat' has a z dimension for the full U and V files, we use the 'ZC' list of depths from this file because it is a subset only down to 485 m, and we do not need the entire 5000m suite of depths. We also specify a subset range of depths for the actual particle interpolations.
	
- 'grid.mat':
	- 'hFacC' - variable that provides x, y, and z dimensions of the U and V flow fields. 	
	- 'XC' - List of all longitude values in CASE model domain. In decimal degrees on 0-360 degree range. 192 values, ranging from 232.0313 - 243.9688. (numeric  value)
	- 'YC' - List of all latitude values in CASE model domain. In decimal degrees. 128 values, ranging from 30.0438 - 37.9813. (numeric value)
	- NOTE: All other variables in this .mat structure are NOT used in this script. 


% Technical details:
MATLAB v.R2018a
