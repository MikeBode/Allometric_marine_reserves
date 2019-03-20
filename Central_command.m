
% This function runs all functions and creates all figures and output spreadsheets.
% Before you run it, make sure you have a subfolder called "Figures" in the directory

clear all

% Step 1: Extract the species parameters from an XLSX spreadsheet and save as a Matlab .MAT file
Extract_species_data

% Step 2: For each of these species, calculate the Beverton-Holt parameter values that yield a particular
%         "steepness" value, defined as the fraction of virgin recruitment that is attained when the 
%         spawning biomass is 20% of virgin spawning biomass.
Steepness_search

% Step 3: Now, for each level of the steepness parameter, calculate what the optimal MPA proportion is
%         (with and without harvests outside the reserves). This function will greate the outputs 
Compare_MPAs_Traditional

