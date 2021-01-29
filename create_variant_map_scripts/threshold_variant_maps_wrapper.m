
clear all

%% Sets up inputs and filepaths for thresholding variant maps
%
% This script sets up the inputs for the function that thresholds variant
% maps. Here, the settings are used to exclude vertices from consideration
% with low SNR and set the threshold(s) to use for defining variants.
%
%﻿OUTPUTS:
%
% -workbenchdir: the filepath to connectome workbench
% -surfdir: the filepath to the directory that has a path to a template for
% the surface of the brain
% -SNRpath: the filepath to the directory that contains the SNR CIFTIs for
% each subject. These files were created by mode 1000 normalizing the
% resting state data for each session within each subject. The mean BOLD
% value for each voxel in each session was calculated, and then this SNR
% value for each voxel was averaged across sessions. This data was then
% transformed to a CIFTI (numvertices x 1) and used for excluding low SNR
% regions (mean signal < 750 after mode 1000 normalization)
% -spatialmapslist: the filepath to the text file that contains a
% filepaths for each variant map and it's associated subject ID (space
% delimited) on each line of the file
% (e.g.,filepath/MSC01.dtseries.nii MSC01)
% -outputdir: the output directory for the thresholded variant files (note
% that subdirectories are created within this directory for the output
% files)
%
% -thresholds: a vector (or single value) of thresholds to define variants
% for
% -SNRexclusion: toggles whether to exclude low signal regions from
% consideration as variants (set to 1), otherwise allow all vertices to be
% defined as variants (set to 0)
% -absolutethresholds: toggles whether to use an absolute correlation (r)
% value to define variants (set to 1), otherwise thresholds are interpreted
% as the xth lowest percentage of correlation values (set to 0)
%
% Written by BK (01-2021)
%

workbenchdir = '/your/path/here/workbench/bin_macosx64/';  %% Connectome workbench filepath
surfdir = '/your/path/here/32k_ConteAtlas_v2_distribute/';  %% path to directory containing left and right hemisphere surfaces
SNRpath = '/your/path/here/';  %% path to SNR maps
spatialmapslist = '/your/path/here/Variant_Maps/Variant_Maps_SubjectList.txt';  %% path of text file with subject/file names
outputdir = '/your/path/here/';  %% top output directory for thresholded maps

thresholds = 5;  %% Thresholds used to calculate variants (lowest % of correlation values unless absolutethresholds set to 1)
SNRexclusion = 1;  %% Toggles whether to exclude vertices based on SNR (set to 1 to exclude)
absolutethresholds = 0;  %% Toggles whether to use absolute or percent thresholds to define variants  (set to 1 for absolute)

threshold_variant_maps(thresholds,SNRexclusion,absolutethresholds,outputdir,spatialmapslist,workbenchdir,surfdir,SNRpath);