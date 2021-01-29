
clear all

%% Load overlap maps and break them up by state, then create uniqueID maps for the maps from each state
%
% This script loads overlap maps (see Make_State_Overlap_Maps.m for more
% documentation on overlap maps) and creates a uniqueID map for the parcels
% that belong to each distinct "state" on the map (task variant vertices
% only, rest variant vertices only, vertices identified in both task and
% rest states). First, the overlap maps are broken up by "state" and 
% binarized (1s indicate a variant is a present at a vertex, 0s indicate no
% variant is present). Then, these binarized maps are converted into
% uniqueID files.
%
% INPUTS:
%
% -outputdir: the output directory for the uniqueID maps for each "state"
% and subject (note that subdirectories are created within this directory 
% for the output files)
% -workbenchdir: the filepath to connectome workbench
% -leftsurf: the filepath and extension for the left hemisphere surface
% midthickness
% -rightsurf: the filepath and extension for the right hemisphere surface
% midthickness
%
% -threshold: a threshold of overlap maps to load and create uniqueID maps
% for
%
% -overlap_files: reads path to a text file containing the paths to each of
% the overlap map files. The format is (pathtofile subjectID) and the file 
% is space-delimited (see below for more formatting details)
%
% OUTPUTS:
%
% -Variant Overlap UniqueID Maps: A CIFTI file (numvertices x 1) with 
% a unique identifier for each contiguous parcel for each "state" (task 
% variant vertices only, rest variant vertices only, vertices identified in
% both task and rest states). The uniqueID map for each state is saved into
% a separate file for each subject and each "state" is saved into its own
% subdirectory
%
% Written by BK (01-2021)
%

%% Initialize Variables

outputdir = '/your/path/here/';  %% output directory for uniqueID files
workbenchdir = '/your/path/here/workbench/bin_macosx64/';  %% path to connectome workbench
leftsurf = '/your/path/here/32k_ConteAtlas_v2_distribute/Conte69.L.midthickness.32k_fs_LR.surf.gii';  %% filepath and extension for left hemisphere midthickness surface
rightsurf = '/your/path/here/32k_ConteAtlas_v2_distribute/Conte69.R.midthickness.32k_fs_LR.surf.gii';  %% filepath and extension for left hemisphere midthickness surface

threshold = 5;  %% threshold of maps to load

%% Set and check output directories

if ~exist([outputdir 'Variant_Overlap_Maps/'], 'dir')
    mkdir([outputdir 'Variant_Overlap_Maps/'])
end

if ~exist([outputdir 'Variant_Overlap_Maps/Overlap_Only/'], 'dir')
        mkdir([outputdir 'Variant_Overlap_Maps/Overlap_Only/'])
end
overlapdir = [outputdir 'Variant_Overlap_Maps/Overlap_Only/'];

if ~exist([outputdir 'Variant_Overlap_Maps/Rest_Only/'], 'dir')
        mkdir([outputdir 'Variant_Overlap_Maps/Rest_Only/'])
end
restdir = [outputdir 'Variant_Overlap_Maps/Rest_Only/'];

if ~exist([outputdir 'Variant_Overlap_Maps/Task_Only/'], 'dir')
        mkdir([outputdir 'Variant_Overlap_Maps/Task_Only/'])
end
taskdir = [outputdir 'Variant_Overlap_Maps/Task_Only/'];

% load filepaths and sub IDs
% load overlap maps for analyses (see Make_State_Overlap_Maps.m for
% additional documentation on overlap maps) using a space-delimited text 
% file in the format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01

[overlap_files, subjects, ~] = textread(['/Users/briankraus/Desktop/Final_Task_Rest_Files/Variant_Overlap_Maps/MSC_overlap_varmaps_' num2str(threshold) '.txt'],'%s%s%s');  %% text file containing filepaths for overlap maps


%% Loop through files, create a separate unique ID file for each type of vertex (overlapping, task only, rest only)


for x = 1:length(overlap_files)
    
    %% Load and binarize maps
    
    cifti_overlap_map = ft_read_cifti_mod(overlap_files{x});  % read each overlap map
    
    overlapvertsdat = zeros(size(cifti_overlap_map.data));  % create empty data vector
    restvertsdat = zeros(size(cifti_overlap_map.data));
    taskvertsdat = zeros(size(cifti_overlap_map.data));

    overlapverts = find(cifti_overlap_map.data == 1);  % find values in each map that correspond to each vertex type
    restverts = find(round(cifti_overlap_map.data,1) == .7);
    taskverts = find(cifti_overlap_map.data == .5);
    
    overlapvertsdat(overlapverts) = 1;  % binarize each data vector
    restvertsdat(restverts) = 1;
    taskvertsdat(taskverts) = 1;
    
    %% Write out spatial maps for each type of vertex
    
    % overlap data

    cifti_overlap_map.data = overlapvertsdat; % replace data for each vertex type
    
    [outfile, outfileuniqueID] = create_filename(subjects{x},threshold,'Overlaponly');  % create filenames for output files
    outfileoverlap = [overlapdir outfile];  % filepath for overlap file
    outfilewboverlap = [overlapdir outfileuniqueID];  % filepath for overlap unique ID file
    
   	ft_write_cifti_mod(outfileoverlap, cifti_overlap_map);  % write overlap output file
    system([workbenchdir 'wb_command -cifti-find-clusters ' outfileoverlap ' 0 0 0 0 COLUMN ' outfilewboverlap ' -left-surface ' leftsurf ' -right-surface ' rightsurf]);  % create overlap unique ID file
    delete(outfileoverlap);  % delete binarized overlap CIFTI
    
    % rest data
    
    cifti_overlap_map.data = restvertsdat; % replace data for each vertex type
    
    [outfile, outfileuniqueID] = create_filename(subjects{x},threshold,'Restonly');  % create filenames for output files
    outfilerest = [restdir outfile];  % filepath for overlap file
 	outfilewbrest = [restdir outfileuniqueID];  % filepath for overlap unique ID file
    
    ft_write_cifti_mod(outfilerest, cifti_overlap_map);  % write rest output file
    system([workbenchdir 'wb_command -cifti-find-clusters ' outfilerest ' 0 0 0 0 COLUMN ' outfilewbrest ' -left-surface ' leftsurf ' -right-surface ' rightsurf]);  % create rest unique ID file
    delete(outfilerest);  % delete binarized rest CIFTI
    
    % task data
    
    cifti_overlap_map.data = taskvertsdat; % replace data for each vertex type
  	
    [outfile, outfileuniqueID] = create_filename(subjects{x},threshold,'Taskonly');  % create filenames for output files
    outfiletask = [taskdir outfile];  % filepath for overlap file
    outfilewbtask = [taskdir outfileuniqueID];  % filepath for overlap unique ID file
    
  	ft_write_cifti_mod(outfiletask, cifti_overlap_map);  % write task output file
    system([workbenchdir 'wb_command -cifti-find-clusters ' outfiletask ' 0 0 0 0 COLUMN ' outfilewbtask ' -left-surface ' leftsurf ' -right-surface ' rightsurf]);  % create task unique ID file
    delete(outfiletask);  % delete binarized task CIFTI
    
end


function [outfile, outfileuniqueID] = create_filename(sub,threshold,comparison)

% set output filename for each map

task = 'All';
ses = 'All';
grp = 'WashU120';
desc = [comparison '-matched'];
thresh_str = num2str(threshold);

outfile = sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_SubunitMap.dtseries.nii',sub,task,ses,grp,desc,thresh_str);
outfileuniqueID = sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_UniqueIDMap.dtseries.nii',sub,task,ses,grp,desc,thresh_str);  

end






  
