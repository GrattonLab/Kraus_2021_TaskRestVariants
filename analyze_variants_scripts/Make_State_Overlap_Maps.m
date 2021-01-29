
clear all

%% Create spatial maps (CIFTIs) showing both spatial location overlap of variants across states and the similarity of variant spatial correlation values across states
%
% This script is for creating maps showing topographic differences in
% variant spatial location and correlation magnitude across states. For 
% creating magnitude maps, a variant map (correlation of individual to 
% group) is necessary. For applying a size exclusion and creating a map of 
% spatial location overlap, a uniqueID map is also necessary. This input 
% data can be excluded by size or by SNR mask. These input maps are used to
% create a spatial location map showing vertices where variants occur 
% either only during tasks, only during rest, or in both states. For the 
% correlation magnitude maps, the correlation values for the vertices that 
% correspond to variants are subtracted across states (e.g. the correlation
% value for a vertex corresponding to a variant in task is subtracted from 
% the correlation value of the corresponding vertex in rest). The absolute 
% value of this subtraction is then calculated and these values are placed 
% on the surface at the corresponding vertices. These maps are then written
% out as CIFTIs in their own respective subdirectories.
%
% INPUTS:
%
% -outputdir: the output directory for the spatial location maps for each
% subject (note that subdirectories are created within this directory for 
% the output files)
% -SNRpath: the filepath to the directory that contains the SNR CIFTIs for
% each subject (see see threshold_variant_maps.m for additional
% documentation)
% -SNR_fstring: the filename and extension for the SNR map
% -templatepath: the path to a template CIFTI to use for outputting spatial
% locations of each bin on the surface
%
% -SNRExclude: toggles whether to exclude low signal regions from
% consideration as variants (set to 1), otherwise allow all vertices to be
% defined as variants (set to 0)
% -minsize: the minimum size threshold (in vertices) that a variant has to
% meet in order to be included in this analysis
% -SizeExclude: toggles whether to exclude variants below a minimum size
% (in vertices) threshold (set to 1), or use all variant locations
% regardless of size (set to 0)
% -threshold: a threshold of uniqueID maps to load and create spatial
% location maps for
%
% -task_variant_maps: reads path to a text file containing the paths to 
% each of the variant files for task data for all sessions. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
% -rest_variant_maps: reads path to a text file containing the paths to 
% each of the variant files for rest data for all sessions. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
% -task_UniqueIDs: reads path to a text file containing the paths to each 
% of the uniqueID files for task data for all sessions. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
% -rest_UniqueIDs: reads path to a text file containing the paths to each 
% of the uniqueID files for rest data for all sessions. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
%
%
% OUTPUTS:
%
% -Variant Overlap Maps: A CIFTI file (numvertices x 1) with values that
% correspond to whether a given vertex is identified as a variant in both
% states (value = 1), identified as a variant in a rest state only (value 
% = .7), identified as a variant in a task state only (value = .5), or not
% identified as a variant in either state (value = 0)
% -Variant Magnitude Maps: A CIFTI file (numvertices x 1) with values that
% correspond to the correlation magnitude differences between states. If a
% vertex is identified as a variant in either task or rest states, its
% value is calculated by the absolute value of the correlation difference
% between states (i.e., abs(taskdatacorrelation-restdatacorrelation)). The
% non-variant vertices are first set to 0, and so have no value after this
% operation. Thus, any non-zero value in this map represents the absolute 
% difference in correlation values across states for vertices identified as
% variants (in either state)
%
% Written by BK (01-2021)
%


%% Set variables for script

templatepath = '/your/path/here/template_cifti.dtseries.nii';  %% Path where template cifti is located
SNRpath = '/your/path/here/';  %% path to SNR maps
SNR_fstring = 'SNRname.dtseries.nii';  %% filename and extension for SNR map
outputdir = '/your/path/here/';  %% output for CIFTI files

SNRExclude = 1;  %% Toggles whether to remove low SNR regions from consideration (set to 1)
SizeExclude = 1;  %% Toggles whether to remove variants below a certain size from consideration (set to 1)
minsize = 50;     %% Minimum size for exclusion of variants (vertices)
threshold = 5;  %% Threshold of variants to create maps from


% Create output directories and set output paths

if ~exist([outputdir 'Variant_Magnitude_Maps'], 'dir')
    mkdir([outputdir 'Variant_Magnitude_Maps']);
end
outputdirmag = [outputdir 'Variant_Magnitude_Maps/'];

if ~exist([outputdir 'Variant_Overlap_Maps'], 'dir')
    mkdir([outputdir 'Variant_Overlap_Maps']);
end
outputdirovrlp = [outputdir 'Variant_Overlap_Maps/'];


% Variant maps
% load variant maps for analyses using a space-delimited text file in the 
% format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% the order of the data files for all of the subjects should be the same in
% all text files

[task_variant_maps, subjects, ~] = textread('/your/path/here/MSC_alltask_varmaps.txt','%s%s%s');

[rest_variant_maps, ~, ~] = textread('/your/path/here/MSC_rest_varmaps.txt','%s%s%s');

% Unique ID maps
% load uniqueID maps for analyses (see threshold_variant_maps.m for
% additional documentation on uniqueID maps) using a space-delimited text 
% file in the format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% the order of the data files for all of the subjects should be the same in
% all text files

[task_UniqueIDs, ~, ~] = textread(['/your/path/here/MSC_alltasks_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');

[rest_UniqueIDs, ~, ~] = textread(['/your/path/here/MSC_rest_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');



%% Loop through subjects, get spatial overlap and correlation magnitude difference across states


subs = subjects;   % Get subject IDs

for s = 1:numel(subs)   % Loop through subjects
    
    subject = subs{s};  % read subject ID from input text file
    
    task_variant_file = task_variant_maps{s};  % loads variant maps for finding magnitude difference
    rest_variant_file = rest_variant_maps{s};
    
    variant_task = ft_read_cifti_mod(task_variant_file);  % load task and rest file for each subject
    variant_rest = ft_read_cifti_mod(rest_variant_file);


  	task_uniqueID_file = task_UniqueIDs{s};  % loads unique ID maps for spatial overlap
 	rest_uniqueID_file = rest_UniqueIDs{s};
        
	task_uniqueid = ft_read_cifti_mod(task_uniqueID_file);    % load unique ID files
 	rest_uniqueid = ft_read_cifti_mod(rest_uniqueID_file);
    
    if SNRExclude     % Performs SNR Exclusion
        
        SNRmap = ft_read_cifti_mod([SNRpath subject SNR_fstring]);  % read in SNR map cifti
        SNRmap.data = SNRmap.data(1:59412,:);   % Select cortical SNR only
        LowSNR = find(SNRmap.data < 750);      % Find vertices with SNR less than 750
        
        variant_task.data(LowSNR) = 0;  %% set low SNR vertices in variant map to 0
        variant_rest.data(LowSNR) = 0;
        
    end
    
    if SizeExclude   %% exclude variants smaller than size threshold
        
        task_uniqueid.data = variant_size_exclude(task_uniqueid,minsize);  % exclude variants by size from rest and task
        rest_uniqueid.data = variant_size_exclude(rest_uniqueid,minsize);
        
        task_uniqueid.data(task_uniqueid.data>0) = 1;  %% binarize data
        rest_uniqueid.data(rest_uniqueid.data>0) = 1;

        nonvariantlocations = find(~task_uniqueid.data & ~rest_uniqueid.data);  % find locations where variants do not occur in either task or rest
        
        variant_task.data(nonvariantlocations) = 0;  % set these values to 0 so that they have no value after subtraction
        variant_rest.data(nonvariantlocations) = 0;

    end
    
    %% Find overlap and magnitude differences
    
    % Find overlap between maps and set thresholds
    
    overlaymapcomp = mean([task_uniqueid.data'; rest_uniqueid.data'],1)';  %% conatenate data files and take mean
    
    overlap = find(overlaymapcomp==1);  %% in vertices where mean equals 1 (a variant occurred in both states), mark as overlapping between states
    
    restunique = setxor(find(rest_uniqueid.data==1), overlap);  %% find vertices that only occur in rest and not in both states
    
    taskunique = setxor(find(task_uniqueid.data==1), overlap);  %% find vertices that only occur in task and not in both states
    
    combinedoutput = zeros(size(rest_uniqueid.data));  %% set different thresholds for each overlap type in output map
    combinedoutput(overlap,:) = 1;
    combinedoutput(restunique,:) = .7;
    combinedoutput(taskunique,:) = .5;
    
    % get mean absolute difference across all for magnitude maps
    
    magnitude = abs(variant_rest.data - variant_task.data);

    
    %% Write out maps
    
    template = ft_read_cifti_mod(templatepath);     %% Read in template file for saving cifti
    template.data = [];
    
    magnitude_fname = create_filename(subject,threshold,'mag');  % write out magnitude map
    
 	template.data = magnitude;
        
   	ft_write_cifti_mod([outputdirmag magnitude_fname],template)
        
    overlap_fname = create_filename(subject,threshold,'ovrlp');  % write out overlap map
    
 	template.data = combinedoutput;
        
   	ft_write_cifti_mod([outputdirovrlp overlap_fname],template)
    
    clear template
    
end

    

function filename = create_filename(sub,threshold,map)

% set output filename for each map

task = 'All';
ses = 'All';
grp = 'WashU120';
desc = 'matched';
thresh_str = num2str(threshold);

if strcmp(map,'mag')

    filename = sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_MagnitudeMap.dtseries.nii',sub,task,ses,grp,desc,thresh_str);
    
elseif strcmp(map,'ovrlp')
    
    filename = sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_OverlapMap.dtseries.nii',sub,task,ses,grp,desc,thresh_str);

end
end
       