
clear all

%% Loads in variants thresholded at a given value and counts the number of vertices
%
% This script loads in uniqueID files in the referenced text files and 
% counts the number of vertices in each spatial location map that contain a
% variant at a given threshold
%
% INPUTS:
%
% -minsize: the minimum size threshold (in vertices) that a variant has to
% meet in order to be included in the output data
% -SizeExclude: toggles whether to exclude variants below a minimum size
% (in vertices) threshold (set to 1), or use all variant locations
% regardless of size (set to 0)
% -threshold: an absolute (r) threshold to load for counting vertices
%
% -task_files: the filepath to the text file that contains a
% filepath for each task uniqueID map (see threshold_variant_maps.m for
% additional documentation on uniqueID maps) and it's associated subject ID
% (space delimited) on each line of the file
% (e.g.,filepath/MSC01.dtseries.nii MSC01)
% -rest_files: the filepath to the text file that contains a
% filepath for each rest uniqueID map (see threshold_variant_maps.m for
% additional documentation on uniqueID maps) and it's associated subject ID
% (space delimited) on each line of the file
% (e.g.,filepath/MSC01.dtseries.nii MSC01)
% the order of the data files for all of the subjects should be the
% same in all text files
%
% OUTPUTS:
%
% -vectors with the count of variant vertices for each subject and state,
% with each row corresponding to each subject
%
% Written by BK (01-2021)
%

%% Initialize Variables

SizeExclude = 1;  %% Toggles whether to exclude small variants
minsize = 50;  %% Minimum size of variants to include in count
threshold = 0.3;  %% Threshold (r) to use for counting vertices
       
% read text files to load data from

[task_files, subjects, ~] = textread(['/your/path/here/MSC_alltasks_varmaps_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
    
[rest_files, ~, ~] = textread(['/your/path/here/MSC_rest_varmaps_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');


%% Loop through subjects and count the number of vertices corresponding to variants in each map

subs = unique(subjects);

RestVerts = [];
TaskVerts = [];

for s = 1:numel(subs)  %% for each subject
    
    subject = subs{s};
    
    cifti_rest = ft_read_cifti_mod(rest_files{s});  %% load data
    cifti_task = ft_read_cifti_mod(task_files{s});


        if SizeExclude    %% Exclude variants below size threshold

            cifti_rest.data = variant_size_exclude(cifti_rest,minsize);
            cifti_task.data = variant_size_exclude(cifti_task,minsize);
            
        end
        
        for h = 1:length(cifti_rest.data)    %% Binarize maps
            
            if cifti_rest.data(h) > 0
                
                cifti_rest.data(h) = 1;
                
            end
        end
        
        for h = 1:length(cifti_task.data)    %% Binarize maps
            
            if cifti_task.data(h) > 0   %% Else binarize
                
                cifti_task.data(h) = 1;
                
            end
        end
    
    RestVerts = [RestVerts; sum(cifti_rest.data)];  %% add sum of binarized maps in each state to vector
    TaskVerts = [TaskVerts; sum(cifti_task.data)];
    
end
