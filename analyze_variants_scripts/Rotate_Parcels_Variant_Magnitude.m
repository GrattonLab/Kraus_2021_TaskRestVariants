
clear all

%% Randomly rotate parcels (variants) within each hemisphere for each subject and calculate absolute difference in correlation values across states for each rotation
%
% This script takes the observed variants for each subject and finds the
% absolute magnitude difference (see Make_State_Overlap_Maps.m for more
% documentation on this calculation) between the correlations observed in
% the rest and task state. This is done both for vertices corresponding to
% variants in rest and task. Then, the uniqueID files from each state are 
% randomly rotated within the same hemisphere (by creating GIFTIs from each
% CIFTI). If after a rotation any individual variant overlaps with the 
% medial wall or another variant, the variant is randomly rotated until it 
% overlaps with neither of these. After performing the rotations, the 
% GIFTIs are recombined into a CIFTI and the absolute magnitude difference 
% in the correlations between states is calculated using the variant maps 
% and vertices selected in each rotation. These values are then compared to
% the observed variant correlation magnitude differences and plotted for 
% each subject as well as aggregated across subjects.
%
% INPUTS:
%
% -outputdir: the output directory for the summary plot and the output 
% spatial and rotation maps for each state and subject (note that 
% subdirectories are created within this directory for the output files)
% -giftitemplatepath: the filepath and filename to the GIFTI used as a
% template for converting GIFTIs to CIFTIs
% -workbenchdir: the filepath to connectome workbench
% -rotationfilepath: the filepath and filename for a file with X number of
% generated rotations for the x, y, and z planes
% -templatepath: the path to a template CIFTI to use for outputting spatial
% locations of each bin on the surface
% -neighborfilepath: path to node neighbor file generated from
% caret '-surface-topology-neighbors' (closest "neighbors" to each vertex)
% -atlasdir: the path to the directory that contains the 32k surface medial
% wall GIFTI files
% -surfdir: the path to the directory that contains the 32k surface medial
% wall GIFTI files
%
% -minsize: the minimum size threshold (in vertices) that a variant has to
% meet in order to be included in this analysis
% -SizeExclude: toggles whether to exclude variants below a minimum size
% (in vertices) threshold (set to 1), or use all variant locations
% regardless of size (set to 0)
% -threshold: a threshold of overlap maps to load and create rotations for
% -iterations: the number of rotations to perform within each hemisphere
% for each subject and state
% -generaterotations: toggles whether to generate new random rotation
% coordinates for each subject (set to 1), or to use the same pre-generated
% rotation coordinates for all subjects
% 
% -task_varmaps: reads path to a text file containing the paths to 
% each of the variant files for task data for all sessions. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
% -rest_varmaps: reads path to a text file containing the paths to 
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
% OUTPUTS:
%
% -CIFTI files: CIFTI spatial maps (numvertices x 1) which contain the 
% uniqueID maps with size excluded variants for each state and each 
% subject. The CIFTI files are then used to create the GIFTI files
% -GIFTI files: GIFTI spatial maps (numvertices x 1) which have been 
% separated into left and right hemisphere files using the CIFTI files 
% output. The GIFTI files are used to create the GIFTI rotation files
% -GIFTI rotation files: GIFTI spatial maps (numvertices x iterations) for
% which each timepoint in the file represents a spatial rotation. The GIFTI
% rotation files are used to make the CIFTI rotation files
% -CIFTI rotation files: CIFTI spatial maps (numvertices x iterations) for
% which each timepoint in the file represents a spatial rotation. These are
% created by recombining the GIFTI rotations into one CIFTI file. These
% files are used to create the final summary plot
% -plots: creates a plot for the magnitude differences in the correlation
% values between states for variant vertices vs. randomly rotated vertices.
% One plot is made for variants in each of the task and rest states. The
% results are plotted by subject as well as aggregated across subjects
%
% Written by BK (01-2021)
%

%% Initialize Variables

outputdir = '/your/path/here/';  %% output directory for plots and spatial maps
giftitemplatepath = '/your/path/here/gifti_template.func.gii'; %% path to GIFTI template for converting to CIFTIs
workbenchdir = '/your/path/here/workbench/bin_macosx64/';  %% Connectome workbench filepath
rotationfilepath = '/your/path/here/rotations1000.mat';  %% path to file containing pre-generated rotations
templatepath = '/your/path/here/template_cifti.dtseries.nii';  %% Path where template cifti is located
neighborfilepath = '/Users/briankraus/Documents/MATLAB/Open Source Functions/node_neighbors.txt';  %% path to nearest neighbors file (used in rotation function)
atlasdir = '/your/path/here/';  %% path to 32k atlas directory with medial wall file
surfdir = '/your/path/here/';  %% path to 32k atlas directory with medial wall file

SizeExclude = 1;    %% Toggles whether to exclude small variants
minsize = 50;       %% sets minimum size for variants (in vertices)
threshold = 5;    %% Threshold of variants to use for maps
iterations = 1000;   %% number of times to rotate variants for each subject
generaterotations = 1;  %% Toggles whether to generate new random rotations for each subject

%% Set and check output directories

if ~exist([outputdir 'Variant_Rotation_Maps/'], 'dir')
    mkdir([outputdir 'Variant_Rotation_Maps/'])
end

if ~exist([outputdir 'Variant_Rotation_Maps/GIFTIs/'], 'dir')
        mkdir([outputdir 'Variant_Rotation_Maps/GIFTIs/'])
end
giftidir = [outputdir 'Variant_Rotation_Maps/GIFTIs/'];

if ~exist([outputdir 'Variant_Rotation_Maps/GIFTI_Rotations/'], 'dir')
        mkdir([outputdir 'Variant_Rotation_Maps/GIFTI_Rotations/'])
end
giftirotdir = [outputdir 'Variant_Rotation_Maps/GIFTI_Rotations/'];

if ~exist([outputdir 'Variant_Rotation_Maps/CIFTIs/'], 'dir')
        mkdir([outputdir 'Variant_Rotation_Maps/CIFTIs/'])
end
ciftidir = [outputdir 'Variant_Rotation_Maps/CIFTIs/'];

if ~exist([outputdir 'Variant_Rotation_Maps/CIFTI_Rotations/'], 'dir')
        mkdir([outputdir 'Variant_Rotation_Maps/CIFTI_Rotations/'])
end
ciftirotdir = [outputdir 'Variant_Rotation_Maps/CIFTI_Rotations/'];

%% Loop through subjects and rotate variants within each hemisphere for each subject

% Unique ID maps
% load uniqueID maps for analyses (see threshold_variant_maps.m for
% additional documentation on uniqueID maps) using a space-delimited text 
% file in the format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% the order of the data files for all of the subjects should be the same in
% all text files

[task_UniqueIDs, subjects, ~] = textread(['/your/path/here/MSC_alltasks_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');

[rest_UniqueIDs, ~, ~] = textread(['/your/path/here/MSC_rest_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');

% Variant maps
% load variant maps for analyses using a space-delimited text file in the 
% format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% the order of the data files for all of the subjects should be the same in
% all text files

[task_varmaps, ~, ~] = textread('/your/path/here/MSC_alltask_varmaps.txt','%s%s%s');

[rest_varmaps, ~, ~] = textread('/your/path/here/MSC_rest_varmaps.txt','%s%s%s');


nfiles = length(rest_UniqueIDs);  %% total number of CIFTI files in text files

TrueRestMags = [];      %% Magnitude of actual variant correlation differences between states
TrueTaskMags = [];

RestRotateMags = [];    %% Magnitude of rotated variant correlation differences between states
TaskRotateMags = [];

for x = 1:nfiles  %% for each text file
    
    subject = subjects{x};  %% get the associated subject number
    
    rotations = get_rotations(rotationfilepath,generaterotations,iterations);  %% get rotations for each subject
            
    disp(['Performing parcel rotation for subject ' subject])
    
    RestRotateMagsSub = [];    %% Magnitude of rotated variant distances for subject
    TaskRotateMagsSub = [];

    cifti_rest_ID = ft_read_cifti_mod(rest_UniqueIDs{x});
    cifti_task_ID = ft_read_cifti_mod(task_UniqueIDs{x});
    
  	if SizeExclude  %% Remove variants smaller than a given threshold
        
        cifti_rest_ID.data = variant_size_exclude(cifti_rest_ID,minsize);
        cifti_task_ID.data = variant_size_exclude(cifti_task_ID,minsize);

    end
    
    % Save files as GIFTIs
    
   	template = ft_read_cifti_mod(templatepath);  %% load template CIFTI file
  	template.data = [];  %% remove data from template
  	template2 = template;
    
    template.data = cifti_rest_ID.data;  %% replace template data with current subject task and rest data
    template2.data = cifti_task_ID.data;
    
    ciftioutfilerest = create_filename(subject,'Rest',threshold,[],0,0);  %% write out CIFTIs, clear template
    ft_write_cifti_mod([ciftidir ciftioutfilerest],template);
    clear template
    
    ciftioutfiletask = create_filename(subject,'Task',threshold,[],0,0);  %% write out CIFTIs, clear template
    ft_write_cifti_mod([ciftidir ciftioutfiletask],template2);
    clear template2
    
    % break up CIFTIs into GIFTIs to rotate their variants within each hemisphere
    
    giftioutfilelhrest = create_filename(subject,'Rest',threshold,'LeftHemisphere',1,0);  %% write out GIFTIs for each state and hemisphere
    
    system([workbenchdir 'wb_command -cifti-separate ' ciftidir ciftioutfilerest ' COLUMN -metric CORTEX_LEFT ' giftidir giftioutfilelhrest]);
    
    giftioutfilerhrest = create_filename(subject,'Rest',threshold,'RightHemisphere',1,0);  %% write out GIFTIs for each state and hemisphere
    
    system([workbenchdir 'wb_command -cifti-separate ' ciftidir ciftioutfilerest ' COLUMN -metric CORTEX_RIGHT ' giftidir giftioutfilerhrest]);
    
    giftioutfilelhtask = create_filename(subject,'Task',threshold,'LeftHemisphere',1,0);  %% write out GIFTIs for each state and hemisphere
    
    system([workbenchdir 'wb_command -cifti-separate ' ciftidir ciftioutfiletask ' COLUMN -metric CORTEX_LEFT ' giftidir giftioutfilelhtask]);
    
    giftioutfilerhtask = create_filename(subject,'Task',threshold,'RightHemisphere',1,0);  %% write out GIFTIs for each state and hemisphere
    
    system([workbenchdir 'wb_command -cifti-separate ' ciftidir ciftioutfiletask ' COLUMN -metric CORTEX_RIGHT ' giftidir giftioutfilerhtask]);

    
    % Run rotation script with GIFTI files
    
    gifti_restLH = gifti([giftidir giftioutfilelhrest]);  %% load GIFTI for each state and hemisphere
    giftirotoutfilelhrest = create_filename(subject,'Rest',threshold,'LeftHemisphere',1,1);  %% name output GIFTI rotation file
    
    vertsLHrest = generate_rotated_parcellation(gifti_restLH.cdata,iterations,rotations,'L',[giftirotdir giftirotoutfilelhrest],neighborfilepath,atlasdir,surfdir);  %% run rotations for each state and hemisphere
    clear gifti_restLH
    
    gifti_restRH = gifti([giftidir giftioutfilerhrest]);  %% load GIFTI for each state and hemisphere
    giftirotoutfilerhrest = create_filename(subject,'Rest',threshold,'RightHemisphere',1,1);  %% name output GIFTI rotation file
    
    vertsRHrest = generate_rotated_parcellation(gifti_restRH.cdata,iterations,rotations,'R',[giftirotdir giftirotoutfilerhrest],neighborfilepath,atlasdir,surfdir);  %% run rotations for each state and hemisphere
    clear gifti_restRH
    
    gifti_taskLH = gifti([giftidir giftioutfilelhtask]);  %% load GIFTI for each state and hemisphere
    giftirotoutfilelhtask = create_filename(subject,'Task',threshold,'LeftHemisphere',1,1);  %% name output GIFTI rotation file
    
    vertsLHtask = generate_rotated_parcellation(gifti_taskLH.cdata,iterations,rotations,'L',[giftirotdir giftirotoutfilelhtask],neighborfilepath,atlasdir,surfdir);  %% run rotations for each state and hemisphere
    clear gifti_restLH
    
    gifti_taskRH = gifti([giftidir giftioutfilerhtask]);  %% load GIFTI for each state and hemisphere
   	giftirotoutfilerhtask = create_filename(subject,'Task',threshold,'RightHemisphere',1,1);  %% name output GIFTI rotation file
    
    vertsRHtask = generate_rotated_parcellation(gifti_taskRH.cdata,iterations,rotations,'R',[giftirotdir giftirotoutfilerhtask],neighborfilepath,atlasdir,surfdir);  %% run rotations for each state and hemisphere
    clear gifti_restRH
    
    
    % Convert output GIFTIs into CIFTIs

    ciftirotoutfilerest = create_filename(subject,'Rest',threshold,[],0,1);  %% write out CIFTIs
    
    attribute_gifti_data_to_cifti_MSC([giftirotdir giftirotoutfilelhrest], [giftirotdir giftirotoutfilerhrest], [ciftirotdir ciftirotoutfilerest], atlasdir, workbenchdir, giftitemplatepath);  %% convert GIFTIs to CIFTI
    
    
    ciftirotoutfiletask = create_filename(subject,'Task',threshold,[],0,1);  %% write out CIFTIs
    
    attribute_gifti_data_to_cifti_MSC([giftirotdir giftirotoutfilelhtask], [giftirotdir giftirotoutfilerhtask], [ciftirotdir ciftirotoutfiletask], atlasdir, workbenchdir, giftitemplatepath);  %% convert GIFTIs to CIFTI
        
end

%% Loop through rotations and calculate magnitude differences between states for each subject

for z = 1:nfiles
    
    RestRotateMagsSub = [];    %% Magnitude of rotated variant distances for subject
    TaskRotateMagsSub = [];
    
    subject = subjects{z};
    
    ciftirotoutfilerest = create_filename(subject,'Rest',threshold,[],0,1);  %% get rest CIFTI rotation map name
    ciftirotoutfiletask = create_filename(subject,'Task',threshold,[],0,1);  %% get task CIFTI rotation map name
    
    % Calculate magnitude difference between task and rest
    
    cifti_rest_ID = ft_read_cifti_mod(rest_UniqueIDs{z});  % load unique ID maps
    cifti_task_ID = ft_read_cifti_mod(task_UniqueIDs{z});

    cifti_rest_vars = ft_read_cifti_mod(rest_varmaps{z});       %% Load spatial correlation variant maps
    cifti_task_vars = ft_read_cifti_mod(task_varmaps{z});
    
  	if SizeExclude         %% Apply exclusion masks for unique ID maps
        
        cifti_rest_ID.data = variant_size_exclude(cifti_rest_ID,minsize);
        cifti_task_ID.data = variant_size_exclude(cifti_task_ID,minsize);
        
    end
    
    cifti_rest_rotations = ft_read_cifti_mod([ciftirotdir ciftirotoutfilerest]);   %% Load CIFTI rotation files for each subject
    cifti_task_rotations = ft_read_cifti_mod([ciftirotdir ciftirotoutfiletask]);
    
    % get mean absolute differences for vertices corresponding to actual variants
    
    TrueRestMags = [TrueRestMags mean(abs(cifti_rest_vars.data(cifti_rest_ID.data > 0) - cifti_task_vars.data(cifti_rest_ID.data > 0)))];  %% subtract the spatial correlations for the task variant map from the rest variant for the vertices corresponding to variants in rest
    TrueTaskMags = [TrueTaskMags mean(abs(cifti_rest_vars.data(cifti_task_ID.data > 0) - cifti_task_vars.data(cifti_task_ID.data > 0)))];  %% subtract the spatial correlations for the task variant map from the rest variant for the vertices corresponding to variants in task
    
    disp(['Calculating mean absolute magnitude difference for rotations for subject ' subject])
    
    for rots = 1:iterations     %% Get mean absolute difference for vertices in each rotation
        
        RestRotateMagsSub = [RestRotateMagsSub; mean(abs(cifti_rest_vars.data(cifti_rest_rotations.data(:,rots) > 0) - cifti_task_vars.data(cifti_rest_rotations.data(:,rots) > 0)))];  %% subtract the spatial correlations for the task variant map from the rest variant for the vertices corresponding to each variant rotation in rest
        TaskRotateMagsSub = [TaskRotateMagsSub; mean(abs(cifti_rest_vars.data(cifti_task_rotations.data(:,rots) > 0) - cifti_task_vars.data(cifti_task_rotations.data(:,rots) > 0)))];  %% subtract the spatial correlations for the task variant map from the rest variant for the vertices corresponding to each variant rotation in task
        
    end
    
    RestRotateMags = [RestRotateMags RestRotateMagsSub];        %% Add subject variable to final output matrix
    TaskRotateMags = [TaskRotateMags TaskRotateMagsSub];
    
end

% Do permutation tests for p-values

RestPerms = [];
TaskPerms = [];

for perm = 1:size(RestRotateMags,2)  %% get proportion of rotations for each subject that have a higher mean difference across states than the actual variants
    
    pvalrest = length(find(RestRotateMags(:,perm) >= TrueRestMags(perm)))/size(RestRotateMags,1);
    pvaltask = length(find(TaskRotateMags(:,perm) >= TrueTaskMags(perm)))/size(RestRotateMags,1);
    
    RestPerms = [RestPerms; pvalrest];
    TaskPerms = [TaskPerms; pvaltask];
    
end
    

%% Plot real and rotated variant magnitude differences for each state

positions = [1:10];      % Rest Data

scatter([sort(repmat([1:1:9],1,1000)) repmat(positions(10),1,size(RestRotateMags,1)*size(RestRotateMags,2))], [RestRotateMags(:,1)' RestRotateMags(:,2)' RestRotateMags(:,3)' RestRotateMags(:,4)' RestRotateMags(:,5)' RestRotateMags(:,6)' RestRotateMags(:,7)' RestRotateMags(:,8)' RestRotateMags(:,9)' reshape(RestRotateMags,1,size(RestRotateMags,1)*size(RestRotateMags,2))], 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
hold on
scatter([positions(1:9) repmat(positions(10),1,length(TrueRestMags))], [TrueRestMags TrueRestMags], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', 'r');
hold on
set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4) positions(5) positions(6) positions(7) positions(8) positions(9) positions(10)])
set(gca,'xticklabel',{'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'All Subjects'}, 'FontSize',18)
ylabel('Average Magnitude Difference', 'FontSize',18)
title('Rest Magnitude Differences')

c = get(gca, 'Children');

[l, hobj, hout, mout] = legend(c(1:2), 'Variant Magnitude Difference', 'Rotated Magnitude Differences', 'Location', 'NorthEast');
l.FontSize = 18;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',18);
ylim([0 .3]);
xlim([0 11]);

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.6, 0.6, 0.45]);


print(gcf,[outputdir 'Rest_Rotation_Maps_' num2str(threshold) '.jpg'],'-dpng','-r300');

close gcf

positions = [1:10];      % Task Data

scatter([sort(repmat([1:1:9],1,1000)) repmat(positions(10),1,size(TaskRotateMags,1)*size(TaskRotateMags,2))], [TaskRotateMags(:,1)' TaskRotateMags(:,2)' TaskRotateMags(:,3)' TaskRotateMags(:,4)' TaskRotateMags(:,5)' TaskRotateMags(:,6)' TaskRotateMags(:,7)' TaskRotateMags(:,8)' TaskRotateMags(:,9)' reshape(TaskRotateMags,1,size(TaskRotateMags,1)*size(TaskRotateMags,2))], 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
hold on
scatter([positions(1:9) repmat(positions(10),1,length(TrueRestMags))], [TrueTaskMags TrueTaskMags], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', 'r');
hold on
set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4) positions(5) positions(6) positions(7) positions(8) positions(9) positions(10)])
set(gca,'xticklabel',{'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'All Subjects'}, 'FontSize',18)
ylabel('Average Magnitude Difference', 'FontSize',16)
title('Task Magnitude Differences')

c = get(gca, 'Children');

% [l, hobj, hout, mout] = legend(c(1:2), 'Variant Magnitude Difference', 'Rotated Magnitude Differences', 'Location', 'NorthEast');
% l.FontSize = 18;
% s = findobj(hobj,'type','patch');
% set(s,'MarkerSize',10);
% x = findobj(hobj,'type','text');
% set(x,'FontSize',18);
ylim([0 .3]);
xlim([0 11]);

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.6, 0.6, 0.45]);

print(gcf,[outputdir 'Task_Rotation_Maps_' num2str(threshold) '.jpg'],'-dpng','-r300')

close gcf


    
function rotations = get_rotations(rotationfilepath,generaterotations,iterations)

%% load rotations and generate new ones if desired

rotations = load(rotationfilepath);  %% load file with 1000 pre-generated rotations

if generaterotations  %% if generating new random rotations
    
    rotationsxrot = [];  %% create vector for each plane (x, y, z) of rotation
    rotationsyrot = [];
    rotationszrot = [];
    
    for rots = 1:iterations  %% for each iteration
        
        rng('shuffle')  %% shuffle random number generator
        
        rotationsxrot = [rotationsxrot min(rotations.rotations.xrot)+(max(rotations.rotations.xrot)-min(rotations.rotations.xrot))*rand(1,1)];  %% create a new random value for each plane of the rotation between the minimum and maximum of all values in the original 1000 rotations for that plane
        rotationsyrot = [rotationsyrot min(rotations.rotations.yrot)+(max(rotations.rotations.yrot)-min(rotations.rotations.yrot))*rand(1,1)];
        rotationszrot = [rotationszrot min(rotations.rotations.zrot)+(max(rotations.rotations.zrot)-min(rotations.rotations.zrot))*rand(1,1)];
        
    end
    
    rotations.rotations.xrot = rotationsxrot;  %% change rotations in loaded file to new rotations
    rotations.rotations.yrot = rotationsyrot;
    rotations.rotations.zrot = rotationszrot;
    
    
end
end


function outfile = create_filename(sub,task,threshold,hemi,gifti,rotfile)

% set output filename for each map

ses = 'All';
grp = 'WashU120';
if isempty(hemi)
    desc = 'matched';
else
    desc = [hemi '-matched'];
end
thresh_str = num2str(threshold);
if gifti
    fileext = '.func.gii';
else
    fileext = '.dtseries.nii';
end
if rotfile
    filetype = 'RotationFileMap';
else
    filetype = 'UniqueIDMap';
end

outfile = sprintf(['sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_' filetype fileext],sub,task,ses,grp,desc,thresh_str);

end
    