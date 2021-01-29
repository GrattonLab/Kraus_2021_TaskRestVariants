
clear all

%% Load uniqueID overlap maps by state, then calculate the average activation for each task within each contiguous variant in each state
%
% This script loads uniqueID overlap maps (see BreakUpOverlapMaps.m for 
% more documentation on uniqueID overlap maps) and finds the mean task
% activation (vs. baseline) within each variant "state" (task variant 
% vertices only, rest variant vertices only, vertices identified in both 
% task and rest states). This also requires mean task activation maps for
% the MSC subjects. For each "state" uniqueID overlap map, the mean task 
% activation (z-score) for each task is calculated within each variant. 
% These activations are then plotted for each task by state, with the 
% values for each subject labeled by color.
%
% INPUTS:
%
% -outputdir: the output directory for the plots of mean task activation
% within each contiguous variant for each task
% -minsize: the minimum size threshold (in vertices) that a variant has to
% meet in order to be included in this analysis
% -SizeExclude: toggles whether to exclude variants below a minimum size
% (in vertices) threshold (set to 1), or use all variant locations
% regardless of size (set to 0)
% -threshold: a threshold of uniqueID maps to load and create spatial
% location maps for
%
% -overlap_only_files: reads path to a text file containing the paths to 
% each of the uniqueID overlap map files for variant vertices that were 
% found in both states. The format is (pathtofile subjectID) and the file 
% is space-delimited (see below for more formatting details)
% -rest_only_files: reads path to a text file containing the paths to 
% each of the uniqueID overlap map files for variant vertices that were 
% found only during rest. The format is (pathtofile subjectID) and the file 
% is space-delimited (see below for more formatting details)
% -task_only_files: reads path to a text file containing the paths to 
% each of the uniqueID overlap map files for variant vertices that were 
% found only during task. The format is (pathtofile subjectID) and the file 
% is space-delimited (see below for more formatting details)
% -mem_activation_files: reads path to a text file containing the paths to 
% each of the activation map files for the memory task. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
% -mixed_activation_files: reads path to a text file containing the paths to 
% each of the activation map files for the mixed task. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
% -motor_activation_files: reads path to a text file containing the paths to 
% each of the activation map files for the motor task. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
%
% OUTPUTS:
%
% -plots: creates a plot for the average task activation for the vertices 
% of each variant ineah "state" by task. The variants in the plot are color
% coded by subject IDs
%
% Written by BK (01-2021)
%

%% Initialize variables

outdir = '/your/path/here/';  %% output directory for files
threshold = 5;  %% the threshold of the maps to use for variant subunits
SizeExclusion = 1;  %% toggles whether to exclude variants based on a minimum size thereshold
minsize = 30;  %% minimum size threshold (in vertices) for each variant subunit

% UniqueID overlap maps
% load uniqueID overlap maps for analyses (see BreakUpOverlapMaps.m for
% additional documentation on uniqueID maps) using a space-delimited text 
% file in the format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% the order of the data files for all of the subjects should be the same in
% all text files

[overlap_only_files, subjects, ~] = textread(['/your/path/here/MSC_overlap_only_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');

[rest_only_files, ~, ~] = textread(['/your/path/here/MSC_overlap_rest_only_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');

[task_only_files, ~, ~] = textread(['/your/path/here/MSC_overlap_task_only_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');

% Task activation maps
% load task activation maps for analyses (preprocessed data and events 
% available on OpenNeuro) using a space-delimited text file in the 
% format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% the order of the data files for all of the subjects should be the same in
% all text files

[mem_activation_files, ~, ~] = textread('/your/path/here/MSC_memtask_activations.txt','%s%s%s');

[mixed_activation_files, ~, ~] = textread('/your/path/here/MSC_mixedtask_activations.txt','%s%s%s');

[motor_activation_files, ~, ~] = textread('/your/path/here/MSC_motortask_activations.txt','%s%s%s');

nfiles = length(overlap_only_files);  %% get number of subjects from text file

mean_mem_activations_overlap = {};  %% store mean task activations for each variant for all subjects
mean_mem_activations_rest = {};
mean_mem_activations_task = {};
mean_mixed_activations_overlap = {};
mean_mixed_activations_rest = {};
mean_mixed_activations_task = {};
mean_motor_activations_overlap = {};
mean_motor_activations_rest = {};
mean_motor_activations_task = {};

for x = 1:nfiles
    
    %% Load Files
    
	subject = subjects{x};
    
    mean_mem_activations_overlap_sub = [];  %% store mean task activations for each variant for each subject
    mean_mem_activations_rest_sub = [];
    mean_mem_activations_task_sub = [];
    mean_mixed_activations_overlap_sub = [];
    mean_mixed_activations_rest_sub = [];
    mean_mixed_activations_task_sub = [];
    mean_motor_activations_overlap_sub = [];
    mean_motor_activations_rest_sub = [];
    mean_motor_activations_task_sub = [];
    
    disp(['Processing data for subject ' subject])
        
	cifti_overlap_only = ft_read_cifti_mod(overlap_only_files{x});  %% load each type of vertex subunit map and task activation map
	cifti_rest_only = ft_read_cifti_mod(rest_only_files{x});
    cifti_task_only = ft_read_cifti_mod(task_only_files{x});
	cifti_mem_activations = ft_read_cifti_mod(mem_activation_files{x});
	cifti_mixed_activations = ft_read_cifti_mod(mixed_activation_files{x});
    cifti_motor_activations = ft_read_cifti_mod(motor_activation_files{x});
    
    %% Size Exclusion for variant subunits smaller than a given threshold
    
    if SizeExclusion
    
        cifti_overlap_only.data = variant_size_exclude(cifti_overlap_only,minsize);  %% remove variants with fewer vertices than threshold
        cifti_rest_only.data = variant_size_exclude(cifti_rest_only,minsize);
        cifti_task_only.data = variant_size_exclude(cifti_task_only,minsize);
    
    end
    
    %% Match Variants to network
    
    vars_overlap_only = unique(cifti_overlap_only.data(cifti_overlap_only.data>0));  %% get unique parcel indices for each vertex type
    vars_rest_only = unique(cifti_rest_only.data(cifti_rest_only.data>0));
    vars_task_only = unique(cifti_task_only.data(cifti_task_only.data>0));

    %%% Overlap Variants (occur in both task and rest states)
    
    %%% loop through variants, get avg seed map, match seed map indices to task activation maps, get mean task activation for each set of vertices %%%
    for var = 1:length(vars_overlap_only)

        inds = find(cifti_overlap_only.data == vars_overlap_only(var));
        
        mean_mem_activations_overlap_sub = [mean_mem_activations_overlap_sub mean(cifti_mem_activations.data(inds))];
        mean_mixed_activations_overlap_sub = [mean_mixed_activations_overlap_sub mean(cifti_mixed_activations.data(inds))];
        mean_motor_activations_overlap_sub = [mean_motor_activations_overlap_sub mean(cifti_motor_activations.data(inds))];
        
    end
    
    
    %%% Rest-Only Variants (only occur during rest)
    
    %%% loop through variants, get avg seed map, match seed map indices to task activation maps, get mean task activation for each set of vertices %%%
    for var = 1:length(vars_rest_only)

        inds = find(cifti_rest_only.data == vars_rest_only(var));
        
        mean_mem_activations_rest_sub = [mean_mem_activations_rest_sub mean(cifti_mem_activations.data(inds))];
        mean_mixed_activations_rest_sub = [mean_mixed_activations_rest_sub mean(cifti_mixed_activations.data(inds))];
        mean_motor_activations_rest_sub = [mean_motor_activations_rest_sub mean(cifti_motor_activations.data(inds))];
        
    end
    
    
    %%% Task-Only Variants (only occur during task)
    
    %%% loop through variants, get avg seed map, match seed map indices to task activation maps, get mean task activation for each set of vertices %%%
    for var = 1:length(vars_task_only)

        inds = find(cifti_task_only.data == vars_task_only(var));
        
        mean_mem_activations_task_sub = [mean_mem_activations_task_sub mean(cifti_mem_activations.data(inds))];
        mean_mixed_activations_task_sub = [mean_mixed_activations_task_sub mean(cifti_mixed_activations.data(inds))];
        mean_motor_activations_task_sub = [mean_motor_activations_task_sub mean(cifti_motor_activations.data(inds))];
        
    end
    
    mean_mem_activations_overlap = [mean_mem_activations_overlap; mean_mem_activations_overlap_sub];  %% add variant task activations for each subject to rest of subjects' data
    mean_mem_activations_rest = [mean_mem_activations_rest; mean_mem_activations_rest_sub];
    mean_mem_activations_task = [mean_mem_activations_task; mean_mem_activations_task_sub];
    mean_mixed_activations_overlap = [mean_mixed_activations_overlap; mean_mixed_activations_overlap_sub];
    mean_mixed_activations_rest = [mean_mixed_activations_rest; mean_mixed_activations_rest_sub];
    mean_mixed_activations_task = [mean_mixed_activations_task; mean_mixed_activations_task_sub];
    mean_motor_activations_overlap = [mean_motor_activations_overlap; mean_motor_activations_overlap_sub];
    mean_motor_activations_rest = [mean_motor_activations_rest; mean_motor_activations_rest_sub];
    mean_motor_activations_task = [mean_motor_activations_task; mean_motor_activations_task_sub];


end

%% Plot Data


% Memory Data


figure;
scatterpos = [.6:1:4];
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{9})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{9})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{9}))], [mean_mem_activations_overlap{9} mean_mem_activations_rest{9} mean_mem_activations_task{9}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0.5, 0]);
xlim([0 4]);
ylim([-4 5]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{8})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{8})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{8}))], [mean_mem_activations_overlap{8} mean_mem_activations_rest{8} mean_mem_activations_task{8}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0.6, 0.6]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{7})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{7})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{7}))], [mean_mem_activations_overlap{7} mean_mem_activations_rest{7} mean_mem_activations_task{7}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{6})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{6})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{6}))], [mean_mem_activations_overlap{6} mean_mem_activations_rest{6} mean_mem_activations_task{6}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{5})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{5})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{5}))], [mean_mem_activations_overlap{5} mean_mem_activations_rest{5} mean_mem_activations_task{5}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{4})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{4})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{4}))], [mean_mem_activations_overlap{4} mean_mem_activations_rest{4} mean_mem_activations_task{4}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{3})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{3})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{3}))], [mean_mem_activations_overlap{3} mean_mem_activations_rest{3} mean_mem_activations_task{3}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{2})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{2})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{2}))], [mean_mem_activations_overlap{2} mean_mem_activations_rest{2} mean_mem_activations_task{2}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mem_activations_overlap{1})) repmat(scatterpos(2),1,numel(mean_mem_activations_rest{1})) repmat(scatterpos(3),1,numel(mean_mem_activations_task{1}))], [mean_mem_activations_overlap{1} mean_mem_activations_rest{1} mean_mem_activations_task{1}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0, 0]);
hold on
set(gca,'xtick',[scatterpos(1) scatterpos(2) scatterpos(3)])
set(gca,'xticklabel',{'Overlapping', 'Rest Only', 'Task Only'}, 'FontSize',24)
ylabel('Mean Task Activation (z)');
m = findobj(gca,'Type','scatter');
[l, hobj, ~, ~] = legend(m(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
l.FontSize = 18;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',18);
ax = gca;
ax.FontSize = 24;
title('Memory Task', 'FontSize',24)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.35, 0.35, 0.5]);


saveas(gcf,[outdir 'Variant_Task_Activation_Memory_ByVariant_' num2str(threshold) '_Percent.tiff'],'tiff');

close gcf


% Mixed Data


figure;
scatterpos = [.6:1:4];
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{9})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{9})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{9}))], [mean_mixed_activations_overlap{9} mean_mixed_activations_rest{9} mean_mixed_activations_task{9}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0.5, 0]);
xlim([0 4]);
ylim([-4 5]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{8})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{8})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{8}))], [mean_mixed_activations_overlap{8} mean_mixed_activations_rest{8} mean_mixed_activations_task{8}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0.6, 0.6]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{7})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{7})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{7}))], [mean_mixed_activations_overlap{7} mean_mixed_activations_rest{7} mean_mixed_activations_task{7}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{6})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{6})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{6}))], [mean_mixed_activations_overlap{6} mean_mixed_activations_rest{6} mean_mixed_activations_task{6}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{5})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{5})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{5}))], [mean_mixed_activations_overlap{5} mean_mixed_activations_rest{5} mean_mixed_activations_task{5}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{4})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{4})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{4}))], [mean_mixed_activations_overlap{4} mean_mixed_activations_rest{4} mean_mixed_activations_task{4}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{3})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{3})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{3}))], [mean_mixed_activations_overlap{3} mean_mixed_activations_rest{3} mean_mixed_activations_task{3}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{2})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{2})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{2}))], [mean_mixed_activations_overlap{2} mean_mixed_activations_rest{2} mean_mixed_activations_task{2}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_mixed_activations_overlap{1})) repmat(scatterpos(2),1,numel(mean_mixed_activations_rest{1})) repmat(scatterpos(3),1,numel(mean_mixed_activations_task{1}))], [mean_mixed_activations_overlap{1} mean_mixed_activations_rest{1} mean_mixed_activations_task{1}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0, 0]);
hold on
set(gca,'xtick',[scatterpos(1) scatterpos(2) scatterpos(3)])
set(gca,'xticklabel',{'Overlapping', 'Rest Only', 'Task Only'}, 'FontSize',24)
ylabel('Mean Task Activation (z)');
m = findobj(gca,'Type','scatter');
[l, hobj, ~, ~] = legend(m(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
l.FontSize = 18;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',18);
ax = gca;
ax.FontSize = 24;
title('Semantic/Coherence Task', 'FontSize',24)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.35, 0.35, 0.5]);


saveas(gcf,[outdir 'Variant_Task_Activation_Mixed_ByVariant_' num2str(threshold) '_Percent.tiff'],'tiff');

close gcf


% Motor Data


figure;
scatterpos = [.6:1:4];
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{9})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{9})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{9}))], [mean_motor_activations_overlap{9} mean_motor_activations_rest{9} mean_motor_activations_task{9}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0.5, 0]);
xlim([0 4]);
ylim([-4 5]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{8})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{8})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{8}))], [mean_motor_activations_overlap{8} mean_motor_activations_rest{8} mean_motor_activations_task{8}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0.6, 0.6]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{7})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{7})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{7}))], [mean_motor_activations_overlap{7} mean_motor_activations_rest{7} mean_motor_activations_task{7}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{6})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{6})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{6}))], [mean_motor_activations_overlap{6} mean_motor_activations_rest{6} mean_motor_activations_task{6}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{5})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{5})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{5}))], [mean_motor_activations_overlap{5} mean_motor_activations_rest{5} mean_motor_activations_task{5}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{4})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{4})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{4}))], [mean_motor_activations_overlap{4} mean_motor_activations_rest{4} mean_motor_activations_task{4}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{3})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{3})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{3}))], [mean_motor_activations_overlap{3} mean_motor_activations_rest{3} mean_motor_activations_task{3}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{2})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{2})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{2}))], [mean_motor_activations_overlap{2} mean_motor_activations_rest{2} mean_motor_activations_task{2}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter([repmat(scatterpos(1),1,numel(mean_motor_activations_overlap{1})) repmat(scatterpos(2),1,numel(mean_motor_activations_rest{1})) repmat(scatterpos(3),1,numel(mean_motor_activations_task{1}))], [mean_motor_activations_overlap{1} mean_motor_activations_rest{1} mean_motor_activations_task{1}], 150, 'filled', 'jitter','on', 'jitterAmount', 0.1, 'MarkerFaceColor', [0, 0, 0]);
hold on
set(gca,'xtick',[scatterpos(1) scatterpos(2) scatterpos(3)])
set(gca,'xticklabel',{'Overlapping', 'Rest Only', 'Task Only'}, 'FontSize',24)
ylabel('Mean Task Activation (z)');
m = findobj(gca,'Type','scatter');
[l, hobj, ~, ~] = legend(m(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
l.FontSize = 24;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',24);
ax = gca;
ax.FontSize = 24;
title('Motor Task', 'FontSize',24)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.35, 0.35, 0.5]);


saveas(gcf,[outdir 'Variant_Task_Activation_Motor_ByVariant_' num2str(threshold) '_Percent.tiff'],'tiff');

close gcf
