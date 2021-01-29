
clear all

%% Compute and plot dice correlations of spatial location overlap for individual tasks, combined tasks, and rest
%
% This script loads uniqueID maps for split-halves of rest and task data
% and calculates the spatial location overlap for variants within states
% (rest vs. rest and each task vs. each task), between states (rest vs.
% each task), and across subjects (e.g., MSC01 each task vs. MSC02 rest).
% This is done separately for network variants for all 3 different tasks,
% rest, and all of the tasks combined. Spatial location overlap is
% calculated according to binarized maps (a 1 represents a vertex where a
% variant exists, and 0 a vertex where a variant does not exist). A
% Dice-Sorenson correlation is then performed on these maps to calculate
% spatial overlap. The values for the within state, between state, and
% across subject comparisons are then plotted by task.
%
% INPUTS:
% 
% -outputdir: the output directory for the plots of spatial location 
% overlap
% -minsize: the minimum size threshold (in vertices) that a variant has to
% meet in order to be included in this analysis
% -SizeExclude: toggles whether to exclude variants below a minimum size
% (in vertices) threshold (set to 1), or use all variant locations
% regardless of size (set to 0)
% -threshold: a threshold of uniqueID maps to load for measuring spatial 
% overlap
%
% (see below for more details on input text files)
%
% OUTPUTS:
%
% -plots: creates a plot for the spatial location overlap of variants
% within states, between states, and across subjects by each subject as
% well as combined for each task separately
%
% Written by BK (01-2021)
%

%% Initialize Variables

threshold = 5;  %% Toggles threshold to create variants
SizeExclude = 0;  %% Toggles whether to exclude variants smaller than criterion
minsize = 50;  %% Minimum variant size for size exclusion

outputdir = '/your/path/here/';  %% Output directory for plots


%% Load all data and assign it to the appropriate variables for comparing spatial overlap

memmemdicecoefs = [];  %% Set empty variables to hold correlation values
mixedmixeddicecoefs = [];
motormotordicecoefs = [];
memmemrestdicecoefs = [];
mixedmixedrestdicecoefs = [];
motormotorrestdicecoefs = [];
memrestmemrestdicecoefs = [];
mixedrestmixedrestdicecoefs = [];
motorrestmotorrestdicecoefs = [];
alltaskalltaskdicecoefs = [];
alltaskalltaskrestdicecoefs = [];
alltaskrestalltaskrestdicecoefs = [];

memfilestaskeven = [];  %% Set empty variables to hold vectors for dice correlation
memfilestaskodd = [];
mixedfilestaskeven = [];
mixedfilestaskodd = [];
motorfilestaskeven = [];
motorfilestaskodd = [];
memfilesresteven = [];
memfilesrestodd = [];
mixedfilesresteven = [];
mixedfilesrestodd = [];
motorfilesresteven = [];
motorfilesrestodd = [];
alltaskfilestaskeven = [];
alltaskfilestaskodd = [];
alltaskfilesresteven = [];
alltaskfilesrestodd = [];

comparisons = {'matched-samedata-task-task', 'matched-samedata-task-rest', 'matched-samedata-taskrest-taskrest', 'matched-samedata-alltask-alltask', 'matched-samedata-alltask-rest', 'matched-samedata-alltaskrest-alltaskrest'};  %% Toggles through comparisons
% 1. individual task to individual task
% 2. individual task to rest
% 3. rest to rest (number of samples per session matched to each individual task)
% 4. all tasks combined to all tasks combined
% 5. all tasks combined to rest
% 6. rest to rest (number of samples per session matched to all tasks combined)


for filetype = 1:4  %% Number of iterations for each file type (see labels below)
    
    % load uniqueID maps for analyses (see threshold_variant_maps.m for
    % additional documentation on uniqueID maps) using a space-delimited
    % text file in the format: pathtofile subID taskid
    % e.g. filepath/MSC01.dtseries.nii MSC01 memtask
    % the order of the data files for all of the subjects should be the
    % same in all text files

    if filetype == 1   %% Individual task files (split-halves of separated memory, mixed, and motor tasks)

        [even_files, subjects, tasks] = textread(['/your/path/here/MSC_task_even_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
        
        [odd_files, ~, ~] = textread(['/your/path/here/MSC_task_odd_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
        
    elseif filetype == 2   %% Rest files with sampling matched to the task data for individual tasks (split-halves of rest matched to separated memory, mixed, and motor tasks)

        [even_files, subjects, tasks] = textread(['/your/path/here/MSC_rest_task_even_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
        
        [odd_files, ~, ~] = textread(['/your/path/here/MSC_rest_task_odd_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
        
    elseif filetype == 3   %% All tasks combined files (split-halves of combined memory, mixed, and motor tasks)

        [even_files, subjects, tasks] = textread(['/your/path/here/MSC_alltask_even_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
        
        [odd_files, ~, ~] = textread(['/your/path/here/MSC_alltask_odd_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
        
    elseif filetype == 4 %% Rest files with sampling matched to the task data for the combined tasks (split-halves of rest matched to combined memory, mixed, and motor tasks)

        [even_files, subjects, tasks] = textread(['/your/path/here/MSC_rest_alltask_even_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
        
        [odd_files, ~, ~] = textread(['/your/path/here/MSC_rest_alltask_odd_' num2str(threshold) '_uniqueID.txt'],'%s%s%s');
        
    end

    for a = 1:length(even_files)  %% Loop through each file in text list

        even_file = even_files{a};
        odd_file = odd_files{a};
        task = tasks{a};

        cifti_even = ft_read_cifti_mod(even_file);      %% Load each file from textlist
        cifti_odd = ft_read_cifti_mod(odd_file);

        if SizeExclude == 1     %% Exclude variants below size threshold

            cifti_even.data = variant_size_exclude(cifti_even,minsize);
            cifti_odd.data = variant_size_exclude(cifti_odd,minsize);
            
        end
        
        cifti_even.data(cifti_even.data>0) = 1;  %% Binarize data
        cifti_odd.data(cifti_odd.data>0) = 1;
        
        % Store data for dice correlations after creating final variants
    
        if filetype == 1   %% Individual task files
            
            if strcmp(strtrim(task),'mem')
                
                memfilestaskeven = [memfilestaskeven; cifti_even.data'];
                memfilestaskodd = [memfilestaskodd; cifti_odd.data'];
                
            elseif strcmp(strtrim(task),'mixed')
                
                mixedfilestaskeven = [mixedfilestaskeven; cifti_even.data'];
                mixedfilestaskodd = [mixedfilestaskodd; cifti_odd.data'];
                
            elseif strcmp(strtrim(task),'motor')
                
                motorfilestaskeven = [motorfilestaskeven; cifti_even.data'];
                motorfilestaskodd = [motorfilestaskodd; cifti_odd.data'];
                
            end
            
        elseif filetype == 2   %% Individual task-matched rest files
                
            if strcmp(strtrim(task),'memrest')
                
                memfilesresteven = [memfilesresteven; cifti_even.data'];
                memfilesrestodd = [memfilesrestodd; cifti_odd.data'];
                
            elseif strcmp(strtrim(task),'mixedrest')
                
                mixedfilesresteven = [mixedfilesresteven; cifti_even.data'];
                mixedfilesrestodd = [mixedfilesrestodd; cifti_odd.data'];
                
            elseif strcmp(strtrim(task),'motorrest')
                
                motorfilesresteven = [motorfilesresteven; cifti_even.data'];
                motorfilesrestodd = [motorfilesrestodd; cifti_odd.data'];
                
            end
            
        elseif filetype == 3   %% All tasks combined files
            
            alltaskfilestaskeven = [alltaskfilestaskeven; cifti_even.data'];
            alltaskfilestaskodd = [alltaskfilestaskodd; cifti_odd.data'];
            
        elseif filetype == 4   %% All tasks combined-matched rest files
            
            alltaskfilesresteven = [alltaskfilesresteven; cifti_even.data'];
            alltaskfilesrestodd = [alltaskfilesrestodd; cifti_odd.data'];
            
        end
    end
end


%% Calculate spatial location overlap for appropriate comparisons

for g = 1:numel(comparisons)
    
    comp = comparisons(g);
    
    for dataset = 1:length(subjects)    %% Loop through each dataset

        if strcmp(comp, 'matched-samedata-task-task')   %% Spatial locations of each individual task versus matched individual task data

            dcorrdatamem = binarize_vectors(memfilestaskeven,memfilestaskodd,dataset);          %% Binarize vectors

            dcmem = dice_coefficient_mod(dcorrdatamem(:,1),dcorrdatamem(:,2));      %% Calculate dice correlation on binary data

            memmemdicecoefs = [memmemdicecoefs; dcmem];        %% Store dice values

            % Mixed Task

            dcorrdatamixed = binarize_vectors(mixedfilestaskeven,mixedfilestaskodd,dataset);

            dcmixed = dice_coefficient_mod(dcorrdatamixed(:,1),dcorrdatamixed(:,2));

            mixedmixeddicecoefs = [mixedmixeddicecoefs; dcmixed];

            % Motor Task

            dcorrdatamotor = binarize_vectors(motorfilestaskeven,motorfilestaskodd,dataset);

            dcmotor = dice_coefficient_mod(dcorrdatamotor(:,1),dcorrdatamotor(:,2));

            motormotordicecoefs = [motormotordicecoefs; dcmotor];

        elseif strcmp(comp, 'matched-samedata-task-rest')       %% Spatial locations of individual tasks versus rest (matched to sampling of individual tasks)

            % Mem Task

            dcorrdatamemresteven = binarize_vectors(memfilesresteven,memfilestaskodd,dataset);          %% Binarize vectors
            dcorrdatamemrestodd = binarize_vectors(memfilestaskeven,memfilesrestodd,dataset);

            dcmemresteven = dice_coefficient_mod(dcorrdatamemresteven(:,1),dcorrdatamemresteven(:,2));      %% Calculate dice correlation on binary data
            dcmemrestodd = dice_coefficient_mod(dcorrdatamemrestodd(:,1),dcorrdatamemrestodd(:,2));

            memmemrestdicecoefs = [memmemrestdicecoefs; dcmemresteven dcmemrestodd];        %% Store dice values

            % Mixed Task

            dcorrdatamixedresteven = binarize_vectors(mixedfilesresteven,mixedfilestaskodd,dataset);
            dcorrdatamixedrestodd = binarize_vectors(mixedfilestaskeven,mixedfilesrestodd,dataset);

            dcmixedresteven = dice_coefficient_mod(dcorrdatamixedresteven(:,1),dcorrdatamixedresteven(:,2));
            dcmixedrestodd = dice_coefficient_mod(dcorrdatamixedrestodd(:,1),dcorrdatamixedrestodd(:,2));

            mixedmixedrestdicecoefs = [mixedmixedrestdicecoefs; dcmixedresteven dcmixedrestodd];

            % Motor Task

            dcorrdatamotorresteven = binarize_vectors(motorfilesresteven,motorfilestaskodd,dataset);
            dcorrdatamotorrestodd = binarize_vectors(motorfilestaskeven,motorfilesrestodd,dataset);

            dcmotorresteven = dice_coefficient_mod(dcorrdatamotorresteven(:,1),dcorrdatamotorresteven(:,2));
            dcmotorrestodd = dice_coefficient_mod(dcorrdatamotorrestodd(:,1),dcorrdatamotorrestodd(:,2));

            motormotorrestdicecoefs = [motormotorrestdicecoefs; dcmotorresteven dcmotorrestodd];

        elseif strcmp(comp, 'matched-samedata-taskrest-taskrest')       %% Spatial locations of rest versus rest (matched to sampling of individual tasks)

            % Mem Task

            dcorrdatamemrest = binarize_vectors(memfilesresteven,memfilesrestodd,dataset);

            dcmemrest = dice_coefficient_mod(dcorrdatamemrest(:,1),dcorrdatamemrest(:,2));

            memrestmemrestdicecoefs = [memrestmemrestdicecoefs; dcmemrest];

            % Mixed Task

            dcorrdatamixedrest = binarize_vectors(mixedfilesresteven,mixedfilesrestodd,dataset);

            dcmixedrest = dice_coefficient_mod(dcorrdatamixedrest(:,1),dcorrdatamixedrest(:,2));

            mixedrestmixedrestdicecoefs = [mixedrestmixedrestdicecoefs; dcmixedrest];

            % Motor Task

            dcorrdatamotorrest = binarize_vectors(motorfilesresteven,motorfilesrestodd,dataset);

            dcmotorrest = dice_coefficient_mod(dcorrdatamotorrest(:,1),dcorrdatamotorrest(:,2));

            motorrestmotorrestdicecoefs = [motorrestmotorrestdicecoefs; dcmotorrest];

        elseif strcmp(comp, 'matched-samedata-alltask-alltask')         %% Spatial locations of all tasks combined versus all tasks combined

            dcorrdataalltasks = binarize_vectors(alltaskfilestaskeven,alltaskfilestaskodd,dataset);

            dcalltasks = dice_coefficient_mod(dcorrdataalltasks(:,1),dcorrdataalltasks(:,2));

            alltaskalltaskdicecoefs = [alltaskalltaskdicecoefs; dcalltasks];

        elseif strcmp(comp, 'matched-samedata-alltask-rest')         %% Spatial locations of all tasks combined versus rest (matched to sampling of all tasks combined)

            dcorrdataalltasksresteven = binarize_vectors(alltaskfilesresteven,alltaskfilestaskodd,dataset);
            dcorrdataalltasksrestodd = binarize_vectors(alltaskfilestaskeven,alltaskfilesrestodd,dataset);

            dcalltasksresteven = dice_coefficient_mod(dcorrdataalltasksresteven(:,1),dcorrdataalltasksresteven(:,2));
            dcalltasksrestodd = dice_coefficient_mod(dcorrdataalltasksrestodd(:,1),dcorrdataalltasksrestodd(:,2));

            alltaskalltaskrestdicecoefs = [alltaskalltaskrestdicecoefs; dcalltasksresteven dcalltasksrestodd];

        elseif strcmp(comp, 'matched-samedata-alltaskrest-alltaskrest')         %% Spatial locations of rest versus rest  (matched to sampling of all tasks combined)

            dcorrdataalltasksrest = binarize_vectors(alltaskfilesresteven,alltaskfilesrestodd,dataset);

            dcalltasksrest = dice_coefficient_mod(dcorrdataalltasksrest(:,1),dcorrdataalltasksrest(:,2));

            alltaskrestalltaskrestdicecoefs = [alltaskrestalltaskrestdicecoefs; dcalltasksrest];

        end
    end
end


%% Get all spatial permutations across subjects and states


nsubs = size(memfilestaskeven,1);

DiceCorrssimmem = [];       %% Create empty variables to store results
DiceCorrSubsmem = zeros(2*(nsubs-1),nsubs);
DiceCorrssimmixed = [];
DiceCorrSubsmixed = zeros(2*(nsubs-1),nsubs);
DiceCorrssimmotor = [];
DiceCorrSubsmotor = zeros(2*(nsubs-1),nsubs);
DiceCorrssimalltasks = [];
DiceCorrSubsalltasks = zeros(2*(nsubs-1),nsubs);

for splithalf = 1:2     %% Do both versions of split-half combinations
    
    for subiter = 1:nsubs   %% Iterate through subjects
        
        loopsubs = [1:nsubs];       %% Loop through each subject
        loopsubs(subiter) = [];     %% Remove each subject from across-subject permutation
        
        Count = 0;      %% Count number of loops to store data in correct place
        
        for subcompare = loopsubs    %% Get dice correlation for each participant to every other participant
            
            Count = Count+1;
            
            % Get binarized map for each participant and loop through
            % other participants to compare spatial location overlap across
            % subjects
            
            if splithalf == 1       %% First set of split-half comparisons
                
                memdatrest = memfilesrestodd(subiter,:);
                memdat = memfilestaskeven(subcompare,:);
                
                mixeddatrest = mixedfilesrestodd(subiter,:);
                mixeddat = mixedfilestaskeven(subcompare,:);
                
                motordatrest = motorfilesrestodd(subiter,:);
                motordat = motorfilestaskeven(subcompare,:);
                
                alltasksdatrest = alltaskfilesrestodd(subiter,:);
                alltasksdat = alltaskfilestaskeven(subcompare,:);
                
            else        %% Second set of split-half comparisons
                
                memdatrest = memfilesresteven(subiter,:);
                memdat = memfilestaskodd(subcompare,:);
                
                mixeddatrest = mixedfilesresteven(subiter,:);
                mixeddat = mixedfilestaskodd(subcompare,:);
                
                motordatrest = motorfilesresteven(subiter,:);
                motordat = motorfilestaskodd(subcompare,:);
                
                alltasksdatrest = alltaskfilesresteven(subiter,:);
                alltasksdat = alltaskfilestaskodd(subcompare,:);
                
            end
            
            % Across-subject permutations for mem task
            
            dcorrdatamem = binarize_vectors_perm(memdat,memdatrest);
            
            if isempty(dcorrdatamem)    %% Catch if no vertices overlap
                
                dcmem = 0;
                
            else
                
                dcmem = dice_coefficient_mod(dcorrdatamem(:,1),dcorrdatamem(:,2));
                
            end
            
            % Across-subject permutations for mixed task
            
            dcorrdatamixed = binarize_vectors_perm(mixeddat,mixeddatrest);
            
            if isempty(dcorrdatamixed)
                
                dcmixed = 0;
                
            else
                
                dcmixed = dice_coefficient_mod(dcorrdatamixed(:,1),dcorrdatamixed(:,2));
                
            end
            
            % Across-subject permutations for motor task
            
            dcorrdatamotor = binarize_vectors_perm(motordat,motordatrest);
            
            if isempty(dcorrdatamotor)
                
                dcmotor = 0;
                
            else
                
                dcmotor = dice_coefficient_mod(dcorrdatamotor(:,1),dcorrdatamotor(:,2));
                
            end
            
            % Across-subject permutations for combined tasks
            
            dcorrdataalltasks = binarize_vectors_perm(alltasksdat,alltasksdatrest);
            
            if isempty(dcorrdataalltasks)
                
                dcalltasks = 0;
                
            else
                
                dcalltasks = dice_coefficient_mod(dcorrdataalltasks(:,1),dcorrdataalltasks(:,2));
                
            end
            
            DiceCorrssimmem = [DiceCorrssimmem; dcmem];
            DiceCorrssimmixed = [DiceCorrssimmixed; dcmixed];
            DiceCorrssimmotor = [DiceCorrssimmotor; dcmotor];
            DiceCorrssimalltasks = [DiceCorrssimalltasks; dcalltasks];
            
            if splithalf == 1       %% Store values from dice correlations
                
                DiceCorrSubsmem(Count,subiter) = dcmem;
                DiceCorrSubsmixed(Count,subiter) = dcmixed;
                DiceCorrSubsmotor(Count,subiter) = dcmotor;
                DiceCorrSubsalltasks(Count,subiter) = dcalltasks;
                
            else
                
                DiceCorrSubsmem(Count+(nsubs-1),subiter) = dcmem;
                DiceCorrSubsmixed(Count+(nsubs-1),subiter) = dcmixed;
                DiceCorrSubsmotor(Count+(nsubs-1),subiter) = dcmotor;
                DiceCorrSubsalltasks(Count+(nsubs-1),subiter) = dcalltasks;
                
            end
            
            
        end
    end
end



%% Plot spatial overlap for each task

jitterAmount = 0.075;
jitterValuesX = 2*(rand(size(alltaskrestalltaskrestdicecoefs))-0.5)*jitterAmount;   % +/-jitterAmount max

positions = [1.5:.5:3];
scatterpos = [1.5:.5:2.5];

figure;
subplot(2,2,1)  %% All Tasks Plot
scatter(scatterpos+jitterValuesX(8), [alltaskrestalltaskrestdicecoefs(8,:) alltaskalltaskdicecoefs(8,:) mean(alltaskalltaskrestdicecoefs(8,:))], 150, 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
hold on
scatter(scatterpos+jitterValuesX(7), [alltaskrestalltaskrestdicecoefs(7,:) alltaskalltaskdicecoefs(7,:) mean(alltaskalltaskrestdicecoefs(7,:))], 150, 'filled', 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(6), [alltaskrestalltaskrestdicecoefs(6,:) alltaskalltaskdicecoefs(6,:) mean(alltaskalltaskrestdicecoefs(6,:))], 150, 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter(scatterpos+jitterValuesX(5), [alltaskrestalltaskrestdicecoefs(5,:) alltaskalltaskdicecoefs(5,:) mean(alltaskalltaskrestdicecoefs(5,:))], 150, 'filled', 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(4), [alltaskrestalltaskrestdicecoefs(4,:) alltaskalltaskdicecoefs(4,:) mean(alltaskalltaskrestdicecoefs(4,:))], 150, 'filled', 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter(scatterpos+jitterValuesX(3), [alltaskrestalltaskrestdicecoefs(3,:) alltaskalltaskdicecoefs(3,:) mean(alltaskalltaskrestdicecoefs(3,:))], 150, 'filled', 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter(scatterpos+jitterValuesX(2), [alltaskrestalltaskrestdicecoefs(2,:) alltaskalltaskdicecoefs(2,:) mean(alltaskalltaskrestdicecoefs(2,:))], 150, 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter(scatterpos+jitterValuesX(1), [alltaskrestalltaskrestdicecoefs(1,:) alltaskalltaskdicecoefs(1,:) mean(alltaskalltaskrestdicecoefs(1,:))], 150, 'filled', 'MarkerFaceColor', [0, 0, 0]);
hold on
scatter(repmat(positions(4),1,length(DiceCorrssimalltasks)), DiceCorrssimalltasks, 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4)])
set(gca,'xticklabel',{'Rest to Rest', 'Task to Task', 'Rest to Task', 'Across Subjects'}, 'FontSize',18)
ylabel('Dice Correlation', 'FontSize',18)
ax = gca;
ax.FontSize = 18;
ylim([0 1]);
xlim([1 3.5]);
title('All Tasks Combined')
subplot(2,2,2)  %% Memory Task Plot
scatter(scatterpos+jitterValuesX(8), [memrestmemrestdicecoefs(8,:) memmemdicecoefs(8,:) mean(memmemrestdicecoefs(8,:))], 150, 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
hold on
scatter(scatterpos+jitterValuesX(7), [memrestmemrestdicecoefs(7,:) memmemdicecoefs(7,:) mean(memmemrestdicecoefs(7,:))], 150, 'filled', 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(6), [memrestmemrestdicecoefs(6,:) memmemdicecoefs(6,:) mean(memmemrestdicecoefs(6,:))], 150, 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter(scatterpos+jitterValuesX(5), [memrestmemrestdicecoefs(5,:) memmemdicecoefs(5,:) mean(memmemrestdicecoefs(5,:))], 150, 'filled', 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(4), [memrestmemrestdicecoefs(4,:) memmemdicecoefs(4,:) mean(memmemrestdicecoefs(4,:))], 150, 'filled', 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter(scatterpos+jitterValuesX(3), [memrestmemrestdicecoefs(3,:) memmemdicecoefs(3,:) mean(memmemrestdicecoefs(3,:))], 150, 'filled', 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter(scatterpos+jitterValuesX(2), [memrestmemrestdicecoefs(2,:) memmemdicecoefs(2,:) mean(memmemrestdicecoefs(2,:))], 150, 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter(scatterpos+jitterValuesX(1), [memrestmemrestdicecoefs(1,:) memmemdicecoefs(1,:) mean(memmemrestdicecoefs(1,:))], 150, 'filled', 'MarkerFaceColor', [0, 0, 0]);
hold on
scatter(repmat(positions(4),1,length(DiceCorrssimmem)), DiceCorrssimmem, 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4)])
set(gca,'xticklabel',{'Rest to Rest', 'Task to Task', 'Rest to Task', 'Across Subjects'}, 'FontSize',18)
ylabel('Dice Correlation', 'FontSize',18)
ax = gca;
ax.FontSize = 18;
ylim([0 1]);
xlim([1 3.5]);
title('Memory Task')
subplot(2,2,3)  %% Mixed Task Plot
scatter(scatterpos+jitterValuesX(8), [mixedrestmixedrestdicecoefs(8,:) mixedmixeddicecoefs(8,:) mean(mixedmixedrestdicecoefs(8,:))], 150, 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
hold on
scatter(scatterpos+jitterValuesX(7), [mixedrestmixedrestdicecoefs(7,:) mixedmixeddicecoefs(7,:) mean(mixedmixedrestdicecoefs(7,:))], 150, 'filled', 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(6), [mixedrestmixedrestdicecoefs(6,:) mixedmixeddicecoefs(6,:) mean(mixedmixedrestdicecoefs(6,:))], 150, 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter(scatterpos+jitterValuesX(5), [mixedrestmixedrestdicecoefs(5,:) mixedmixeddicecoefs(5,:) mean(mixedmixedrestdicecoefs(5,:))], 150, 'filled', 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(4), [mixedrestmixedrestdicecoefs(4,:) mixedmixeddicecoefs(4,:) mean(mixedmixedrestdicecoefs(4,:))], 150, 'filled', 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter(scatterpos+jitterValuesX(3), [mixedrestmixedrestdicecoefs(3,:) mixedmixeddicecoefs(3,:) mean(mixedmixedrestdicecoefs(3,:))], 150, 'filled', 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter(scatterpos+jitterValuesX(2), [mixedrestmixedrestdicecoefs(2,:) mixedmixeddicecoefs(2,:) mean(mixedmixedrestdicecoefs(2,:))], 150, 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter(scatterpos+jitterValuesX(1), [mixedrestmixedrestdicecoefs(1,:) mixedmixeddicecoefs(1,:) mean(mixedmixedrestdicecoefs(1,:))], 150, 'filled', 'MarkerFaceColor', [0, 0, 0]);
hold on
scatter(repmat(positions(4),1,length(DiceCorrssimmixed)), DiceCorrssimmixed, 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4)])
set(gca,'xticklabel',{'Rest to Rest', 'Task to Task', 'Rest to Task', 'Across Subjects'}, 'FontSize',18)
ylabel('Dice Correlation', 'FontSize',18)
ax = gca;
ax.FontSize = 18;
ylim([0 1]);
xlim([1 3.5]);
title('Semantic/Coherence Task')
subplot(2,2,4)  %% Motor Task Plot
scatter(scatterpos+jitterValuesX(8), [motorrestmotorrestdicecoefs(8,:) motormotordicecoefs(8,:) mean(motormotorrestdicecoefs(8,:))], 150, 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
hold on
scatter(scatterpos+jitterValuesX(7), [motorrestmotorrestdicecoefs(7,:) motormotordicecoefs(7,:) mean(motormotorrestdicecoefs(7,:))], 150, 'filled', 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(6), [motorrestmotorrestdicecoefs(6,:) motormotordicecoefs(6,:) mean(motormotorrestdicecoefs(6,:))], 150, 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter(scatterpos+jitterValuesX(5), [motorrestmotorrestdicecoefs(5,:) motormotordicecoefs(5,:) mean(motormotorrestdicecoefs(5,:))], 150, 'filled', 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(4), [motorrestmotorrestdicecoefs(4,:) motormotordicecoefs(4,:) mean(motormotorrestdicecoefs(4,:))], 150, 'filled', 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter(scatterpos+jitterValuesX(3), [motorrestmotorrestdicecoefs(3,:) motormotordicecoefs(3,:) mean(motormotorrestdicecoefs(3,:))], 150, 'filled', 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter(scatterpos+jitterValuesX(2), [motorrestmotorrestdicecoefs(2,:) motormotordicecoefs(2,:) mean(motormotorrestdicecoefs(2,:))], 150, 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter(scatterpos+jitterValuesX(1), [motorrestmotorrestdicecoefs(1,:) motormotordicecoefs(1,:) mean(motormotorrestdicecoefs(1,:))], 150, 'filled', 'MarkerFaceColor', [0, 0, 0]);
hold on
scatter(repmat(positions(4),1,length(DiceCorrssimmotor)), DiceCorrssimmotor, 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4)])
set(gca,'xticklabel',{'Rest to Rest', 'Task to Task', 'Rest to Task', 'Across Subjects'}, 'FontSize',18)
ylabel('Dice Correlation', 'FontSize',18)
ax = gca;
ax.FontSize = 18;
ylim([0 1]);
xlim([1 3.5]);
title('Motor Task')

set(gcf, 'Units', 'Normalized', 'OuterPosition', [.1, .1, .9, .7]);

subiter = findobj(gca,'Type','scatter');

[l, hobj, hout, mout] = legend([subiter(2:9)], 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC10', 'Location', 'NorthEast');
l.FontSize = 18;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',18);

filename = ['/AllSubjects_Subplots_AllTasksComp_' num2str(threshold) '_Percent.jpg'];

print(gcf,[outputdir filename],'-dpng','-r300');

close gcf



%% Functions for binarizing variant maps

function outmat = binarize_vectors(evenvector,oddvector,dataset)

outmat = [];

for vertices = 1:size(oddvector,2)
    
    if oddvector(dataset,vertices) == 1 && evenvector(dataset,vertices) == 1
        
        outmat = [outmat;1 1];
        
    elseif oddvector(dataset,vertices) == 1
        
        outmat = [outmat;1 0];
        
    elseif evenvector(dataset,vertices) == 1
        
        outmat = [outmat;0 1];
        
    end
end


end


function outmatperm = binarize_vectors_perm(taskvector,restvector)

outmatperm = [];

for vertices = 1:size(taskvector,2)
    
    if restvector(vertices) == 1 && taskvector(vertices) == 1
        
        outmatperm = [outmatperm;1 1];
        
    elseif restvector(vertices) == 1
        
        outmatperm = [outmatperm;1 0];
        
    elseif taskvector(vertices) == 1 
        
        outmatperm = [outmatperm;0 1];
        
    end
end

end
