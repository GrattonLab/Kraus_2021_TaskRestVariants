
clear all

%% Plot reliability of spatial correlation maps and variant spatial locations for "true" vs. "test" split-halves
%
% This script loads variant maps from "true" and "test" split halves for 
% different times (in minutes) from the provided text files. One text file
% corresponds to the "true" halves and the other to the "test" halves. For 
% each "test" half, a spatial correlation is calculated between that map 
% and the "true" half. This is repeated for all "test" maps in order to 
% create a reliability curve for each subject. This same process is also 
% repeated for a binarized map after applying a threshold to the spatial 
% correlation maps in the "true" and "test" halves. The reliability curves 
% for the binarized variant locations and spatial map correlations are then
% plotted.
%
% INPUTS:
%
% -outputdir: the output directory for the plots of the reliability curves
% -SNRpath: the filepath to the directory that contains the SNR CIFTIs for
% each subject (see see threshold_variant_maps.m for additional
% documentation)
% -SNR_fstring: the filename and extension for the SNR map
%
% -times:﻿a vector of times (in minutes) to plot the reliability curves for
% -absthresh: toggles whether to use an absolute correlation (r)
% value to define variants (set to 1), otherwise thresholds are interpreted
% as the xth lowest percentage of correlation values (set to 0)
% -threshold: a value to use to define variants in each map (values lower
% than threshold)
% -TaskResids: toggles whether the underlying data was created using task
% data that are the residuals of a GLM (set to 1), or task data that still
% include the task activations (set to 0)
% -SNRExclude: toggles whether to exclude low signal regions from
% consideration as variants (set to 1), otherwise allow all vertices to be
% defined as variants (set to 0)
%
% -test_files: reads path to a text file containing the paths to each of
% the "test" halves from least amount of data to the most amount of data
% (ascending order) for each subject. The format is (pathtofile subjectID)
% and the file is space-delimited (see below for more formatting details)
% -true_files: reads path to a text file containing the paths to each of
% the "true" halves for each subject. The format is (pathtofile subjectID)
% and the file is space-delimited (see below for more formatting details)
%
% OUTPUTS:
%
% -plots: creates reliability plots for both the spatial correlation
% between each set of variant maps (split-halves) as well as the binary
% correlation for spatial location
%
% Written by BK (01-2021)
%

%% Initialize Variables

outputdir = '/your/path/here/';  %% output directory for reliability plots
SNRpath = '/your/path/here/';  %% filepath to SNR maps
SNR_fstring = 'SNRname.dtseries.nii';  %% filename and extension for SNR map

times = [5:5:100];  %% Vector of times to plot reliability for
absthresh = 1;  %% Toggles whether to use absolute (=1) or percent threshold for variants
threshold = 0.3;  %% Threshold to use for reliability
TaskResids = 1;  %% Toggles whether to use task residuals (=1) or pre-GLM task data
SNRExclude = 1;  %% Toggles whether to apply SNR Exclusion to maps

CiftiCorrs = {};  %% Store Pearson correlations for each subject
DiceCorrs = {};  %% Dice Dice correlations for each subject


%% Load Data
% load variant maps for analyses
% a space-delimited text file in the format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% times for each subject in the test_files text file should be listed in
% ascending order
% (e.g. [5;10;15])

if TaskResids

    [test_files, subjects1, ~] = textread('/your/path/here/MSC_task_odd_reliability_varmaps.txt','%s%s%s');

    [true_files, subjects2, ~] = textread('/your/path/here/MSC_task_even_full_varmaps.txt','%s%s%s');
    
else
    
    [test_files, subjects1, ~] = textread('/your/path/here/MSC_task_odd_preGLM_reliability_varmaps.txt','%s%s%s');

    [true_files, subjects2, ~] = textread('/your/path/here/MSC_task_even_preGLM_full_varmaps.txt','%s%s%s');
    
end


%% Loop through "true" half for each subject, compare it to each "test" half for that subject

for x = 1:length(subjects2)    %% Get unique subject IDs
    
    CiftiCorrSub = [];
    DiceCorrsSub = [];
    
    cifti_template = ft_read_cifti_mod(true_files{x});
    
    % exclude low SNR for "true" half
    
    if SNRExclude
        
        SNRmask = ft_read_cifti_mod([SNRpath subjects2{x} SNR_fstring]); 
            
     	SNRmask.data = SNRmask.data(1:59412,:);
                
     	LowSNR = find(SNRmask.data < 750);
        
        cifti_template.data(LowSNR) = NaN;
        
    end
    
    % threshold variant maps for "true" half
    
    if absthresh
        
        template_dice_inds = find(cifti_template.data < threshold);
        
    else
    
        template_thresh = prctile(cifti_template.data,threshold);

        template_dice_inds = find(cifti_template.data < template_thresh);
        
    end
    
    % binarize "true" half
    
    cifti_template_dice = zeros(size(cifti_template.data));
    
    cifti_template_dice(template_dice_inds) = 1;
    
    % match subject files for each "true" half to incremental files from
    % "test" half
    
    ciftiInds = [];
    
    for y = 1:length(subjects1)
    
        if strcmp(subjects2{x},subjects1{y})
        
            ciftiInds = [ciftiInds y];
        
        end
    end
    
    % read in "test" files one at a time and find reliability with "true"
    % half
    
    for z = 1:length(ciftiInds)
        
        cifti_file = ft_read_cifti_mod(test_files{ciftiInds(z)});
        
        % exclude low SNR for "test" half
        
        if SNRExclude
        
            cifti_file.data(LowSNR) = NaN;
        
        end
        
        % spatial correlation between "true" and "test" split-halves
        
        corr = corrcoef(cifti_template.data,cifti_file.data,'Rows','complete');
        
        CiftiCorrSub = [CiftiCorrSub corr(1,2)];
        
        % threshold variant maps for "test" half
        
        if absthresh
            
            cifti_dice_inds = find(cifti_file.data < threshold);
            
        else
        
            cifti_thresh = prctile(cifti_file.data,threshold);
    
            cifti_dice_inds = find(cifti_file.data < cifti_thresh);
            
        end
        
        % binarize "test" half
        
        cifti_dice = zeros(size(cifti_file.data));
        
        cifti_dice(cifti_dice_inds) = 1;
        
        % calculate spatial location overlap between "true" and "test"
        % split-halves
        
        DiceCorrsSub = [DiceCorrsSub dice_coefficient_mod(cifti_template_dice,cifti_dice)];
        
    end
    
    CiftiCorrs = [CiftiCorrs; CiftiCorrSub];
    DiceCorrs = [DiceCorrs; DiceCorrsSub];
    
end

%% Plot reliability data

% Plot Pearson correlations

figure;
plot(times,CiftiCorrs{9,:},'Color',[1, 0.5, 0],'LineWidth', 3)
ylim([0 1]);
hold on
plot(times(1:11),CiftiCorrs{8,:},'Color',[0, 0.6, 0.6],'LineWidth', 3)
hold on
plot(times,CiftiCorrs{7,:},'Color',[1, 0, 1],'LineWidth', 3)
hold on
plot(times,CiftiCorrs{6,:},'Color',[0.2, 1, 1],'LineWidth', 3)
hold on
plot(times,CiftiCorrs{5,:},'Color',[0, 0, 1],'LineWidth', 3)
hold on
plot(times,CiftiCorrs{4,:},'Color',[1, 0, 0],'LineWidth', 3)
hold on
plot(times(1:19),CiftiCorrs{3,:},'Color',[0, 1, 0],'LineWidth', 3)
hold on
plot(times,CiftiCorrs{2,:},'Color',[0.9, 0.9, 0],'LineWidth', 3)
hold on
plot(times,CiftiCorrs{1,:},'Color',[0, 0, 0],'LineWidth', 3)
hold on

ylabel('Pearson Correlation (r)');
xlabel('Time (Minutes)');
title('Spatial Correlation Map Reliability', 'FontSize',18)

m = findobj(gca,'Type','line');

hleg1 = legend(m(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'SouthEast');
hleg1.FontSize = 18;
ax = gca;
ax.FontSize = 18;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.5, 0.3, 0.45]);

print(gcf,[outputdir 'Spatial Correlation Reliability Task Data Pearson.jpg'],'-dpng','-r300');

close gcf

% Plot Dice Correlations

figure;
plot(times,DiceCorrs{9,:},'Color',[1, 0.5, 0],'LineWidth', 3)
ylim([0 1]);
hold on
plot(times(1:11),DiceCorrs{8,:},'Color',[0, 0.6, 0.6],'LineWidth', 3)
hold on
plot(times,DiceCorrs{7,:},'Color',[1, 0, 1],'LineWidth', 3)
hold on
plot(times,DiceCorrs{6,:},'Color',[0.2, 1, 1],'LineWidth', 3)
hold on
plot(times,DiceCorrs{5,:},'Color',[0, 0, 1],'LineWidth', 3)
hold on
plot(times,DiceCorrs{4,:},'Color',[1, 0, 0],'LineWidth', 3)
hold on
plot(times(1:19),DiceCorrs{3,:},'Color',[0, 1, 0],'LineWidth', 3)
hold on
plot(times,DiceCorrs{2,:},'Color',[0.9, 0.9, 0],'LineWidth', 3)
hold on
plot(times,DiceCorrs{1,:},'Color',[0, 0, 0],'LineWidth', 3)
hold on

ylabel('Dice Correlation');
xlabel('Time (Minutes)');
title('Variant Spatial Location Reliability', 'FontSize',18)

m = findobj(gca,'Type','line');

hleg1 = legend(m(1:9), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'SouthEast');
hleg1.FontSize = 18;
ax = gca;
ax.FontSize = 18;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.5, 0.3, 0.45]);

print(gcf,[outputdir 'Variant Location Reliability Task Data Dice ' num2str(threshold) '.jpg'],'-dpng','-r300');


close gcf

