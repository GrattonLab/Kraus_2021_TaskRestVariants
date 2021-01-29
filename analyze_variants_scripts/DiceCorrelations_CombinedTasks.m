
clear all

%% Calculates the spatial location overlap of variants within states, between states, and across subjects
%
% This script loads uniqueID maps for split-halves of rest and task data
% and calculates the spatial location overlap for variants within states
% (rest vs. rest and task vs. task), between states (rest vs. task), and 
% across subjects (e.g., MSC01 task vs. MSC02 rest). Spatial location
% overlap is calculated according to binarized maps (a 1 represents a
% vertex where a variant exists, and 0 a vertex where a variant does not
% exist). A Dice-Sorenson correlation is then performed on these maps to
% calculate spatial overlap. The values for the within state, between
% state, and across subject comparisons are then plotted by subject.
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
% -TaskResids: toggles whether the underlying data was created using task
% data that are the residuals of a GLM (set to 1), or task data that still
% include the task activations (set to 0)
% -thresholds: a vector (or single value) of threshold(s) to load for
% measuring spatial overlap
%
% -task_files_even: reads path to a text file containing the paths to each
% of the uniqueID split-half files for task data in even-numbered sessions.
% The format is (pathtofile subjectID) and the file is space-delimited (see
% below for more formatting details)
% -task_files_odd: reads path to a text file containing the paths to each
% of the uniqueID split-half files for task data in odd-numbered sessions.
% The format is (pathtofile subjectID) and the file is space-delimited (see
% below for more formatting details)
% -rest_files_even: reads path to a text file containing the paths to each
% of the uniqueID split-half files for rest data in even-numbered sessions.
% The format is (pathtofile subjectID) and the file is space-delimited (see
% below for more formatting details)
% -rest_files_odd: reads path to a text file containing the paths to each
% of the uniqueID split-half files for rest data in odd-numbered sessions.
% The format is (pathtofile subjectID) and the file is space-delimited (see
% below for more formatting details)
%
% OUTPUTS:
%
% -plots: creates a plot for the spatial location overlap of variants
% within states, between states, and across subjects by each subject as
% well as combined
%
% Written by BK (01-2021)
%


%% Initialize Variables

outputdir = '/your/path/here/';  %% output directory for summary plot of spatial overlap

SizeExclude = 1;  %% Toggles whether to exclude small variants
minsize = 50;     %% Minimum size variants for size exclusion
TaskResids = 1;   %% Toggles whether to use residualized task data (=1) or pre-GLM data
thresholds = 5;  %% Sets the thresholds of files to load

%% Loop through all input thresholds 

for v = 1:numel(thresholds)
    
    %% Load in Unique ID files for each split-half in both states
    % load uniqueID maps for analyses (see threshold_variant_maps.m for
    % additional documentation on uniqueID maps) using a space-delimited
    % text file in the format: pathtofile subID
    % e.g. filepath/MSC01.dtseries.nii MSC01
    % the order of the data files for all of the subjects should be the
    % same in all text files

    if TaskResids
        
        [task_files_even, ~, ~] = textread(['/your/path/here/MSC_alltask_varmaps_' num2str(thresholds(v)) '_even_uniqueID.txt'],'%s%s%s');
        
        [rest_files_even, ~, ~] = textread(['/your/path/here/MSC_rest_alltask_varmaps_' num2str(thresholds(v)) '_even_uniqueID.txt'],'%s%s%s');
        
        [task_files_odd, ~, ~] = textread(['/your/path/here/MSC_alltask_varmaps_' num2str(thresholds(v)) '_odd_uniqueID.txt'],'%s%s%s');
        
        [rest_files_odd, subjects, ~] = textread(['/your/path/here/MSC_rest_alltask_varmaps_' num2str(thresholds(v)) '_odd_uniqueID.txt'],'%s%s%s');
        
    else
        
        [task_files_even, ~, ~] = textread(['/your/path/here/MSC_alltask_varmaps_' num2str(thresholds(v)) '_even_preGLM_uniqueID.txt'],'%s%s%s');
        
        [rest_files_even, ~, ~] = textread(['/your/path/here/MSC_rest_varmaps_' num2str(thresholds(v)) '_even_preGLM_uniqueID.txt'],'%s%s%s');
        
        [task_files_odd, ~, ~] = textread(['/your/path/here/MSC_alltask_varmaps_' num2str(thresholds(v)) '_odd_preGLM_uniqueID.txt'],'%s%s%s');
        
        [rest_files_odd, subjects, ~] = textread(['/your/path/here/MSC_rest_varmaps_' num2str(thresholds(v)) '_odd_preGLM_uniqueID.txt'],'%s%s%s');
        
    end
    
    %% Loop through all subjects and find spatial location overlap within and between states
    
    alltaskfilestaskeven = [];
    alltaskfilestaskodd = [];
    alltaskfilesresteven = [];
    alltaskfilesrestodd = [];
    
    DiceCorrsTaskRest = [];
    DiceCorrsTaskTask = [];
    DiceCorrsRestRest = [];
    
    nfiles = length(rest_files_even);
    
    for x = 1:nfiles
        
        subject = subjects{x};
        
        cifti_rest_even = ft_read_cifti_mod(rest_files_even{x});
        cifti_task_even = ft_read_cifti_mod(task_files_even{x});
        cifti_rest_odd = ft_read_cifti_mod(rest_files_odd{x});
        cifti_task_odd = ft_read_cifti_mod(task_files_odd{x});
        
        if SizeExclude == 1         %% apply size exclusion for variants smaller than threshold
            
            cifti_rest_even.data = variant_size_exclude(cifti_rest_even,minsize);
            cifti_rest_odd.data = variant_size_exclude(cifti_rest_odd,minsize);
            cifti_task_odd.data = variant_size_exclude(cifti_task_odd,minsize);
            cifti_task_even.data = variant_size_exclude(cifti_task_even,minsize);
            
        end
        
        % binarize vectors for spatial correlation
        
        dcorrdatataskrest1 = binarize_vectors(cifti_rest_even.data,cifti_task_odd.data);  %% Task-Rest Comparison 1
        dcorrdatataskrest2 = binarize_vectors(cifti_rest_odd.data,cifti_task_even.data);  %% Task-Rest Comparison 2
        dcorrdatarestrest = binarize_vectors(cifti_rest_odd.data,cifti_rest_even.data);  %% Rest-Rest Comparison
        dcorrdatatasktask = binarize_vectors(cifti_task_odd.data,cifti_task_even.data);  %% Task-Task Comparison
        
        % check if any binarized vectors are empty, else dice correlation
        % throws an error
        
        if isempty(dcorrdatataskrest1)
            dctaskrest1 = 0;
        else
            dctaskrest1 = dice_coefficient_mod(dcorrdatataskrest1(:,1),dcorrdatataskrest1(:,2));
        end
        
        if isempty(dcorrdatataskrest2)
            dctaskrest2 = 0;
        else
            dctaskrest2 = dice_coefficient_mod(dcorrdatataskrest2(:,1),dcorrdatataskrest2(:,2));
        end
        
        if isempty(dcorrdatarestrest)
            dcrestrest = 0;
        else
            dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
        end
        
        if isempty(dcorrdatatasktask)
            dctasktask = 0;
        else
            dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
        end
        
        DiceCorrsTaskRest = [DiceCorrsTaskRest; [dctaskrest1 dctaskrest2]];  % add data from each subject to vector
        DiceCorrsTaskTask = [DiceCorrsTaskTask; dctasktask];
        DiceCorrsRestRest = [DiceCorrsRestRest; dcrestrest];
        
        alltaskfilestaskeven = [alltaskfilestaskeven; cifti_task_even.data'];
        alltaskfilestaskodd = [alltaskfilestaskodd; cifti_task_odd.data'];
        alltaskfilesresteven = [alltaskfilesresteven; cifti_rest_even.data'];
        alltaskfilesrestodd = [alltaskfilesrestodd; cifti_rest_odd.data'];
        
        
    end
    
    %% Run permutations for spatial location overlap across subjects
    
    nsubs = nfiles;
    
    DiceCorrssim = [];  % store dice coefficients across all subjects
    DiceCorrSubs = zeros(2*(nsubs-1),nsubs);  % store comparison for each subject with every other subject
    
    for l = 1:2     %% loop through both versions of split-half combinations
        
        for m = 1:nsubs  %% for each subject
            
            loopvals = [1:nsubs];  % get indices of all subjects
            loopvals(m) = [];  % leave out current subject
            
            for n = 1:length(loopvals)  %% for the indices of all other subjects

                if l == 1  % if first set of split-halves
                    
                    restdat = alltaskfilesrestodd(m,:);  % select rest split-half for current subject
                    taskdat = alltaskfilestaskeven(loopvals(n),:);  % select a different subject's task split-half
                    
                else  % else if second set of split-halves
                    
                    restdat = alltaskfilesresteven(m,:);  % select rest split-half for current subject
                    taskdat = alltaskfilestaskodd(loopvals(n),:);  % select a different subject's task split-half
                    
                end

                dcorrdata = binarize_vectors(restdat,taskdat);  % binarize task and rest data
                
                % check if any binarized vectors are empty, else dice correlation
                % throws an error
                
                if isempty(dcorrdata)
                    dc = 0;
                else
                    dc = dice_coefficient_mod(dcorrdata(:,1),dcorrdata(:,2));
                end
                
                DiceCorrssim = [DiceCorrssim; dc];
                
                if l == 1   %% store across subject correlation in the appropriate column for plotting
                    
                    DiceCorrSubs(n,m) = dc;
                    
                else
                    
                    DiceCorrSubs(n+(nsubs-1),m) = dc;
                    
                end

            end
        end
    end
    
    DiceCorrspval = prctile(DiceCorrssim, 95);  %% 95th percentile of all permutations across subjects
    
    %% Plot spatial location overlap within each subject and across subjects
    
    x = [];
    group = [];
    positions = [1 2 3 4 5 6 7 8 9 10 10.25 10.5];
    positionssubs = [1 2 3 4 5 6 7 8 9 10];
    scatterposbetween = [sort(repmat([1:1:9],1,2))];
    scatterposwithin = [sort(repmat([1:1:9],1,2))];
    scatterposacross = [sort(repmat([1:1:9],1,16))];
    
    for g = 1:nsubs+1
        
        if g < nsubs+1

            x = [x 0];
            
        else
            
            x = [x [reshape(DiceCorrsTaskRest,1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2)) DiceCorrsTaskTask' DiceCorrsRestRest' DiceCorrssim']];
            
        end
    end
    
    for i = 1:length(positionssubs)
        
        if i < length(positionssubs)
            
            group = [group positionssubs(i)];
            
        else
            
            group = [group [repmat(positionssubs(i),1,size(DiceCorrsTaskRest,1)*size(DiceCorrsTaskRest,2)) repmat(positionssubs(i)+.25,1,size(DiceCorrsTaskTask,1)+size(DiceCorrsRestRest,1)) repmat(positionssubs(i)+.5,1,size(DiceCorrssim,1),1)]];
            
        end
    end
    
    boxplot(x,group, 'BoxStyle', 'outline', 'positions', positions, 'Symbol', '');
    hold on
    scatter(scatterposbetween, [DiceCorrsTaskRest(1,:) DiceCorrsTaskRest(2,:) DiceCorrsTaskRest(3,:) DiceCorrsTaskRest(4,:) DiceCorrsTaskRest(5,:) DiceCorrsTaskRest(6,:) DiceCorrsTaskRest(7,:) DiceCorrsTaskRest(8,:) DiceCorrsTaskRest(9,:)], 150, 'filled', 'MarkerFaceColor', 'b');
    hold on
    scatter(scatterposwithin, [DiceCorrsRestRest(1,:) DiceCorrsTaskTask(1,:) DiceCorrsRestRest(2,:) DiceCorrsTaskTask(2,:) DiceCorrsRestRest(3,:) DiceCorrsTaskTask(3,:) DiceCorrsRestRest(4,:) DiceCorrsTaskTask(4,:) DiceCorrsRestRest(5,:) DiceCorrsTaskTask(5,:) DiceCorrsRestRest(6,:) DiceCorrsTaskTask(6,:) DiceCorrsRestRest(7,:) DiceCorrsTaskTask(7,:) DiceCorrsRestRest(8,:) DiceCorrsTaskTask(8,:) DiceCorrsRestRest(9,:) DiceCorrsTaskTask(9,:)], 150, 'filled', 'MarkerFaceColor', 'g');
    hold on
    scatter(scatterposacross, [DiceCorrSubs(:,1)' DiceCorrSubs(:,2)' DiceCorrSubs(:,3)' DiceCorrSubs(:,4)' DiceCorrSubs(:,5)' DiceCorrSubs(:,6)' DiceCorrSubs(:,7)' DiceCorrSubs(:,8)' DiceCorrSubs(:,9)'], 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);
    hold on
    set(gca,'xtick',[positions(1) positions(2) positions(3) positions(4) positions(5) positions(6) positions(7) positions(8) positions(9) mean(positions(10:12))])
    set(gca,'xticklabel',{'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Group Average'}, 'FontSize',18)
    ylabel('Dice Correlation', 'FontSize',18)
    title(['Task Residuals ' num2str(thresholds(v)) ' Threshold'], 'FontSize',18)
    
    color = ['k', 'g', 'b'];
    
    h = findobj(gca,'Tag','Box');
    
    for j=1:3
        
        patch(get(h(j),'XData'),get(h(j),'YData'),color(j), 'FaceColor', 'none', 'EdgeColor', color(j), 'LineWidth',3);
        
    end
    
    c = get(gca, 'Children');
    
    [hleg1, hobj1, ~, ~] = legend(c(1:3), 'Between State Comparison', 'Within State Comparison', 'Across Subject Comparison', 'Location', 'NorthEast');
    hleg1.FontSize = 18;
    x = findobj(hobj1,'type','text');
    set(x,'FontSize',18);
    ylim([0 1]);
    
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.5, 0.5, 0.55]);
    
    
    filename = ['/AllSubjects_AllTasks_SpatialOverlap_' num2str(thresholds(v)) '_Percent.jpg'];
    
    print(gcf,[outputdir filename],'-dpng','-r300');
    
    
    close gcf
    
    
    
end




function outmat = binarize_vectors(evenvector,oddvector)

% binarize vectors of unique ID maps

outmat = [];

for vertices = 1:length(oddvector)
    
    if oddvector(vertices) > 0 && evenvector(vertices) > 0
        
        outmat = [outmat;1 1];
        
    elseif oddvector(vertices) > 0
        
        outmat = [outmat;1 0];
        
    elseif evenvector(vertices) > 0
        
        outmat = [outmat;0 1];
        
    end
end


end



    