
clear all

%% Plot data for spatial overlap within states, between states, and across subjects within predetermined bins
%
% This script loads data from variant maps and thresholds them at 
% predetermined values. Then it finds the spatial location overlap for
% for variants within states (rest vs. rest and task vs. task), between
% states (rest vs. task), and across subjects (e.g., MSC01 task vs. MSC02 
% rest). Spatial location overlap is calculated according to binarized maps
% (a 1 represents a vertex where a variant exists, and 0 a vertex where a
% variant does not exist). A Dice-Sorenson correlation is then performed on
% these maps to calculate spatial overlap. The values for the within state,
% between state, and across subject comparisons are then plotted with
% standard error bars for each bin. The spatial maps for each bin can also
% be output as CIFTIs.
%
% INPUTS:
% 
% -outputdir: the output directory for the plots of spatial location 
% overlap in each bin
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
% -absthresh: toggles whether to use an absolute correlation (r)
% value to define variants (set to 1), otherwise thresholds in the bins
% variable are interpreted as the xth lowest percentage of correlation
% values (set to 0)
% -SaveSpatialMaps: toggles whether to save out spatial maps (CIFTIs) 
% corresponding to each bin for each file and subject (set to 1), otherwise
% does not save out spatial maps (set to 0)
% -bins: an x by 2 matrix where each row corresponds to the lower and upper
% thresholds (in the first and second columns, respectively) for each bin
% in which to compute a spatial correlation
%
% -task_files_even: reads path to a text file containing the paths to each
% of the variant split-half files for task data in even-numbered sessions.
% The format is (pathtofile subjectID) and the file is space-delimited (see
% below for more formatting details)
% -task_files_odd: reads path to a text file containing the paths to each
% of the variant split-half files for task data in odd-numbered sessions.
% The format is (pathtofile subjectID) and the file is space-delimited (see
% below for more formatting details)
% -rest_files_even: reads path to a text file containing the paths to each
% of the variant split-half files for rest data in even-numbered sessions.
% The format is (pathtofile subjectID) and the file is space-delimited (see
% below for more formatting details)
% -rest_files_odd: reads path to a text file containing the paths to each
% of the variant split-half files for rest data in odd-numbered sessions.
% The format is (pathtofile subjectID) and the file is space-delimited (see
% below for more formatting details)
%
% OUTPUTS:
%
% -Binarized Maps: A CIFTI file (numvertices x 1) with binarized (0s and 
% 1s) values. The 1s correspond to vertices that fell within each bin at
% each threshold, and the 0s are vertices which did not fall within a given
% bin
% -plots: creates a plot for the spatial location overlap of variants
% within states, between states, and across subjects by each subject as
% well as combined
%
% Written by BK (01-2021)
%

%% Initialize Variables

outputdir = '/your/path/here/';  %% output directory for plot of spatial overlap
SNRpath = '/your/path/here/';  %% path to SNR maps
SNR_fstring = 'SNRname.dtseries.nii';  %% filename and extension for SNR map
templatepath = '/your/path/here/template_cifti.dtseries.nii';  %% Path where template cifti is located

SNRExclude = 1;        %% Toggles whether to exclude low SNR areas
AbsoluteThresholds = 0;        %% Toggles whether threshold is set by percentile or absolute value
SaveSpatialMaps = 1;  %% Toggles whether to save out spatial maps of each percentile bin


bins = [0 2.5;      %% Percentage bins
           2.5 5;
           5 7.5;
           7.5 10;
           10 12.5;
           12.5 15;
           15 17.5;
           17.5 20;
           20 22.5;
           22.5 25
           25 27.5
           27.5 30
           30 32.5
           32.5 35
           35 37.5
           37.5 40
           40 42.5
           42.5 45
           45 47.5
           47.5 50
           50 52.5
           52.5 55
           55 57.5
           57.5 60
           60 62.5
           62.5 65
           65 67.5
           67.5 70
           70 72.5
           72.5 75
           75 77.5
           77.5 80
           80 82.5
           82.5 85
           85 87.5
           87.5 90
           90 92.5
           92.5 95
           95 97.5
           97.5 100];
       
if SaveSpatialMaps  %% create output directory for spatial maps
     
    if ~exist([outputdir 'Binned_Spatial_Maps'], 'dir')
        mkdir([outputdir 'Binned_Spatial_Maps']);
    end
    outputdirmaps = [outputdir 'Binned_Spatial_Maps/'];

end


%% Loop through data, load variant map and threshold for each bin


DiceCorrsTaskRest = [];  % Spatial overlap between state
DiceCorrsTaskTask = [];  % Spatial overlap within rest
DiceCorrsRestRest = [];  % Spatial overlap within task

alltaskfilestaskeven = {};  % Store data from each condition for permutation across subjects
alltaskfilestaskodd = {};
alltaskfilesresteven = {};
alltaskfilesrestodd = {};

for g = 1:size(bins,1)
    
    % load variant maps for analyses using a space-delimited text file in 
    % the format: pathtofile subID
    % e.g. filepath/MSC01.dtseries.nii MSC01
    % the order of the data files for all of the subjects should be the
    % same in all text files
   
    [task_files_even, ~, ~] = textread('/your/path/here/MSC_alltask_even_varmaps.txt','%s%s%s');
    
    [rest_files_even, ~, ~] = textread('/your/path/here/MSC_rest_even_varmaps.txt','%s%s%s');
    
    [task_files_odd, subjects, ~] = textread('/your/path/here/MSC_alltask_odd_varmaps.txt','%s%s%s');
    
    [rest_files_odd, ~, ~] = textread('/your/path/here/MSC_rest_odd_varmaps.txt','%s%s%s');

    subs = unique(subjects);

    DiceCorrsTaskRestTemp = [];  %% Store data for each subject temporarily
    DiceCorrsTaskTaskTemp = [];
    DiceCorrsRestRestTemp = [];
    
    alltaskfilestaskeventemp = [];
    alltaskfilestaskoddtemp = [];
    alltaskfilesresteventemp = [];
    alltaskfilesrestoddtemp = [];

    for s = 1:numel(subs)
    
        subject = subs{s};

        cifti_rest_even = ft_read_cifti_mod(rest_files_even{s});  %% Read data from filepaths for each subejct for each split-half
        cifti_task_even = ft_read_cifti_mod(task_files_even{s});
        cifti_rest_odd = ft_read_cifti_mod(rest_files_odd{s});
        cifti_task_odd = ft_read_cifti_mod(task_files_odd{s});

        if SNRExclude  %% Exclude low SNR areas
            
            SNRmap = ft_read_cifti_mod([SNRpath subject SNR_fstring]);

            SNRmap.data = SNRmap.data(1:59412,:);
            
            LowSNR = find(SNRmap.data < 750);

            cifti_rest_even.data(LowSNR) = NaN;
            cifti_task_even.data(LowSNR) = NaN;
            cifti_rest_odd.data(LowSNR) = NaN;
            cifti_task_odd.data(LowSNR) = NaN;
        
        end
        
        if AbsoluteThresholds  %% get thresholds for each bin
            
            cifti_rest_even_thresh = bins(g,:);
            cifti_task_even_thresh = bins(g,:);
            cifti_rest_odd_thresh = bins(g,:);
            cifti_task_odd_thresh = bins(g,:);
            
        else
            
            cifti_rest_even_thresh = prctile(cifti_rest_even.data,bins(g,:));
            cifti_task_even_thresh = prctile(cifti_task_even.data,bins(g,:));
            cifti_rest_odd_thresh = prctile(cifti_rest_odd.data,bins(g,:));
            cifti_task_odd_thresh = prctile(cifti_task_odd.data,bins(g,:));
            
        end
        
        % find indices in each variant map that are outside each bin
        
        cifti_rest_even_remove = find(cifti_rest_even.data < cifti_rest_even_thresh(1) | cifti_rest_even.data >= cifti_rest_even_thresh(2));
        cifti_task_even_remove = find(cifti_task_even.data < cifti_task_even_thresh(1) | cifti_task_even.data >= cifti_task_even_thresh(2));
        cifti_rest_odd_remove = find(cifti_rest_odd.data < cifti_rest_odd_thresh(1) | cifti_rest_odd.data >= cifti_rest_odd_thresh(2));
        cifti_task_odd_remove = find(cifti_task_odd.data < cifti_task_odd_thresh(1) | cifti_task_odd.data >= cifti_task_odd_thresh(2));
        
        % set these indices that don't fall within current bin to NaN
        
        cifti_rest_even.data(cifti_rest_even_remove) = NaN;
        cifti_task_even.data(cifti_task_even_remove) = NaN;
        cifti_rest_odd.data(cifti_rest_odd_remove) = NaN;
        cifti_task_odd.data(cifti_task_odd_remove) = NaN;

        % loop through data, create binarized vectors based off of non-NaN
        % values in data
        
        dcorrdatataskrest = create_nan_vecs(cifti_rest_even.data,cifti_task_odd.data);    %% Task-Rest Comparison
        dcorrdatarestrest = create_nan_vecs(cifti_rest_odd.data,cifti_rest_even.data);    %% Rest-Rest Comparison
        dcorrdatatasktask = create_nan_vecs(cifti_task_odd.data,cifti_task_even.data);    %% Task-Task Comparison
 
        % exclude empty vectors from spatial correlation (otherwise throws errror)
        
        if isempty(dcorrdatataskrest) && isempty(dcorrdatarestrest) && isempty(dcorrdatatasktask)
            
            dctaskrest = 0;
            dcrestrest = 0;
            dctasktask = 0;
            
        elseif isempty(dcorrdatataskrest) && isempty(dcorrdatarestrest)
            
            dctaskrest = 0;
            dcrestrest = 0;
            dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
            
        elseif isempty(dcorrdatataskrest) && isempty(dcorrdatatasktask)
            
            dctaskrest = 0;
            dctasktask = 0;
            dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
            
        elseif isempty(dcorrdatarestrest) && isempty(dcorrdatatasktask)
            
            dcrestrest = 0;
            dctasktask = 0;
            dctaskrest = dice_coefficient_mod(dcorrdatataskrest(:,1),dcorrdatataskrest(:,2));
            
        elseif isempty(dcorrdatarestrest)
            
            dcrestrest = 0;
            dctaskrest = dice_coefficient_mod(dcorrdatataskrest(:,1),dcorrdatataskrest(:,2));
            dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
            
        elseif isempty(dcorrdatatasktask)
            
            dctasktask = 0;
            dctaskrest = dice_coefficient_mod(dcorrdatataskrest(:,1),dcorrdatataskrest(:,2));
            dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
            
        elseif isempty(dcorrdatataskrest)
            
            dctaskrest = 0;
            dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
            dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
            
        else
            
            dctaskrest = dice_coefficient_mod(dcorrdatataskrest(:,1),dcorrdatataskrest(:,2));
            dcrestrest = dice_coefficient_mod(dcorrdatarestrest(:,1),dcorrdatarestrest(:,2));
            dctasktask = dice_coefficient_mod(dcorrdatatasktask(:,1),dcorrdatatasktask(:,2));
            
        end
        
        % Save correlation values and store data for across subject
        % permutations
        
        DiceCorrsTaskRestTemp = [DiceCorrsTaskRestTemp; dctaskrest];
        DiceCorrsTaskTaskTemp = [DiceCorrsTaskTaskTemp; dctasktask];
        DiceCorrsRestRestTemp = [DiceCorrsRestRestTemp; dcrestrest];
           
  
        alltaskfilestaskeventemp = [alltaskfilestaskeventemp; cifti_task_even.data'];
        alltaskfilestaskoddtemp = [alltaskfilestaskoddtemp; cifti_task_odd.data'];
        alltaskfilesresteventemp = [alltaskfilesresteventemp; cifti_rest_even.data'];
        alltaskfilesrestoddtemp = [alltaskfilesrestoddtemp; cifti_rest_odd.data'];

        
        % Write spatial maps for each decile
        
        if SaveSpatialMaps
        
            template = ft_read_cifti_mod(templatepath);
            template.data = [];
            template2 = template;
            template3 = template;
            template4 = template;
            
            cifti_rest_even.data(~isnan(cifti_rest_even.data)) = 1;  %% Binarize data
            cifti_rest_odd.data(~isnan(cifti_rest_odd.data)) = 1;
            cifti_task_even.data(~isnan(cifti_task_even.data)) = 1;
            cifti_task_odd.data(~isnan(cifti_task_odd.data)) = 1;
            
            template.data = cifti_rest_even.data;  % replace template data with data from each bin
            template2.data = cifti_rest_odd.data;
            template3.data = cifti_task_even.data;
            template4.data = cifti_task_odd.data;
            
            resteven_fname = create_filename(subject,bins(g,:),'Rest','Even',AbsoluteThresholds);  % generate filename for each map
            restodd_fname = create_filename(subject,bins(g,:),'Rest','Odd',AbsoluteThresholds);
            taskeven_fname = create_filename(subject,bins(g,:),'Task','Even',AbsoluteThresholds);
            taskodd_fname = create_filename(subject,bins(g,:),'Task','Odd',AbsoluteThresholds);
            
            
            ft_write_cifti_mod([outputdirmaps resteven_fname],template);  % write out each spatial map
            ft_write_cifti_mod([outputdirmaps restodd_fname],template2);
            ft_write_cifti_mod([outputdirmaps taskeven_fname],template3);
            ft_write_cifti_mod([outputdirmaps taskodd_fname],template4);
            
            clear template
            clear template2
            clear template3
            clear template4
        
        end
        
    end
    
    DiceCorrsTaskRest = [DiceCorrsTaskRest DiceCorrsTaskRestTemp];
    DiceCorrsTaskTask = [DiceCorrsTaskTask DiceCorrsTaskTaskTemp];
    DiceCorrsRestRest = [DiceCorrsRestRest DiceCorrsRestRestTemp];
        
      
    alltaskfilestaskeven = [alltaskfilestaskeven alltaskfilestaskeventemp];
    alltaskfilestaskodd = [alltaskfilestaskodd alltaskfilestaskoddtemp];
    alltaskfilesresteven = [alltaskfilesresteven alltaskfilesresteventemp];
    alltaskfilesrestodd = [alltaskfilesrestodd alltaskfilesrestoddtemp];


end




%% Across subject permutation

SubjectDiceCorrsFinal = [];

% loop across subjects and calculate the spatial overlap of each bin

for a = 1:size(bins,1)

    SubjectDiceCorrs = [];
    
    for j = 1:length(subs)
        
        subject = subs{j};
        
        % for each split-half pair, get data for each threshold for each subject

        alltasktasksubject = alltaskfilestaskeven{:,a}(j,:);
        alltaskrestsubject = alltaskfilesresteven{:,a}(j,:);
        
        % leave that subject out of the comparison data
        
        alltasktaskcomp = alltaskfilestaskodd{:,a};
        alltasktaskcomp(j,:) = [];
        alltaskrestcomp = alltaskfilesrestodd{:,a};
        alltaskrestcomp(j,:) = [];
        
        % loop through all other subjects and compare them to the
        % comparison subject's data

        for d = 1:size(alltasktaskcomp,1)

            SubjectDiceCorrstemp = [];
            
            % compare rest split-halves to all other subjects'
            % corresponding task split-halves
   
            dcorrdatataskcomp = create_nan_vecs(alltaskrestsubject,alltasktaskcomp(d,:));
            dcorrdatarestcomp = create_nan_vecs(alltasktasksubject,alltaskrestcomp(d,:));
            
            %  exclude empty vectors from spatial correlation (otherwise throws errror)
   
           	if isempty(dcorrdatataskcomp) && isempty(dcorrdatarestcomp)
                        
               	dctaskcomp = 0;
               	dcrestcomp = 0;
                        
           	elseif isempty(dcorrdatataskcomp)
                        
               	dctaskcomp = 0;
               	dcrestcomp = dice_coefficient_mod(dcorrdatarestcomp(:,1),dcorrdatarestcomp(:,2));
                        
           	elseif isempty(dcorrdatarestcomp)
                        
               	dcrestcomp = 0;
              	dctaskcomp = dice_coefficient_mod(dcorrdatataskcomp(:,1),dcorrdatataskcomp(:,2));
                        
            else
                        
               	dctaskcomp = dice_coefficient_mod(dcorrdatataskcomp(:,1),dcorrdatataskcomp(:,2));
              	dcrestcomp = dice_coefficient_mod(dcorrdatarestcomp(:,1),dcorrdatarestcomp(:,2));
                        
            end
        
          	SubjectDiceCorrstemp = [SubjectDiceCorrstemp dctaskcomp dcrestcomp];

        end
        
        SubjectDiceCorrs = [SubjectDiceCorrs; SubjectDiceCorrstemp];

    end
    
    % average across iterations for final across subject permutation values
    
    SubjectDiceCorrsFinal = [SubjectDiceCorrsFinal mean(SubjectDiceCorrs,2)];

end

%% Plot Data

% Plot data for spatial overlap of each bin

if AbsoluteThresholds
    
    plotthresholds = [.1:.025:.725];
    
else
    
    plotthresholds = [2.5:2.5:100];
    
end



withinstate = [DiceCorrsRestRest; DiceCorrsTaskTask];
betweenstate = DiceCorrsTaskRest;
acrosssubs = SubjectDiceCorrsFinal;

shadedErrorBar(plotthresholds,mean(withinstate,1),std(withinstate,1)/sqrt(size(withinstate,1)),'g');
hold on
shadedErrorBar(plotthresholds,mean(betweenstate,1),std(betweenstate,1)/sqrt(size(betweenstate,1)),'b');
hold on
shadedErrorBar(plotthresholds,mean(acrosssubs,1),std(acrosssubs,1)/sqrt(size(acrosssubs,1)),'k');
c = get(gca, 'Children');
hleg1 = legend(c(1:3), 'Across Subjects Overlap', 'Between State Overlap', 'Within State Overlap', 'Location', 'North');
hleg1.FontSize = 14;
title('Decile Overlap Plot', 'fontsize',18)
ylabel('Dice Correlations')
xlabel('Similarity to Group Functional Connectivity')
ax = gca;
ax.FontSize = 14;
ylim([0 .75])

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.3]);

print(gcf,[outputdir '/AllSubjects_Bin_Spatial_Overlap_Plot.jpg'],'-dpng','-r300');

close gcf





function binarizevec = create_nan_vecs(ciftidata1,ciftidata2)

% binarize vectors for spatial correlation in each bin

binarizevec = [];

for q = 1:length(ciftidata1)
    
    if ~isnan(ciftidata1(q)) && ~isnan(ciftidata2(q))
        
        binarizevec = [binarizevec;1 1];
        
    elseif ~isnan(ciftidata1(q))
        
        binarizevec = [binarizevec;1 0];
        
    elseif ~isnan(ciftidata2(q))
        
        binarizevec = [binarizevec;0 1];
        
    end
    
end

end


function filename = create_filename(sub,thresholds,task,splithalf,AbsoluteThresholds)

% set output filename for each map

ses = splithalf;
grp = 'WashU120';
desc = 'matched';

if AbsoluteThresholds
    thresh_str = [num2str(thresholds(1)) 'r-' num2str(thresholds(2)) 'r'];
else
    thresh_str = [num2str(thresholds(1)) 'pct-' num2str(thresholds(2)) 'pct'];
end


filename = sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_BinnedMap.dtseries.nii',sub,task,ses,grp,desc,thresh_str);

end
