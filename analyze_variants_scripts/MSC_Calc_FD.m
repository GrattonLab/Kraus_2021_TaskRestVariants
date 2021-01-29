
clear all

%% Calculate FD for MSC participants taking into account only points which fall within each task period
%
% This script requires QC structures which contain data on the FD values
% for each frame of each run, task, and session within each subject. In
% addition, for the task data, temporal masks demarcating which frames in
% each run occur during task performance and not fixation periods are also
% used to calculate FD only during task performance. The script
% concatenates FD across sessions and calculates the summary values. These
% values are output into a structure in the first part of the script. In 
% the second part, the data in this structure is used to summarize the FD
% values for each task and subject.
%
% INPUTS:
% -outdir: an output directory for the created structure
% -QC_Task_Root_Dir: root of the filepath where MSC task data tmasks and FD
% values are located
% -QC_Task_Root_Dir: root of the filepath where MSC rest data tmasks and FD
% values are located
% -rest_folder_stem: part of the filepath where the MSC rest data is
% located
% -subs: a cell array of subject IDs to calculate FD for
% -tasks: a cell array of tasks to calculate FD for
%
% OUTPUTS:
% -MSC_FD_Summary: A structure for each subject containing the mean FD and
% the percentage of bad frames (FD < .2 mm) per usable run of data. The
% second part of the script summarizes these metrics for each subject
%
% Written by BK (01-2021)
%


%% Initialize variables

% set up directories and output

QC_Task_Root_Dir = '/your/path/here/';  %% root path to task QC files
QC_Rest_Root_Dir = '/your/path/here/';  %% root path to rest QC files
rest_folder_stem = '/your/path/here/';  %% location of rest data and temporal masks
outputdir = '/your/path/here/';  %% output directory for FD summary structure

subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed','rest'};

% set up output structure for summary metrics

MSC_FD_Summary = struct('Subject',{},'Task',{},'Mean_FD',{},'Percent_Bad_FD',{},'Mean_FD_Condindices',{},'Percent_Bad_FD_Condindices',{});

%% Loop through subjects, get mean FD and % Bad Frames (FD < .2 mm)

for sub = 1:numel(subs)  %% loop through each subject
    
    for task = 1:numel(tasks)  %% loop through each task
        
        if ~strcmp(tasks{task}, 'rest')  %% if current task is not rest
            
            % Load task tmasks from correct directories (FD filtered for
            % subs MSC03 and MSC10)
            
            if strcmp(subs{sub}, 'MSC03') || strcmp(subs{sub}, 'MSC10')
                
                load([QC_Task_Root_Dir 'FCProc_' subs{sub} '_' tasks{task} '_pass2_FDfilt/QC.mat']);
                
                load([QC_Task_Root_Dir 'FCProc_' subs{sub} '_' tasks{task} '_pass2_FDfilt/condindices.mat']);
                
            else
        
                load([QC_Task_Root_Dir 'FCProc_' subs{sub} '_' tasks{task} '_pass2/QC.mat']);
                
                load([QC_Task_Root_Dir 'FCProc_' subs{sub} '_' tasks{task} '_pass2/condindices.mat']);
                
            end
                
        else  %% else, load rest tmasks
            
            load([QC_Rest_Root_Dir subs{sub} rest_folder_stem 'QC.mat'])
            
        end
                
    
        % concatenate data across runs
        
        FD_concat = [];
        
        for runs = 1:length(QC)
            
            FD_concat = [FD_concat; QC(runs).FD];
            
        end
        
        % calculate FD for each task
        
        PercentFD_condindices = [];
        MeanFD_condindices = [];
        
        if ~strcmp(tasks{task}, 'rest')  %% if current task is not rest
            
            runcount = 0;  %% count number of runs for missing task runs

            for runs = 1:length(FDtot)  %% for each run
                
                if strcmp(tasks{task}, 'mem')
                    
                    if sum(TIndTot(runs).AllMem) > 0  %% if any high-quality data exists for a session and task
                        
                        runcount=runcount+1;  %% add 1 to run count
                        FD_cond = QC(runcount).FD(logical(TIndTot(runs).AllMem));  %% get FD for frames that fall inside the run
                        PercentFD_condindices = [PercentFD_condindices; length(find(FD_cond > .2))/length(FD_cond)];  %% calculate the percentage of frames within each task session where FD exceeds .2 mm
                        
                    end
                    
                    if ~isnan(FDtot(runs).AllMem)  %% if FD exists for current subject, session, and task

                        MeanFD_condindices = [MeanFD_condindices; FDtot(runs).AllMem];  %% add mean FD during task to calculation
                        
                    end
                    
              	elseif strcmp(tasks{task}, 'mixed')
                    
                    if sum(cell2mat(struct2cell(TIndTot(runs)))) > 0
                        
                        runcount=runcount+1;  %% add 1 to run count
                        FD_cond = QC(runcount).FD(logical(sum([TIndTot(runs).AllGlass TIndTot(runs).AllSemantic],2)));  %% get FD for frames that fall inside the run
                        PercentFD_condindices = [PercentFD_condindices; length(find(FD_cond > .2))/length(FD_cond)];  %% calculate the percentage of frames within each task session where FD exceeds .2 mm
                        
                    end
                    
                    if ~isnan(mean(cell2mat(struct2cell(FDtot(runs)))))
                    
                        MeanFD_condindices = [MeanFD_condindices; mean([FDtot(runs).AllGlass FDtot(runs).AllSemantic])];  %% add mean FD during task to calculation
                    
                    end
                        
                elseif strcmp(tasks{task}, 'motor')
                    
                    if sum(TIndTot(runs).AllMotor) > 0
                        
                        runcount=runcount+1;  %% add 1 to run count
                        FD_cond = QC(runcount).FD(logical(TIndTot(runs).AllMotor));  %% get FD for frames that fall inside the run
                        PercentFD_condindices = [PercentFD_condindices; length(find(FD_cond > .2))/length(FD_cond)];  %% calculate the percentage of frames within each task session where FD exceeds .2 mm
                        
                    end
                    
                    if ~isnan(FDtot(runs).AllMotor)
                    
                        MeanFD_condindices = [MeanFD_condindices; FDtot(runs).AllMotor];  %% add mean FD during task to calculation

                    end
            
                end
            end
        end
        
        % add data from each task to the output structure
        
        StructAdd = size(MSC_FD_Summary,2)+1;
        
        MSC_FD_Summary(StructAdd).Subject = subs{sub};  %% sub ID
        MSC_FD_Summary(StructAdd).Task = tasks{task};  %% task name
        MSC_FD_Summary(StructAdd).Mean_FD = mean(FD_concat);  %% Mean of concatenated FD (regardless of task indices)
        MSC_FD_Summary(StructAdd).Percent_Bad_FD = length(find(FD_concat > .2))/length(FD_concat);  %% percent of all bad frames (regardless of task indices)
        
        if ~strcmp(tasks{task}, 'rest')  %% if current task is not rest
            
            MSC_FD_Summary(StructAdd).Mean_FD_Condindices = mean(MeanFD_condindices);  %% mean of concatenated FD (within task indices)
            MSC_FD_Summary(StructAdd).Percent_Bad_FD_Condindices = mean(PercentFD_condindices);  %% percent of all bad frames (within task indices)
            
        else % else, if task is rest
            
            MSC_FD_Summary(StructAdd).Mean_FD_Condindices = mean(FD_concat);  %% mean of concatenated FD (regardless of task indices)
            MSC_FD_Summary(StructAdd).Percent_Bad_FD_Condindices = length(find(FD_concat > .2))/length(FD_concat);  %% percent of all bad frames (regardless of task indices)
            
        end
        
    end
end

% write output structure

save([outputdir 'MSC_FD_Summary_Struct.mat'], 'MSC_FD_Summary');

%% Calculate FD summary statistics

clear all

%% Initialize Variables

FDdir = outputdir;  %% load structure from output directory

load([FDdir 'MSC_FD_Summary_Struct.mat'])

restFD = [];
motorFD = [];
mixedFD = [];
memoryFD = [];

restframes = [];
motorframes = [];
mixedframes = [];
memoryframes = [];

motorFDCond = [];
mixedFDCond = [];
memoryFDCond = [];

motorframesCond = [];
mixedframesCond = [];
memoryframesCond = [];

%% Loop through tasks, get summary statistics for each subject and task

for tasks = 1:size(MSC_FD_Summary,2)
    
    if strcmp(MSC_FD_Summary(tasks).Task, 'rest')
        
        restFD = [restFD; MSC_FD_Summary(tasks).Mean_FD];
        restframes = [restframes; MSC_FD_Summary(tasks).Percent_Bad_FD];
        
    elseif strcmp(MSC_FD_Summary(tasks).Task, 'motor')
        
        motorFD = [motorFD; MSC_FD_Summary(tasks).Mean_FD];
        motorframes = [motorframes; MSC_FD_Summary(tasks).Percent_Bad_FD];
        motorFDCond = [motorFDCond; MSC_FD_Summary(tasks).Mean_FD_Condindices];
        motorframesCond = [motorframesCond; MSC_FD_Summary(tasks).Percent_Bad_FD_Condindices];
        
    elseif strcmp(MSC_FD_Summary(tasks).Task, 'mixed')
        
        mixedFD = [mixedFD; MSC_FD_Summary(tasks).Mean_FD];
        mixedframes = [mixedframes; MSC_FD_Summary(tasks).Percent_Bad_FD];
      	mixedFDCond = [mixedFDCond; MSC_FD_Summary(tasks).Mean_FD_Condindices];
        mixedframesCond = [mixedframesCond; MSC_FD_Summary(tasks).Percent_Bad_FD_Condindices];
        
    elseif strcmp(MSC_FD_Summary(tasks).Task, 'mem')
        
        memoryFD = [memoryFD; MSC_FD_Summary(tasks).Mean_FD];
        memoryframes = [memoryframes; MSC_FD_Summary(tasks).Percent_Bad_FD];
        memoryFDCond = [memoryFDCond; MSC_FD_Summary(tasks).Mean_FD_Condindices];
        memoryframesCond = [memoryframesCond; MSC_FD_Summary(tasks).Percent_Bad_FD_Condindices];
        
    end
end

%% Summary statistics

restFDsum = mean(restFD);
restframessum = mean(restframes);
motorFDsum = mean(motorFD);
motorframessum = mean(motorframes);
mixedFDsum = mean(mixedFD);
mixedframessum = mean(mixedframes);
memoryFDsum = mean(memoryFD);
memoryframessum = mean(memoryframes);
motorFDsumCond = mean(motorFDCond);
motorframessumCond = mean(motorframesCond);
mixedFDsumCond = mean(mixedFDCond);
mixedframessumCond = mean(mixedframesCond);
memoryFDsumCond = mean(memoryFDCond);
memoryframessumCond = mean(memoryframesCond);
alltasksFDsum = mean([motorFDsum mixedFDsum memoryFDsum]);
alltasksframessum = mean([motorframessum mixedframessum memoryframessum]);
alltasksFDsumCond = mean([motorFDsumCond mixedFDsumCond memoryFDsumCond]);
alltasksframessumCond = mean([motorframessumCond mixedframessumCond memoryframessumCond]);
