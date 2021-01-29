
clear all

%% Creates a data structure for each subject which contains the number of high-quality sample points for each task in each session
%
% This script requires temporal masks (tmasks) that indicate which
% timepoints in a given run contain "high-quality" data (i.e., low motion).
% These files should be structured as a binarized (timepoints x 1) vector
% (containing 1s and 0s with 1s representing "high-quality" timepoints. The
% values in this binarized file should correspond to the timepoints in a
% CIFTI file for a given subject/task/session combination which have usable
% data. These masks are then used to create a summary output structure of
% "high-quality" sample points that can be used to equate data across
% subjects and tasks.
%
% INPUTS:
% -outdir: an output directory for the created structure
% -MSCTaskdir: filepath where MSC task data tmasks are located
% -rest_pathroot/rest_folder_stem: parts of the filepath where the MSC rest
% data is located
% -subs: a cell array of subject IDs to calculate "high-quality" sample
% points for
% -tasks: a cell array of tasks to calculate "high-quality" sample
% points for
%
% OUTPUTS:
% -SubStruct: A structure for each subject containing the amount of
% "high-quality" sample points for each task and session, as well as
% whether the session is odd or even numbered for creating split-halves of
% the data.
%
% Written by BK (01-2021)
%

%% Initialize Variables

outdir = '/your/path/here/';  %% specify output directory
MSCTaskdir = '/your/path/here/';  %% specify location of task data
rest_pathroot = '/your/path/here/'; %% specify location of rest data
rest_folder_stem = '/your/path/here/';  %% location of rest data and temporal masks (tmasks)

subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed','rest'};

if ~exist(outdir , 'dir')
mkdir(outdir)
end

%% Loop through subjects, add number of high-quality sample points per session to structure

for i=1:numel(subs)

    % Create Structure to output for each subject
    
    SubStruct = struct('SessionFile',{},'MemExists',[],'MemSampPts',[], 'MixedExists',[],'MixedSampPts',[], 'MotorExists',[],'MotorSampPts',[], 'AllTaskSampPts', [], 'RestExists', [], 'RestSampPts', [], 'OddEven', []);
    
    for j=1:length(tasks)
        
        if ~strcmp(tasks{j},'rest')  %% if task not rest
        
            % Load task tmasks from correct directories (FD filtered for
            % subs MSC03 and MSC10)

            if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')

                MSCcondidir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2_FDfilt/']);

            else

                MSCcondidir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2/']);

            end

            load ([MSCcondidir 'condindices.mat']);
            
        end
        
        % Load rest tmasks (even during tasks for vcid names)
        load([rest_pathroot subs{i} rest_folder_stem 'QC.mat'])
        
        % Loop through sessions (vcids) and tabulate good TRs in each
        % session
        for k=1:size(QC,2)

            SubStruct(k).SessionFile = QC(k).vcnum;  %% add session name to structure
            
            tmask = 0;  % Create empty tmask variable
            
            if strcmp(tasks{j},'motor')  %% select correct tmask from current task
                tmask = TIndFin(k).AllMotor;
            elseif strcmp(tasks{j},'mem')
                tmask = TIndFin(k).AllMem;
            elseif strcmp(tasks{j},'mixed')
                tmask = TIndFin(k).AllGlass + TIndFin(k).AllSemantic;
            elseif strcmp(tasks{j},'rest')
                
                % Load rest tmasks from correct variables (FD filtered for
                % subs MSC03 and MSC10)
                
                if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')
                    tmask = QC(k).filtertmask;
                else
                    tmask = QC(k).tmask;
                end
            end
            
            fprintf('tmask for task %s file has %i good sample points, %s\n', tasks{j}, sum(tmask), datestr(now));
            
            if sum(tmask) == 0      % Skip file if no task task data
                
                if strcmp(tasks{j},'motor')  %% set correct structure field for each task
                    SubStruct(k).MotorExists = 0;  %% mark task as not existing
                    SubStruct(k).MotorSampPts = 0;  %% set number of high-quality sample points in this session to 0
                elseif strcmp(tasks{j},'mem')
                    SubStruct(k).MemExists = 0;
                    SubStruct(k).MemSampPts = 0;
                elseif strcmp(tasks{j},'mixed')
                    SubStruct(k).MixedExists = 0;
                    SubStruct(k).MixedSampPts = 0;
             	elseif strcmp(tasks{j},'rest')
                    SubStruct(k).RestExists = 0;  %% mark rest as not existing
                    SubStruct(k).RestSampPts = 0;  %% set number of high-quality sample points in this session to 0    
                end
                
                fprintf('Vcid %s for subject %s has no usable data for task %s\n', QC(k).vcnum, subs{i}, tasks{j});
                
            else                    % Else get first xth points in resting file corresponding to number of task data points
                
                if strcmp(tasks{j},'motor')    %% set correct structure field for each task
                    SubStruct(k).MotorExists = 1;  %% mark task as existing
                    SubStruct(k).MotorSampPts = sum(tmask);    %% set number of high-quality TRs to number of good sample points
                elseif strcmp(tasks{j},'mem')
                    SubStruct(k).MemExists = 1;
                    SubStruct(k).MemSampPts = sum(tmask);
                elseif strcmp(tasks{j},'mixed')
                    SubStruct(k).MixedExists = 1;
                    SubStruct(k).MixedSampPts = sum(tmask);
                elseif strcmp(tasks{j},'rest')    
                    SubStruct(k).RestExists = 1;  %% mark session as existing
                    SubStruct(k).RestSampPts = sum(tmask);  %% set number of high-quality TRs to number of good sample points
                end
            end
            
        	if mod(k,2) == 1  %% if session is odd numbered
                    
              	SubStruct(k).OddEven = 1;  %% mark session as odd numbered (= 1)
                    
            else
                    
            	SubStruct(k).OddEven = 2;  %% mark session as even numbered (= 2)

           	end
            
        end
    end
    
    for l = 1:size(SubStruct,2)  %% Calculate data points for all tasks
        
        SubStruct(l).AllTaskSampPts = (SubStruct(l).MotorSampPts + SubStruct(l).MemSampPts + SubStruct(l).MixedSampPts);
        
    end
    
    % Save final struct for each subject
    
    save([outdir subs{i} '_NumSamps.mat'], 'SubStruct');
    
end

    