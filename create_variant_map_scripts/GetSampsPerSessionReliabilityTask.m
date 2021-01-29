function [sampspersessionmem, sampspersessionmixed, sampspersessionmotor, memmissing, mixedmissing, motormissing] = GetSampsPerSessionReliabilityTask(sub,time,TR,NumSampsFiles_path)

%% Samples the desired amount of data from odd-numbered sessions equally across tasks starting with the first session
%
% This function calculates the amount of task data to sample from each
% session to both attempt to equally sample data from each task and sample
% the correct number of sample points of data to equal X minutes. It
% attempts to equally sample the amount of data across all 3 tasks to reach
% the target amount of data. Each task is sampled starting with the first
% session. When data in this session runs out, data from the next
% odd-numbered session is used. This continues until either the target
% amount of data per task is hit or no more data exists. If one task has no
% more high-quality data, any remaining data from the other tasks is used
% to try and hit the target.
%
% INPUTS:
% -NumSampsFiles_path: path to structures that contain the number of sample
% points per task and session (see CreateTMaskStruct_MSC.m)
% -sub: a subject ID to calculate data sampling for
% -time: a time (in minutes) to match the length of the concatenated task 
% data to
% -TR: The TR (sample rate) of the data (in seconds), used for matching the
% number of sample points to the desired time
%
% OUTPUTS:
% -sampspersessionmem: a vector of the amount of data points from the
% memory task to sample from each session
% -sampspersessionmixed: a vector of the amount of data points from the
% mixed (semantic and coherence) task to sample from each session
% -sampspersessionmotor: a vector of the amount of data points from the
% motor task to sample from each session
% -memmissing: a vector of numbers that indicates which sessions lacked
% enough "high-quality" memory task data to sample any data from
% -mixedmissing: a vector of numbers that indicates which sessions lacked
% enough "high-quality" mixed (semantic and coherence) task data to sample
% any data from
% -motormissing: a vector of numbers that indicates which sessions lacked
% enough "high-quality" motor task data to sample any data from
%
% Written by BK, edited by DP & ZL (01-2021)
%

%% Initialize Variables

load([NumSampsFiles_path sub '_NumSamps.mat']);

totalsamppts = round(time*(60/TR));  %% Calculate the number of sample points for a time (in minutes) for a given TR
tasksamppts = ceil(totalsamppts/3);  %% Divide that time by the number of tasks, and sample that amount of data consecutively for each task

%% Loops through odd sessions for each task and tries to sample an equal amount of data from each until the target amount is reached

tasksampptstmp = tasksamppts;

ptscountmem = 0;       %% Running sample point count across sessions
ptscountmixed = 0;
ptscountmotor = 0;
enoughptsmem = 0;      %% Toggles whether enough sample points have been taken consecutively from each file
enoughptsmixed = 0;
enoughptsmotor = 0;
maxptsmem = 0;         %% Toggles whether the maximum amount of data for each task has been sampled
maxptsmixed = 0;
maxptsmotor = 0;
sampspersessionmem = zeros(length(SubStruct),1);
sampspersessionmixed = zeros(length(SubStruct),1);
sampspersessionmotor = zeros(length(SubStruct),1);

memmissing = [];
mixedmissing = [];
motormissing = [];

remaindertaskpts = 0;

while true
    
    for u = 1:length(SubStruct)
        
        % record runs with missing data for later
        
        if SubStruct(u).MemExists == 0
            
            memmissing = [memmissing u];
            
        end
        
        if SubStruct(u).MixedExists == 0
            
            mixedmissing = [mixedmissing u];
            
        end
        
        if SubStruct(u).MotorExists == 0
            
            motormissing = [motormissing u];
            
        end
        
        %% loop through all odd sessions until enough data is sampled or no more data exists
        
        if SubStruct(u).OddEven == 1 && (~enoughptsmem && ~maxptsmem)
            
            if SubStruct(u).MemSampPts + ptscountmem >= tasksampptstmp
                
                sampspersessionmem(u) = sampspersessionmem(u) + (tasksampptstmp - ptscountmem);
                ptscountmem = sampspersessionmem(u) + ptscountmem;
                enoughptsmem = 1;
                
            elseif SubStruct(u).MemSampPts + ptscountmem < tasksampptstmp
                
                sampspersessionmem(u) = sampspersessionmem(u) + SubStruct(u).MemSampPts;
                ptscountmem = SubStruct(u).MemSampPts + ptscountmem;
                
            end
        end
        
        if u == length(SubStruct) && ~enoughptsmem
            
            maxptsmem = 1;
            
        end
        
        if SubStruct(u).OddEven == 1 && (~enoughptsmixed && ~maxptsmixed)
            
            if SubStruct(u).MixedSampPts + ptscountmixed >= tasksampptstmp
                
                sampspersessionmixed(u) = sampspersessionmixed(u) + (tasksampptstmp - ptscountmixed);
                ptscountmixed = sampspersessionmixed(u) + ptscountmixed;
                enoughptsmixed = 1;
                
            elseif SubStruct(u).MixedSampPts + ptscountmixed < tasksampptstmp
                
                sampspersessionmixed(u) = sampspersessionmixed(u) + SubStruct(u).MixedSampPts;
                ptscountmixed = SubStruct(u).MixedSampPts + ptscountmixed;
                
            end
        end
        
        if u == length(SubStruct) && ~enoughptsmixed
            
            maxptsmixed = 1;
            
        end
        
        if SubStruct(u).OddEven == 1 && (~enoughptsmotor && ~maxptsmotor)
            
            if SubStruct(u).MotorSampPts + ptscountmotor >= tasksampptstmp
                
                sampspersessionmotor(u) = sampspersessionmotor(u) + (tasksampptstmp - ptscountmotor);
                ptscountmotor = sampspersessionmotor(u) + ptscountmotor;
                enoughptsmotor = 1;
                
            elseif SubStruct(u).MotorSampPts + ptscountmotor < tasksampptstmp
                
                sampspersessionmotor(u) = sampspersessionmotor(u) + SubStruct(u).MotorSampPts;
                ptscountmotor = SubStruct(u).MotorSampPts + ptscountmotor;
                
            end
        end
        
        if u == length(SubStruct) && ~enoughptsmotor
            
            maxptsmotor = 1;
            
        end
        
    end
    
    % if the target amount of points are sampled or all the data points for
    % a subject have been sampled, end the loop
    
    if (ptscountmem + ptscountmixed + ptscountmotor) >= totalsamppts || (maxptsmem && maxptsmixed && maxptsmotor)
        
        break
        
    else
        
        % if a task runs out of sample points, try to sample the remaining
        % needed points from the other tasks

        if maxptsmem && maxptsmixed
            
            remaindertaskpts = totalsamppts - (ptscountmem + ptscountmixed + tasksamppts);
            
            enoughptsmotor = 0;
            ptscountmotor = 0;
            sampspersessionmotor = zeros(length(SubStruct),1);
            
        elseif maxptsmem && maxptsmotor
            
            remaindertaskpts = totalsamppts - (ptscountmem + ptscountmotor + tasksamppts);
            
            enoughptsmixed = 0;
            ptscountmixed = 0;
            sampspersessionmixed = zeros(length(SubStruct),1);
            
        elseif maxptsmixed && maxptsmotor
            
            remaindertaskpts = totalsamppts - (ptscountmixed + ptscountmotor + tasksamppts);
            
            enoughptsmem = 0;
            ptscountmem = 0;
            sampspersessionmem = zeros(length(SubStruct),1);
            
        elseif maxptsmem
            
            remaindertaskpts = ceil((tasksamppts - ptscountmem)/2);
            
            enoughptsmixed = 0;
            ptscountmixed = 0;
            sampspersessionmixed = zeros(length(SubStruct),1);
            enoughptsmotor = 0;
            ptscountmotor = 0;
            sampspersessionmotor = zeros(length(SubStruct),1);
            
        elseif maxptsmixed
            
            remaindertaskpts = ceil((tasksamppts - ptscountmixed)/2);
            
            enoughptsmem = 0;
            ptscountmem = 0;
            sampspersessionmem = zeros(length(SubStruct),1);
            enoughptsmotor = 0;
            ptscountmotor = 0;
            sampspersessionmotor = zeros(length(SubStruct),1);
            
        elseif maxptsmotor
            
            remaindertaskpts = ceil((tasksamppts - ptscountmotor)/2);
            
            enoughptsmem = 0;
            ptscountmem = 0;
            sampspersessionmem = zeros(length(SubStruct),1);
            enoughptsmixed = 0;
            ptscountmixed = 0;
            sampspersessionmixed = zeros(length(SubStruct),1);
            
        end
    end
    
    tasksampptstmp = tasksamppts + remaindertaskpts;
    
end

% delete missing sessions for each task from data sampling

if ~isempty(memmissing)
    
    sampspersessionmem(memmissing) = [];
    
    
end

if ~isempty(mixedmissing)
    
    sampspersessionmixed(mixedmissing) = [];
    
end

if ~isempty(motormissing)
    
    sampspersessionmotor(motormissing) = [];
    
end

% if not enough data, do not create a variant map. else, create one with
% the sampled data

if maxptsmem && maxptsmixed && maxptsmotor
    
    fprintf('Subject %s does not have %i minutes of total task data. There are %i data points for the memory task, %i data points for the mixed task, and %i data points for the motor task (%.2f minutes total): %s\n', sub, time, ptscountmem, ptscountmixed, ptscountmotor, ((ptscountmem + ptscountmixed + ptscountmotor)*TR)/60, datestr(now));
    
    sampspersessionmem = [];
    sampspersessionmixed = [];
    sampspersessionmotor = [];
    
else
    
    fprintf('Subject %s has enough data for %i minutes of total task data. There are %i data points for the memory task, %i data points for the mixed task, and %i data points for the motor task (%.2f minutes total): %s\n', sub, time, ptscountmem, ptscountmixed, ptscountmotor, ((ptscountmem + ptscountmixed + ptscountmotor)*TR)/60, datestr(now));
    
end