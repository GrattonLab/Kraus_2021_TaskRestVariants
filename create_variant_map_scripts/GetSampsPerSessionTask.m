function sampspersession = GetSampsPerSessionTask(mintasksamps,task,sub,NumSampsFiles_path,SplitHalf,MatchData,ConcatenateTasks)

%% Calculates the amount of data to sample from each session for the current task in order to match data across subjects, tasks, and sessions
%
% This function calculates the amount of data to sample from each session 
% for each task and subject given the input parameters. Using the amount of
% data calculated previously (mintasksamps), the function determines how
% much data must be sampled from each session of each task. The function 
% attempts to first equally sample the required amount of data from each 
% session and task. When this is not possible due a lack of "high-quality"
% data in a session, the remaining sessions are used to equally sample the
% required amount of data points from. This process repeats until either
% all task data for a subject has been sampled, or the desired total amount
% of data is reached. 
%
% INPUTS:
% -mintasksamps: the total target amount of data to sample across sessions
% for each task (calculated previously in script if matching data)
% -NumSampsFiles_path: path to structures that contain the number of sample
% points per task and session (see CreateTMaskStruct_MSC.m)
% -sub: a subject ID to calculate data sampling for
% -task: a task to calculate data sampling for
% -SplitHalf: toggles whether data should be concatenated according to
% whether the sessions are odd/even numbered (set to 1), otherwise all
% session data is combined (set to 0)
% -MatchData: toggles whether the amount of data sampled should be matched
% as well as possible across subjects, tasks, and sessions (set to 1), or
% if all possible data should be sampled (set to 0)
% -ConcatenateTasks: toggles whether task data should be concatenated for
% the output files (set to 1), or if each task should concatenated
% separately across sessions (set to 1)
%
% OUTPUTS:
% -sampspersession: a vector of the amount of data points from the current
% task to sample from each session. Is an empty vector if not matching data
%
% Written by BK, edited by DP & ZL (01-2021)
%

fprintf('Matching sampling of task data across sessions for subject %s and task %s: %s\n', sub, task, datestr(now));

if MatchData && SplitHalf
    %% load and extract data from QC structure
    load ([NumSampsFiles_path sub '_NumSamps.mat']);
    
    memptsodd = [];
    motorptsodd = [];
    mixedptsodd = [];
    memptseven = [];
    motorptseven = [];
    mixedptseven = [];
    for u = 1:length(SubStruct)
        if SubStruct(u).OddEven == 1
            memptsodd = [memptsodd; SubStruct(u).MemSampPts];
            mixedptsodd = [mixedptsodd; SubStruct(u).MixedSampPts];
            motorptsodd = [motorptsodd SubStruct(u).MotorSampPts];
        elseif SubStruct(u).OddEven == 2
            memptseven = [memptseven; SubStruct(u).MemSampPts];
            mixedptseven = [mixedptseven; SubStruct(u).MixedSampPts];
            motorptseven = [motorptseven; SubStruct(u).MotorSampPts];
        end
    end
    
    %% calculate the amount to equally sample from each session
    
    % gets minimum number of sample points and then mean needed per
    % session (10 sessions total, 5 in each split-half)
    if strcmp(sub, 'MSC09')
        minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(memptseven) sum(mixedptseven)]);
    else
        minsampspts = min([sum(memptsodd) sum(mixedptsodd) sum(motorptsodd) sum(memptseven) sum(mixedptseven) sum(motorptseven)]);
    end
    
    meansamppts = floor(mintasksamps/5);
    fprintf('For subject %s, the minimum number of sample points for a task in a split-half is %i. the current sampling target is %.0f points per sesssion per task for concatenation: %s\n', sub, minsampspts, meansamppts/3, datestr(now));
    
    meansampptstempodd = meansamppts;
    meansampptstempeven = meansamppts;
    
    % compensates for lack of motor data in MSC09 by oversampling data from
    % the other 2 tasks
    if strcmp(sub, 'MSC09') && strcmp(task, 'motor') && ConcatenateTasks
        meansampptstempodd = floor(sum(motorptsodd)/5);
        meansampptstempeven = floor(sum(motorptseven)/5);
    elseif strcmp(sub, 'MSC09') && ConcatenateTasks
        meansampptstempodd = meansamppts + round((meansamppts - round(sum(motorptsodd)/5))/2);
        meansampptstempeven = meansamppts + round((meansamppts - round(sum(motorptseven)/5))/2);
    end
    
    
    %% Loop through sessions, account for missing sessions and sessions without enough data for sampling in each task
    
    % get the data for all the other subjects
    remaindertotalodd = 0;
    remaindertotaleven = 0;
    notenoughdata = zeros(length(SubStruct),1);
    
    while true
        remainderodd = 0;     %% update amount of points that still need to be sampled on each iteration
        remaindereven = 0;
        sampspersession = zeros(length(SubStruct),1);
        
        notenoughdataeven = sum(notenoughdata(2:2:end));
        notenoughdataodd = sum(notenoughdata(1:2:end));
        
        meansampptstempodd = (meansampptstempodd + round(remaindertotalodd/(5-notenoughdataodd)));
        meansampptstempeven = (meansampptstempeven + round(remaindertotaleven/(5-notenoughdataeven)));
        
        % loop attempts to sample data equally per session. where this is
        % not possible due to a lack of high-quality data, it oversamples
        % equally across the remaining sessions. this process iteratively
        % repeats until either all of the data is sampled, or the target
        % amount of data (mintasksamps as previously calculated) is reached
        
        for v = 1:length(SubStruct)
            if strcmp(task, 'mem')
                if SubStruct(v).OddEven == 1
                    if SubStruct(v).MemSampPts < meansampptstempodd && notenoughdata(v) == 0
                        notenoughdata(v) = 1;
                        remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MemSampPts));
                        sampspersession(v) = SubStruct(v).MemSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    elseif notenoughdata(v) == 0
                        sampspersession(v) = meansampptstempodd;
                        fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    else
                        sampspersession(v) = SubStruct(v).MemSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    end
                else
                    if SubStruct(v).MemSampPts < meansampptstempeven && notenoughdata(v) == 0
                        notenoughdata(v) = 1;
                        remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MemSampPts));
                        sampspersession(v) = SubStruct(v).MemSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    elseif notenoughdata(v) == 0
                        sampspersession(v) = meansampptstempeven;
                        fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    else
                        sampspersession(v) = SubStruct(v).MemSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    end
                end
            elseif strcmp(task, 'mixed')
                if SubStruct(v).OddEven == 1
                    if SubStruct(v).MixedSampPts < meansampptstempodd && notenoughdata(v) == 0
                        notenoughdata(v) = 1;
                        remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MixedSampPts));
                        sampspersession(v) = SubStruct(v).MixedSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    elseif notenoughdata(v) == 0
                        sampspersession(v) = meansampptstempodd;
                        fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    else
                        sampspersession(v) = SubStruct(v).MixedSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    end
                else
                    if SubStruct(v).MixedSampPts < meansampptstempeven && notenoughdata(v) == 0
                        notenoughdata(v) = 1;
                        remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MixedSampPts));
                        sampspersession(v) = SubStruct(v).MixedSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    elseif notenoughdata(v) == 0
                        sampspersession(v) = meansampptstempeven;
                        fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    else
                        sampspersession(v) = SubStruct(v).MixedSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    end
                end
            elseif strcmp(task, 'motor')
                if SubStruct(v).OddEven == 1
                    if SubStruct(v).MotorSampPts < meansampptstempodd && notenoughdata(v) == 0
                        notenoughdata(v) = 1;
                        remainderodd = (remainderodd + (meansampptstempodd - SubStruct(v).MotorSampPts));
                        sampspersession(v) = SubStruct(v).MotorSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    elseif notenoughdata(v) == 0
                        sampspersession(v) = meansampptstempodd;
                        fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    else
                        sampspersession(v) = SubStruct(v).MotorSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    end
                else
                    if SubStruct(v).MotorSampPts < meansampptstempeven && notenoughdata(v) == 0
                        notenoughdata(v) = 1;
                        remaindereven = (remaindereven + (meansampptstempeven - SubStruct(v).MotorSampPts));
                        sampspersession(v) = SubStruct(v).MotorSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    elseif notenoughdata(v) == 0
                        sampspersession(v) = meansampptstempeven;
                        fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    else
                        sampspersession(v) = SubStruct(v).MotorSampPts;
                        fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                    end
                end
            end
        end
        
        % if needed, samples data from other sessions. stops loop when
        % enough data sampled
        if remainderodd == 0 && remaindereven == 0
            fprintf('Data sampling calculated for subject %s and task %s: %s\n', sub, task, datestr(now));
            break
        elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
            fprintf('Not enough data is present for subject %s and task %s: %s\n', sub, task, datestr(now));
            break
        else
            remaindertotalodd = remainderodd;
            remaindertotaleven = remaindereven;
            fprintf('Remainder of %i calculated for odd sessions for subject %s and task %s: %s\n', remainderodd, sub, task, datestr(now));
            fprintf('Remainder of %i calculated for even sessions for subject %s and task %s: %s\n', remaindereven, sub, task, datestr(now));
        end
    end
    
elseif MatchData
    %% load and extract data from QC structure
    load ([NumSampsFiles_path sub '_QCFile.mat']);
    
    mempts = [];
    motorpts = [];
    mixedpts = [];
    
    % get number of data points from QC files for each task
    for u = 1:length(SubStruct)
        mempts = [mempts; SubStruct(u).MemSampPts];
        mixedpts = [mixedpts; SubStruct(u).MixedSampPts];
        motorpts = [motorpts; SubStruct(u).MotorSampPts];
    end
    
    %% calculate the amount to equally sample from each session
    
    minsampspts = min([sum(mempts) sum(mixedpts) sum(motorpts)]);
    
    % gets the mean number of samples per session to sample (10
    % sessions total)
    
    meansamppts = floor(mintasksamps/10);
    meansampptstemp = meansamppts;

    fprintf('For subject %s the minimum number of sample points is %i, with a mean of %.0f points per sesssion per task for concatenation: %s\n', sub, minsampspts, meansamppts/3, datestr(now));
    
    % compensates for lack of motor data in MSC09 by sampling data from
    % the other 2 tasks
    if strcmp(sub, 'MSC09') && strcmp(task, 'motor') && ConcatenateTasks  %% Sets motor task data points for subject MSC09
        meansampptstemp = floor(sum(motorpts)/10);
    elseif strcmp(sub, 'MSC09') && ConcatenateTasks
        meansampptstemp = meansamppts + round((meansamppts - round(sum(motorpts)/10))/2);
    end
    
    %% Loop through sessions, account for missing sessions and sessions without enough data for sampling in each task
    
    remaindertotal = 0;
    notenoughdata = zeros(length(SubStruct),1);
    
    % if not enough data, will add sample points to match other
    % tasks/subjects
    while true
        remainder = 0;   %% update amount of points that still need to be sampled on each iteration
        sampspersession = zeros(length(SubStruct),1);
        notenoughdatatotal = sum(notenoughdata);
        meansampptstemp = (meansampptstemp + round(remaindertotal/(10-notenoughdatatotal)));
        
        % loop attempts to sample data equally per session. where this is
        % not possible due to a lack of high-quality data, it oversamples
        % equally across the remaining sessions. this process iteratively
        % repeats until either all of the data is sampled, or the target
        % amount of data (mintasksamps as previously calculated) is reached
        
        for v = 1:length(SubStruct)
            if strcmp(task, 'mem')
                if SubStruct(v).MemSampPts < meansampptstemp && notenoughdata(v) == 0
                    notenoughdata(v) = 1;
                    remainder = (remainder + (meansampptstemp - SubStruct(v).MemSampPts));
                    sampspersession(v) = SubStruct(v).MemSampPts;
                    fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                elseif notenoughdata(v) == 0
                    sampspersession(v) = meansampptstemp;
                    fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                else
                    sampspersession(v) = SubStruct(v).MemSampPts;
                    fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                end
            elseif strcmp(task, 'mixed')
                if SubStruct(v).MixedSampPts < meansampptstemp && notenoughdata(v) == 0
                    notenoughdata(v) = 1;
                    remainder = (remainder + (meansampptstemp - SubStruct(v).MixedSampPts));
                    sampspersession(v) = SubStruct(v).MixedSampPts;
                    fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                elseif notenoughdata(v) == 0
                    sampspersession(v) = meansampptstemp;
                    fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                else
                    sampspersession(v) = SubStruct(v).MixedSampPts;
                    fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                end
            elseif strcmp(task, 'motor')
                if SubStruct(v).MotorSampPts < meansampptstemp && notenoughdata(v) == 0
                    notenoughdata(v) = 1;
                    remainder = (remainder + (meansampptstemp - SubStruct(v).MotorSampPts));
                    sampspersession(v) = SubStruct(v).MotorSampPts;
                    fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                elseif notenoughdata(v) == 0
                    sampspersession(v) = meansampptstemp;
                    fprintf('For subject %s and session %s there are enough sample points, sampling %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                else
                    sampspersession(v) = SubStruct(v).MotorSampPts;
                    fprintf('For subject %s and session %s there are not enough sample points, sampling only %i sample points from this session: %s\n', sub, SubStruct(v).SessionFile, sampspersession(v), datestr(now));
                end
            end
        end
        
        % if needed, samples data from other sessions. stops loop when
        % enough data sampled
        if remainder == 0
            fprintf('Data sampling calculated for subject %s for task %s: %s\n', sub, task, datestr(now));
            break
        elseif (sum(notenoughdata) == 10) || (sum(notenoughdata) == 10 && remainder == 0)
            fprintf('Not enough data is present for this task to have enough data for subject %s for task %s: %s\n', sub, task, datestr(now));
            break
        else
            remaindertotal = remainder;
            fprintf('Remainder of %i calculated for all sessions for subject %s and task %s: %s\n', remaindertotal, sub, task, datestr(now));
        end
    end
    
else       %% if not matching data
    
    sampspersession = [];
    
end



