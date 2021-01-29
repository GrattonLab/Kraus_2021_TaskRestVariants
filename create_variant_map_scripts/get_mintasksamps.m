function mintasksamps = get_mintasksamps(NumSampsFiles_path, subs, SplitHalf, MatchData)

%% Find the lowest amount of "high-quality" data that exists across all tasks and subjects
%
% If matching data, this function requires the output structures from
% CreateTMaskStruct_MSC.m, a subject ID, and binary inputs for whether the
% data should be treated as split-halves and whether data should be matched
% across subjects. If not matching data, set MatchedData to 0 and all other
% inputs are ignored.
%
% Given the input parameters, this function finds the lowest amount of task
% data that a participant has for any given analysis (e.g., the person with
% the least amount of task data in their odd numbered sessions). This value
% is then used as a target for data matching across sessions for each
% participant in the rest of the script
%
% INPUTS:
% -NumSampsFiles_path: path to structures that contain the number of sample
% points per task and session (see CreateTMaskStruct_MSC.m)
% -subs: a cell array of subject IDs to calculate data sampling for
% -SplitHalf: toggles whether data should be concatenated according to
% whether the sessions are odd/even numbered (set to 1), otherwise all
% session data is combined (set to 0)
% -MatchData: toggles whether the amount of data sampled should be matched
% as well as possible across subjects, tasks, and sessions (set to 1), or
% if all possible data should be sampled (set to 0)
%
% OUTPUTS:
% -mintasksamps: The lowest amount of available sample points across all 3
% tasks and subjects that can be sampled given the input parameters for
% sampling. Is an empty vector if not matching data
%
% Written by BK, edited by DP & ZL (01-2021)
%


%% Matches data points for all tasks across subjects
if MatchData
    
    % sets up variables for sum of sample points in each task of even- and 
    % odd-numbered sessions to be used later
    memptsoddsum = [];
    motorptsoddsum = [];
    mixedptsoddsum = [];
    memptsevensum = [];
    motorptsevensum = [];
    mixedptsevensum = [];
    
    %determines the number of sample points to match
    for i = 1:numel(subs)
        load([NumSampsFiles_path subs{i} '_NumSamps.mat']);
        
        % sets up variables for number of sample points available for
        % current subject
        memptsodd = [];
        motorptsodd = [];
        mixedptsodd = [];
        memptseven = [];
        motorptseven = [];
        mixedptseven = [];
        
        % QC files contain SubStruct with # of sample points for each task
        % and whether the session is an even- or odd-numbered session
        for u = 1:length(SubStruct)
            
            % divides sample points between even and odd sessions
            % 1 = Odd Session, 2 = Even Session
            if SubStruct(u).OddEven == 1
                memptsodd = [memptsodd; SubStruct(u).MemSampPts];
                mixedptsodd = [mixedptsodd; SubStruct(u).MixedSampPts];
                if ~strcmp(subs{i}, 'MSC09') % different for MSC09 because not enough motor sample points
                    motorptsodd = [motorptsodd; SubStruct(u).MotorSampPts];
                end
            elseif SubStruct(u).OddEven == 2
                memptseven = [memptseven; SubStruct(u).MemSampPts];
                mixedptseven = [mixedptseven; SubStruct(u).MixedSampPts];
                if ~strcmp(subs{i}, 'MSC09')
                    motorptseven = [motorptseven; SubStruct(u).MotorSampPts];
                end
                
            end
        end
        
        memptsoddsum = [memptsoddsum; sum(memptsodd)];
        mixedptsoddsum = [mixedptsoddsum; sum(mixedptsodd)];
        memptsevensum = [memptsevensum; sum(memptseven)];
        mixedptsevensum = [mixedptsevensum; sum(mixedptseven)];
        
        % sums all sample points for all even- or odd-numbered sessions for
        % current subject. sets motor points for MSC09 too high so that
        % they are not included in the determination for the minimum amount
        % of available data points for a given task (this is compensated
        % for by oversampling of the mixed and memory tasks).
        if strcmp(subs{i}, 'MSC09')
            motorptsoddsum = [motorptsoddsum; 9999];
            motorptsevensum = [motorptsevensum; 9999];
        else
            motorptsoddsum = [motorptsoddsum; sum(motorptsodd)];
            motorptsevensum = [motorptsevensum; sum(motorptseven)];
        end
    end
    
    %finds minimum number of sample points to match the other tasks across
    %subjects
    if SplitHalf
        mintasksamps = min(min([memptsoddsum motorptsoddsum mixedptsoddsum memptsevensum motorptsevensum mixedptsevensum]));
        fprintf('For all subjects, the minimum number of sample points available per task in a split-half is %i: %s\n', mintasksamps, datestr(now));
    else
        mintasksamps = min(min([(memptsoddsum + memptsevensum) (motorptsoddsum + motorptsevensum) (mixedptsoddsum + mixedptsevensum)]));
        fprintf('For all subjects, the minimum number of total sample points available per task is %i: %s\n', mintasksamps, datestr(now));
    end
    
else  %% if not matching data
    
    mintasksamps = [];

end
