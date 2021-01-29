function restsampspersession = GetSampsPerSessionRest(tasksampspersession,NumSampsFiles_path,SplitHalf,MatchData,sub,task)

%% Calculates the amount of data to sample from each session for rest in order to match data across subjects, tasks, and sessions
%
% This function calculates the amount of data to sample from each session 
% for rest for each subject given the input parameters. Using the amount of
% data calculated previously (tasksampspersession), the function determines
% how much data must be sampled from each session of rest based on how
% much data was sampled from each session for the task data. The function 
% attempts to first sample the rest data exactly the same as the 
% corresponding task data was sampled. When this is not possible due a lack
% of "high-quality" data in a session, the remaining sessions are used to 
% equally sample the required amount of rest data points from. This process 
% repeats until either all rest data for a subject has been sampled, or the
% desired total amount of data is reached.
%
% INPUTS:
% -tasksampspersession: the total target amount of data to sample across 
% sessions for rest (calculated previously in script if matching data)
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
%
% OUTPUTS:
% -restsampspersession: a vector of the amount of data points from rest to 
% sample from each session. Is an empty vector if not matching data
%
% Written by BK, edited by DP & ZL (01-2021)
%


if MatchData
    %% load and extract data from QC structure
    
    fprintf('Matching sampling of resting state data for subject %s to data for task %s: %s\n', sub, task, datestr(now));
    
    load ([NumSampsFiles_path sub '_NumSamps.mat']);
    
    if SplitHalf
        
        % Calculate rest data for each session based on task data
        remaindertotalodd = 0;
        remaindertotaleven = 0;
        notenoughdatarest = zeros(length(SubStruct),1);
        
        % checks if rest runs have enough data to match task data. If not
        % enough, samples from other rest sessions until total number of
        % data points matches
        
        while true
            
            remainderodd = 0;
            remaindereven = 0;
            restsampspersession = zeros(length(SubStruct),1);
            
            notenoughdata = sum(notenoughdatarest);
            notenoughdataeven = sum(notenoughdatarest(2:2:end));
            notenoughdataodd = sum(notenoughdatarest(1:2:end));
            
            if notenoughdataeven > 0 || notenoughdataodd > 0
                for b = 1:length(SubStruct)
                    if SubStruct(b).OddEven == 1 && notenoughdatarest(b) == 0
                        tasksampspersession(b) = (tasksampspersession(b) + round(remaindertotalodd/(5-notenoughdataodd)));
                    elseif SubStruct(b).OddEven == 2 && notenoughdatarest(b) == 0
                        tasksampspersession(b) = (tasksampspersession(b) + round(remaindertotaleven/(5-notenoughdataeven)));
                    end
                end
            end
            
            for c = 1:length(SubStruct)
                if SubStruct(c).OddEven == 1
                    if SubStruct(c).RestSampPts <  tasksampspersession(c) && notenoughdatarest(c) == 0
                        notenoughdatarest(c) = 1;
                        remainderodd = (remainderodd + (tasksampspersession(c) - SubStruct(c).RestSampPts));
                        restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                        fprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                    elseif notenoughdatarest(c) == 0
                        restsampspersession(c) = tasksampspersession(c);
                        
                        fprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                    else
                        restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                        fprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                    end
                else
                    if SubStruct(c).RestSampPts < tasksampspersession(c) && notenoughdatarest(c) == 0
                        notenoughdatarest(c) = 1;
                        remaindereven = (remaindereven + (tasksampspersession(c) - SubStruct(c).RestSampPts));
                        restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                        fprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                    elseif notenoughdatarest(c) == 0
                        restsampspersession(c) = tasksampspersession(c);
                        
                        fprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                    else
                        restsampspersession(c) = SubStruct(c).RestSampPts;
                        
                        fprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                    end
                end
            end
            
            if remainderodd == 0 && remaindereven == 0
                fprintf('Sampling calculated for subject %s rest data matched with task %s: %s\n', sub, task, datestr(now));
                break
            elseif (sum(notenoughdata) == 10) || (sum(notenoughdataeven) == 5 && remainderodd == 0) || (sum(notenoughdataodd) == 5 && remaindereven == 0)
                fprintf('Not enough resting data is present to have enough data for matching subject %s and task %s: %s\n', sub, task, datestr(now));
                break
            else
                remaindertotalodd = remainderodd;
                remaindertotaleven = remaindereven;
                fprintf('Remainder of %i calculated for odd sessions for subject %s rest data matched with task %s: %s\n', remainderodd, sub, task, datestr(now));
                fprintf('Remainder of %i calculated for even sessions for subject %s rest data matched with task %s: %s\n', remaindereven, sub, task, datestr(now));
            end
        end
        
    else
        
        remaindertotal = 0;
        notenoughdatarest = zeros(length(SubStruct),1);
        
        % checks if rest runs have enough data to match task data. If not
        % enough, samples from other rest sessions until number of data
        % points matches
        while true
            remainder = 0;
            restsampspersession = zeros(length(SubStruct),1);
            notenoughdata = sum(notenoughdatarest);
            if notenoughdata > 0
                for b = 1:length(SubStruct)
                    tasksampspersession(b) = (tasksampspersession(b) + round(remaindertotal/(10-notenoughdata)));
                end
            end
            
            for c = 1:length(SubStruct)
                if SubStruct(c).RestSampPts <  tasksampspersession(c) && notenoughdatarest(c) == 0
                    notenoughdatarest(c) = 1;
                    remainder = (remainder + (tasksampspersession(c) - SubStruct(c).RestSampPts));
                    restsampspersession(c) = SubStruct(c).RestSampPts;
                    
                    fprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                elseif notenoughdatarest(c) == 0
                    restsampspersession(c) = tasksampspersession(c);
                    
                    fprintf('For subject %s and session %s there are enough sample points in the resting data to match the task, sampling %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                else
                    restsampspersession(c) = SubStruct(c).RestSampPts;
                    
                    fprintf('For subject %s and session %s there are not enough sample points in the resting data to match the task, sampling only %i sample points from this session: %s\n', sub, SubStruct(c).SessionFile, restsampspersession(c), datestr(now));
                end
            end
            
            if remainder == 0
                fprintf('Sampling calculated for subject %s rest data matched with task %s: %s\n', sub, task, datestr(now));
                break
            elseif (sum(notenoughdata) == 10)
                fprintf('Not enough resting data is present to have enough data for matching subject %s and task %s: %s\n', sub, task, datestr(now));
                break
            else
                remaindertotal = remainder;
                fprintf('Remainder of %i calculated for all sessions for subject %s rest data matched with task %s: %s\n', remainder, sub, task, datestr(now));
            end
        end
    end
    
else  %% else, if not matching data
    
    restsampspersession= [];  %% set output if not matching data
    
end


