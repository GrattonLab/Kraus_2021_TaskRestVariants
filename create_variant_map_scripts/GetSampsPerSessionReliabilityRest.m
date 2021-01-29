function sampspersession = GetSampsPerSessionReliabilityRest(sub,time,TR,NumSampsFiles_path)

%% Samples the desired amount of data from odd-numbered sessions equally for rest starting with the first session
%
% This function calculates the amount of rest data to sample from each
% session to sample the correct number of sample points of data to equal X 
% minutes. Rest data is sampled starting with the first session. When data 
% in this session runs out, data from the next odd-numbered session is 
% used. This continues until either the target amount of rest data is hit 
% or no more rest data exists.
%
% INPUTS:
% -NumSampsFiles_path: path to structures that contain the number of sample
% points per task and session (see CreateTMaskStruct_MSC.m)
% -sub: a subject ID to calculate data sampling for
% -time: a time (in minutes) to match the length of the conctenated rest 
% data to
% -TR: The TR (sample rate) of the data (in seconds), used for matching the
% number of sample points to the desired time
%
% OUTPUTS:
% -sampspersession: a vector of the amount of data points from rest to
% sample from each session
%
% Written by BK, edited by DP & ZL (01-2021)
%

%% Set up data sampling variables

load([NumSampsFiles_path sub '_NumSamps.mat']);

totalsamppts = round(time*(60/TR));  %% Calculate the number of sample points for a time (in minutes) for a given TR

ptscount = 0;  %% Running sample point count across sessions
enoughpts = 0;  %% Toggles whether enough sample points have been taken consecutively from each file
maxpts = 0;  %% Toggles whether the maximum amount of data for each task has been sampled

sampspersession = zeros(length(SubStruct),1);


%% Loops through odd rest sessions and samples data consecutively from each until the target amount is reached


% for rest, data is sampled consecutively starting with the first
% session. if all the rest data in a session is exhausted, data
% is then sampled from the next session until the target amount is reached.
% or no more data is present


for u = 1:length(SubStruct)
    
    %% loop through all odd sessions until enough data is sampled or no more data exists
    
    if SubStruct(u).OddEven == 1 && (~enoughpts)
        
        if SubStruct(u).RestSampPts + ptscount >= totalsamppts
            
            sampspersession(u) = sampspersession(u) + (totalsamppts - ptscount);
            ptscount = sampspersession(u) + ptscount;
            enoughpts = 1;
            
        elseif SubStruct(u).RestSampPts + ptscount < totalsamppts
            
            sampspersession(u) = sampspersession(u) + SubStruct(u).RestSampPts;
            ptscount = SubStruct(u).RestSampPts + ptscount;
            
        end
    end
    
    if u == length(SubStruct) && ~enoughpts
        
        maxpts = 1;
        
    end
end


% if not enough data, do not create a variant map. else, create one with
% the sampled data

if maxpts
    
    fprintf('Subject %s does not have %i minutes of total rest data. There are %i data points for rest (%.2f minutes total): %s\n', sub, time, ptscount, (ptscount*TR)/60, datestr(now));
    
    sampspersession = [];
    
else
    
    fprintf('Subject %s has enough data for %i minutes of total rest data. There are %i data points for rest (%.2f minutes total): %s\n', sub, time, ptscount, (ptscount*TR)/60, datestr(now));
    
end