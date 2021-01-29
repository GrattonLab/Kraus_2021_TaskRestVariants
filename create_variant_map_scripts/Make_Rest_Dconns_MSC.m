
clear all

%% Creates a CIFTI variant map for each subject which contains the specified number of high-quality sample points for each task in each session
%
% This script requires CIFTI timeseries that have been converted to the
% surface (numvertices x timepoints), temporal masks for "high-quality"
% data points (timepoints x 1 binarized files), the output structures from
% CreateTMaskStruct_MSC.m (if matching data), and a group-average
% comparison dconn (numvertices x numvertices). For each subject, an
% individual dconn is created (sampling of data will vary according to the
% settings chosen, see below) from input CIFTIs by sampling the data from
% each session after applying the temporal mask. In this script, the task
% data is matched according to the settings first, then the rest data
% sampling is matched as close to possible to the task data sampling. Once
% all the data has been concatenated, pairwise temporal correlations are
% calculated between all pairs of vertices. This generates a dconn which
% contains a seed map for each vertex (numvertices x numvertices). This
% dconn is then correlated row-wise with a group-average dconn to create a
% variant map (numvertices x 1), which shows the similarity of the seed map
% of a given vertex to the group average. The script can also output the
% dconns for each individual and the masked, concatenated timeseries.
%
% INPUTS:
% -outdir: an output directory for the created files (a separate directory
% is created within this directory for each of the outputs)
% -rest_pathroot/rest_folder_stem: parts of the filepath where the MSC rest
% data is located
% -NumSampsFiles_path: path to structures that contain the number of sample
% points per task and session (see CreateTMaskStruct_MSC.m)
% -data_folder: folder containing FC processed surface data to load
% -template_path: path to template CIFTI file from which all output CIFTI
% files are written
% -subs: a cell array of subject IDs to calculate temporal correlations for
% -tasks: a cell array of tasks to calculate temporal correlations for
%
% -SplitHalf: toggles whether data should be concatenated according to
% whether the sessions are odd/even numbered (set to 1), otherwise all
% session data is combined (set to 0)
% -MatchData: toggles whether the amount of data sampled should be matched
% as well as possible across subjects, tasks, and sessions (set to 1), or
% if all possible data should be sampled (set to 0)
% -CortexOnly: toggles whether the output CIFTI files and dconns should
% contain only cortical vertices on the surface (set to 1), or if they
% should also include subcortical voxels (set to 0)
% -ConcatenateTasks: toggles whether task data (for matching to rest)
% should be concatenated for the output files (set to 1), or if each task
% should be concatenated and matched separately across sessions (set to 1)
% -SubsampleRest: this is only used for a special case. toggles whether
% the amount of matched rest data for all tasks combined should be equal to
% the matched amount of data for individual tasks (set to 1; see Kraus et
% al., 2021, NIMG - Figure 7), otherwise set to 0
% -MakeDconn: toggles whether the dconn for an individual subject should be
% saved out (set to 1), or only used for a spatial correlation with the
% group average (set to 0)
% -SaveTimeseries: toggles whether the masked, concatenated timeseries for
% an individual subject should be saved out (set to 1), or whether it
% should only be used to create the dconn for each subject (set to 0)
%
% OUTPUTS:
% -Variant Map: a spatial correlation map (numvertices x 1) that shows the
% similarity of the seed map for an individual at any given vertex to the
% group average comparison
% -Dconn (optional): a map containing the spatial correlations of every
% vertex to every other vertex for an individual (numvertices x
% numvertices; warning: very large!)
% -Timeseries (optional): a file containing the concatenated FC-processed
% timeseries containing only high-quality data points (numvertices x
% timepoints). this is used to match variants to network templates in later
% processing steps
%
% Written by BK, edited by DP & ZL (01-2021)
%

%% Initialize Variables

rest_pathroot = '/your/path/here/'; %specify location of rest data
outdir = '/your/path/here/'; %specify output directory
NumSampsFiles_path = '/your/path/here/';  %% path to the structure containing the number of good sample points per task
rest_folder_stem = '/your/path/here/';  %% location of rest data and temporal masks
data_folder = '/your/path/here/';  %% folder name corresponding to data to use
template_path = '/your/path/here/template_cifti.dtseries.nii';

SplitHalf = 1;  %% Toggles whether to calculate rest data sampling for split-half of sessions
MatchData = 0;  %% Toggles whether to match the amount of data rest data to the corresponding lowest amount of task data across subjects
CortexOnly = 1;  %% Toggles whether to include subcortical voxels in the correlation maps
ConcatenateTasks = 0;  %% Toggles whether to calculate rest data matched to the task data for concatenated tasks or individual tasks
SubsampleRest = 0;  %% Toggle for special case; creates concatenated rest data (matched equally to all of the task data) at the same length as the rest data matched to the individual tasks
MakeDconn = 0;  %% Toggles whether to write dconn to disk
SaveTimeseries = 0;  %% Toggles whether to save concatenated timeseries for subject

subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed'};


if CortexOnly   % set the number of voxels for template
    voxnum = 59412;
else
    voxnum = 65625;
end

if ~ConcatenateTasks || SubsampleRest  %% removes MSC09 from individual task comparisons
    subs(strcmp(subs,'MSC09')) = [];
end

%% Start Analysis 

if ~exist([outdir 'Variant_Maps'], 'dir')  %% create directories that don't exist
    mkdir([outdir 'Variant_Maps']);
end
if MakeDconn && ~exist([outdir 'Dconns'], 'dir')
    mkdir([outdir 'Dconns']);
end
if SaveTimeseries && ~exist([outdir 'Masked_Timeseries'], 'dir')
    mkdir([outdir 'Masked_Timeseries']);
end

%% Matches data points for all tasks across subjects
mintasksamps = get_mintasksamps(NumSampsFiles_path, subs, SplitHalf, MatchData);
if SubsampleRest
    mintasksamps = round(mintasksamps/numel(tasks));
end
%% Main for-loop: matches number of data points to sample from rest to task, makes Dconns
for i=1:numel(subs)
    
    fprintf('Running spatial correlation to group for subject %s for resting data: %s\n', subs{i}, datestr(now));
    
    % get the number of samples across all sessions and set whether rest
    % should be calculated for individual or concatenated tasks
    
    tasksampspersession = [];
    for a=1:numel(tasks)
        tasksampspersession = [tasksampspersession GetSampsPerSession(mintasksamps, tasks{a}, subs{i}, NumSampsFiles_path, SplitHalf, MatchData, ConcatenateTasks)];
    end
    
    if ConcatenateTasks
        taskloop = {'AllTasks'};
    elseif ~MatchData
        taskloop = {'AllTasks'};
    else
        taskloop = tasks;
    end
    %% Make Dconns
    
    % Initialize concatenated session data (if concatenating tasks)
    if SplitHalf && (ConcatenateTasks || ~MatchData)
        catDataodd = [];
        catDataeven = [];
    elseif ConcatenateTasks || ~MatchData
        catData = [];
    end
    
    for j=1:length(taskloop)  %% loop over individual tasks or just once for all tasks
        
        fname_desc = variant_file_desc(MatchData,CortexOnly,taskloop{j});  %% create filename description for naming outputs
        
        % Initialize concatenated session data (if individual tasks)
        if SplitHalf && ~ConcatenateTasks && MatchData
            catDataodd = [];
            catDataeven = [];
        elseif ~ConcatenateTasks && MatchData
            catData = [];
        end
        
        % if concatenating tasks, add the samples from each task across
        % sessions. otherwise, use the number of samples across sessions
        % for each task individually
        if ConcatenateTasks && MatchData
            matchtasksamps = sum(tasksampspersession,2);
        elseif MatchData
            matchtasksamps = tasksampspersession(:,j);
        else
            matchtasksamps = [];
        end
        
        
        % matches rest data sampling per session as closely as possible to
        % the corresponding task data
	    restsampspersession = GetSampsPerSessionRest(matchtasksamps, NumSampsFiles_path, SplitHalf, MatchData, subs{i}, taskloop{j});

        % Load rest vcids (sessions)
        restdir = dir([rest_pathroot subs{i} rest_folder_stem  data_folder]);
        vcidlist = restdir(contains({restdir.name},'dtseries') & contains({restdir.name},'vc'));
        
        fprintf('%i resting sessions found for subject %s: %s\n', size(vcidlist,1), subs{i}, datestr(now));
        
        % Load rest tmasks
        load([rest_pathroot subs{i} rest_folder_stem 'QC.mat'])
        
        % Load and concatentate data
        for k=1:length(vcidlist)
            vcid = vcidlist(k).name;
            fprintf('Loading rest timeseries for session %s and subject %s: %s\n', vcid, subs{i}, datestr(now));
            
            data = ft_read_cifti_mod([rest_pathroot subs{i} rest_folder_stem  data_folder vcid]);
            data = data.data;
            
            resttmask = QC(k).tmask;
            fprintf('tmask for rest file has %i high-quality data points, %s\n', sum(resttmask), datestr(now));
            
            data = data(1:voxnum,logical(resttmask));        %% select good timepoints from rest
            if MatchData == 1     %% Subsample data from the first voxel to the desired number of voxels (in this case the number of task data points per session)
                data = data(:,1:restsampspersession(k));          %% Get correct amount of data points per session corresponding to task data
            end
            
            if SplitHalf
                if mod(k,2) == 1
                    fprintf('Adding %i data points from session %s for subject %s to odd data\n', size(data,2), vcid, subs{i});
                    catDataodd = [catDataodd data];
                else
                    fprintf('Adding %i data points from session %s for subject %s to even data\n', size(data,2), vcid, subs{i});
                    catDataeven = [catDataeven data];
                end
            else
                catData = [catData data];
            end
        end
        
        %% Saves rest data matched versus specified tasks
        
        if SplitHalf
            
            if SaveTimeseries         %% Save concatenated FC processed timeseries
                
                % creates template for timeseries from MSC01 data
                timeseriestemplate = ft_read_cifti_mod(template_path);
                timeseriestemplate.data = [];
                timeseriestemplate2 = timeseriestemplate;
                
                % putting concatenated data in template
                timeseriestemplate.data = catDataodd;
                timeseriestemplate2.data = catDataeven;
                
                %write timeseries as cifti file
                
                fprintf('Data timseries is size %i by %i, %s\n', size(catDataodd,1), size(catDataodd,2), datestr(now));
                
                timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},'Rest','Odd','WashU120',fname_desc)];
                
                ft_write_cifti_mod(timeseries_fname,timeseriestemplate)
                
                fprintf('Data timseries is size %i by %i, %s\n', size(catDataeven,1), size(catDataeven,2), datestr(now));
                
                timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},'Rest','Even','WashU120',fname_desc)];
                
                ft_write_cifti_mod(timeseries_fname,timeseriestemplate2)
                
            end
            
            % writes Dconn output files and creates output file names for variant maps
            
            fprintf('Running Correlations: on data size %i by %i, %s\n', size(catDataodd,1), size(catDataodd,2), datestr(now));
            
            variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},'Rest','Odd','WashU120',fname_desc)];
            
            dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},'Rest','Odd','WashU120',fname_desc)];
            
            make_dconn_wrapper_MSC(catDataodd,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)
            
            fprintf('Running Correlations: on data size %i by %i, %s\n', size(catDataeven,1), size(catDataeven,2), datestr(now));
            
            variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},'Rest','Even','WashU120',fname_desc)];
            
            dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},'Rest','Even','WashU120',fname_desc)];
            
            make_dconn_wrapper_MSC(catDataeven,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)
            
        else
            
            if SaveTimeseries          %% Save concatenated FC processed timeseries
                % creating template from MSC01 data
                timeseriestemplate = ft_read_cifti_mod(template_path);
                timeseriestemplate.data = [];
                
                %adding data to template
                timeseriestemplate.data = catData;
                
                fprintf('Data timseries is size %i by %i, %s\n', size(catData,1), size(catData,2), datestr(now));
                
                % write timeseries as cifti file
                timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},'Rest','All','WashU120',fname_desc)];
                
            end
            
            % writes Dconn output files and creates output file names for variant maps
            
            fprintf('Running Correlations: on data size %i by %i, %s\n', size(catData,1), size(catData,2), datestr(now));
            
            variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},'Rest','All','WashU120',fname_desc)];
            
            dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},'Rest','Even','WashU120',fname_desc)];
            
            make_dconn_wrapper_MSC(catData,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)
            
        end
        
    end
end



function fname_desc = variant_file_desc(MatchData,CortexOnly,task)

% create filename description

if MatchData && CortexOnly
    
    fname_desc = ['Matchedvs' task '-cortexonly'];
    
elseif MatchData
    
    fname_desc = ['Matchedvs' task '-allvoxels'];
    
elseif CortexOnly
    
    fname_desc = 'Alldata-cortexonly';
    
else
    
    fname_desc = 'Alldata-allvoxels';
    
end
end
