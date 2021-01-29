
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
% each session after applying the temporal mask. Once all the data has been
% concatenated, pairwise temporal correlations are calculated between all
% pairs of vertices. This generates a dconn which contains a seed map for
% each vertex (numvertices x numvertices). This dconn is then correlated
% row-wise with a group-average dconn to create a variant map (numvertices
% x 1), which shows the similarity of the seed map of a given vertex to the
% group average. The script can also output the dconns for each individual
% and the masked, concatenated timeseries.
%
% INPUTS:
% -outdir: an output directory for the created files (a separate directory
% is created within this directory for each of the outputs)
% -MSCTaskdir: filepath where MSC task data tmasks are located
% -NumSampsFiles_path: path to structures that contain the number of sample
% points per task and session (see CreateTMaskStruct_MSC.m)
% -cifti_filename: end of filename and extension for MSC timeseries files
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
% -ConcatenateTasks: toggles whether task data should be concatenated for
% the output files (set to 1), or if each task should concatenated
% separately across sessions (set to 1)
% -SubsampleAllTasks: this is only used for a special case. toggles whether
% the amount of matched data for all tasks combined should be equal to the
% matched amount of data for individual tasks (set to 1; see Kraus et al.,
% 2021, NIMG - Figure 7), otherwise set to 0
% -TaskResids: toggles whether file outputs should be created using task
% data that are the residuals of a GLM (set to 1), or task data that still
% include the task activations (set to 0)
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

outdir = '/your/path/here/'; %specify output directory
MSCTaskdir = '/your/path/here/'; %specify location of data
NumSampsFiles_path = '/your/path/here/';  %% path to the structure containing the number of good sample points per task
cifti_filename = 'yourfilenamehere.dtseries.nii';  %% file name and extension for MSC CIFTI timeseries
data_folder = '/your/path/here/';  %% folder name corresponding to data to use
template_path = '/your/path/here/template_cifti.dtseries.nii';  %% path to .dtseries.nii file to use as a template


SplitHalf = 1;  %% Toggles whether to create a separate file for odd/even sessions (set = 1)
MatchData = 0; %% Toggles whether to match the amount of data per task as the lowest amount across subjects (set = 1)
CortexOnly = 1; %% Toggles whether to include only cortical voxels in the correlation maps (set = 1)
ConcatenateTasks = 1;   %% Toggles whether to concatenate data for all tasks (set = 1)
SubsampleAllTasks = 0;  %% Toggle for special case; creates concatenated task data at the same length as the individual tasks (set = 1)
TaskResids = 0;  %% Toggles whether to use FC of task residuals (set = 1) or data with task activations
MakeDconn = 0;  %% Toggles whether to write a dconn (set = 1)
SaveTimeseries = 0;     %% Toggles whether to save concatenated timeseries for subject (set = 1)

subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed'};

if CortexOnly   % set the number of voxels for template
    voxnum = 59412;
else
    voxnum = 65625;
end

if ~ConcatenateTasks || SubsampleAllTasks  %% removes MSC09 from individual task comparisons
    subs(strcmp(subs,'MSC09')) = [];
end

%% Start analysis

fname_desc = variant_file_desc(MatchData,CortexOnly,TaskResids);  %% create filename description for naming outputs

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
if SubsampleAllTasks  %% special case for matching data across all tasks with data amount for individual tasks
    mintasksamps = round(mintasksamps/numel(tasks));
end
%% Main for-loop: makes Dconns
for i=1:numel(subs)
    
    fprintf('Running spatial correlation to group for subject %s, data description %s: %s\n', subs{i}, fname_desc, datestr(now));
  
    % Initialize concatenated session data (if concatenating tasks)
    if SplitHalf && ConcatenateTasks
    	catDataodd = [];
    	catDataeven = [];
    elseif ConcatenateTasks
    	catData = [];
    end
    
    for j=1:length(tasks)
        
        % Initialize concatenated session data (if individual tasks)
        if SplitHalf && ~ConcatenateTasks
            catDataodd = [];
            catDataeven = [];
        elseif ~ConcatenateTasks
            catData = [];
        end

        sampspersession = GetSampsPerSessionTask(mintasksamps, tasks{j}, subs{i}, NumSampsFiles_path, SplitHalf, MatchData, ConcatenateTasks);
        
        fprintf('Loading data for subject %s task %s: %s\n', subs{i}, tasks{j}, datestr(now));
        
        % Load tmasks from correct directories (FD filtered for
        % subs MSC03 and MSC10)
        if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')
            MSCcondidir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_pass2_FDfilt/']);
        else
            MSCcondidir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_pass2/']);
        end
        
        load ([MSCcondidir 'condindices.mat']);
        
        % Load vcids - DATALIST.txt is a file with paths to subjects'
        % sessions, other MSC-specific paths added below for data files and
        % temporal masks for high-quality data

        [~,vcids,~,~,~] = textread([MSCTaskdir '/' subs{i} '_' tasks{j} '_DATALIST.txt'],'%s%s%s%s%s');

        fprintf('%i sessions found for subject %s task %s: %s\n', size(vcids,1), subs{i}, tasks{j}, datestr(now));

        % Load data from correct directories (FD filtered for
        % subs MSC03 and MSC10) for either residuals or pre-GLM data
        if TaskResids
            if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')                    
                MSCdatadir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2_FDfilt/']);                    
            else                
                MSCdatadir = strcat(strcat(MSCTaskdir, '/'), ['FCProc_' subs{i} '_' tasks{j} '_pass2/']);                    
            end
        else
            if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')                    
                MSCdatadir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_preGLM_pass2_FDfilt/']);                    
            else                
                MSCdatadir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_preGLM_pass2/']);                    
            end
        end

        
        % Load and concatentate session (vcid) data
        for k=1:length(vcids)            
            fprintf('Loading timeseries for session %s for subject %s task %s: %s\n', vcids{k}, subs{i}, tasks{j}, datestr(now));            

            MSCciftidir = strcat(MSCdatadir, data_folder);  %% navigate to correct folder with processed surface data

            % try to load data file for each task/session, if no data file exists (from too little
            % high-quality data) then move on to the next session (vcid)
            try
                data = ft_read_cifti_mod([MSCciftidir vcids{k} cifti_filename]);  %% read data file with timeseries data
                data = data.data;
                if strcmp(tasks{j},'motor')   %% read temporal mask for each task from structure
                    tmask = TIndFin(k).AllMotor;
                elseif strcmp(tasks{j},'mem')
                    tmask = TIndFin(k).AllMem;
                else
                    tmask = TIndFin(k).AllGlass + TIndFin(k).AllSemantic;
                end

                if sum(tmask) == 0      % Skip file if no task task data                    
                    data = [];
                   fprintf('Vcid %s for subject %s has no usable data for task %s\n', vcids{k}, subs{i}, tasks{j});                    
                elseif ~MatchData               
                    data = data(1:voxnum,logical(tmask));                                       
                else                    
                    tmasknum = sampspersession(k);                
                    data = data(1:voxnum,logical(tmask));        %% Get good data from task
                    data = data(:,1:tmasknum);          %% Get correct amount of points corresponding to task                    
                end


                if SplitHalf               
                    fprintf('Data is %i by %i\n', size(data,1), size(data,2));                   
                    if mod(k,2) == 1                                
                        fprintf('Adding %i data points from session %s for subject %s to odd data\n', size(data,2), vcids{k}, subs{i});                        
                        catDataodd = [catDataodd data];                            
                    else                                
                        fprintf('Adding %i data points from session %s for subject %s to even data\n', size(data,2), vcids{k}, subs{i});                                
                        catDataeven = [catDataeven data];                            
                    end                         
                else
                    fprintf('Adding %i data points from session %s for subject %s to data\n', size(data,2), vcid, subs{i});
                    catData = [catData data];                        
                end                        
            catch ME
                if strcmp(ME.message,'Invalid file identifier.  Use fopen to generate a valid file identifier.')
                    continue
                end
            end
        end
        
        %% Saves task data either for individual tasks or concatenated tasks
        
        if ~ConcatenateTasks

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
                    
                    timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},tasks{j},'Odd','WashU120',fname_desc)];
                    
                    ft_write_cifti_mod(timeseries_fname,timeseriestemplate)
                    
                    fprintf('Data timseries is size %i by %i, %s\n', size(catDataeven,1), size(catDataeven,2), datestr(now));
                    
                    timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},tasks{j},'Even','WashU120',fname_desc)];
                    
                    ft_write_cifti_mod(timeseries_fname,timeseriestemplate2)
                    
                end
                
                % writes Dconn output files and creates output file names for variant maps
                
                fprintf('Running Correlations: on data size %i by %i, %s\n', size(catDataodd,1), size(catDataodd,2), datestr(now));
                
                variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},tasks{j},'Odd','WashU120',fname_desc)];
                
                dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},tasks{j},'Odd','WashU120',fname_desc)];
                
                make_dconn_wrapper_MSC(catDataodd,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)
                
                fprintf('Running Correlations: on data size %i by %i, %s\n', size(catDataeven,1), size(catDataeven,2), datestr(now));
                
                variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},tasks{j},'Even','WashU120',fname_desc)];
                
                dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},tasks{j},'Even','WashU120',fname_desc)];
                
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
                    timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},tasks{j},'All','WashU120',fname_desc)];
                    
                end
                
                % writes Dconn output files and creates output file names for variant maps
                
                fprintf('Running Correlations: on data size %i by %i, %s\n', size(catData,1), size(catData,2), datestr(now));
            
                variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},tasks{j},'All','WashU120',fname_desc)];
            
                dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},tasks{j},'All','WashU120',fname_desc)];
            
                make_dconn_wrapper_MSC(catData,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)
                
            end
        end
            
        
    end
    
    if ConcatenateTasks
        if SplitHalf
            if SaveTimeseries          %% Save concatenated FC processed timeseries
                
                % creates template for timeseries from MSC01 data
                timeseriestemplate = ft_read_cifti_mod(template_path);
                timeseriestemplate.data = [];
                timeseriestemplate2 = timeseriestemplate;
                
                % putting concatenated data in template
                timeseriestemplate.data = catDataodd;
                timeseriestemplate2.data = catDataeven;
                
                %write timeseries as cifti file
                
                fprintf('Data timseries is size %i by %i, %s\n', size(catDataodd,1), size(catDataodd,2), datestr(now));
                
                timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},'AllTasks','Odd','WashU120',fname_desc)];
                
                ft_write_cifti_mod(timeseries_fname,timeseriestemplate)
                
                fprintf('Data timseries is size %i by %i, %s\n', size(catDataeven,1), size(catDataeven,2), datestr(now));
                
                timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},'AllTasks','Even','WashU120',fname_desc)];
                
                ft_write_cifti_mod(timeseries_fname,timeseriestemplate2)
                
            end
            
            % writes Dconn output files and creates output file names for variant maps
            
            fprintf('Running Correlations: on data size %i by %i, %s\n', size(catDataodd,1), size(catDataodd,2), datestr(now));
            
            variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},'AllTasks','Odd','WashU120',fname_desc)];
            
            dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},'AllTasks','Odd','WashU120',fname_desc)];
            
            make_dconn_wrapper_MSC(catDataodd,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)
            
            fprintf('Running Correlations: on data size %i by %i, %s\n', size(catDataeven,1), size(catDataeven,2), datestr(now));
            
            variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},'AllTasks','Even','WashU120',fname_desc)];
            
            dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},'AllTasks','Even','WashU120',fname_desc)];
            
            make_dconn_wrapper_MSC(catDataeven,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)
            
        else
            if SaveTimeseries

                % creating template from MSC01 data
                timeseriestemplate = ft_read_cifti_mod(template_path);
                timeseriestemplate.data = [];
                
                %adding data to template
                timeseriestemplate.data = catData;
                
                fprintf('Data timseries is size %i by %i, %s\n', size(catData,1), size(catData,2), datestr(now));
                
                % write timeseries as cifti file
                timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},'AllTasks','All','WashU120',fname_desc)];
                
            end
            
            % writes Dconn output files and creates output file names for variant maps
            
            fprintf('Running Correlations: on data size %i by %i, %s\n', size(catData,1), size(catData,2), datestr(now));
            
            variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},'AllTasks','All','WashU120',fname_desc)];
            
            dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},'AllTasks','All','WashU120',fname_desc)];
            
            make_dconn_wrapper_MSC(catData,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)
            
        end
    end         
end




function fname_desc = variant_file_desc(MatchData,CortexOnly,TaskResids)

% create filename description

if MatchData && CortexOnly && TaskResids
    
    fname_desc = 'Resids-matched-cortexonly';

elseif MatchData && CortexOnly
    
    fname_desc = 'PreGLM-matched-cortexonly';
    
elseif MatchData && TaskResids
    
    fname_desc = 'Resids-matched-allvoxels';
    
elseif CortexOnly && TaskResids
    
    fname_desc = 'Resids-alldata-cortexonly';
    
elseif MatchData
    
    fname_desc = 'PreGLM-matched-allvoxels';
    
elseif CortexOnly
    
    fname_desc = 'PreGLM-alldata-cortexonly';
    
elseif TaskResids
    
    fname_desc = 'Resids-alldata-allvoxels';
    
else
    
    fname_desc = 'PreGLM-alldata-allvoxels';
    
end
end