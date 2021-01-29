
clear all

%% Create variant maps in X minute increments with odd numbered sessions for test split-half using combined task data
%
% This script requires CIFTI timeseries that have been converted to the
% surface (numvertices x timepoints), temporal masks for "high-quality"
% data points (timepoints x 1 binarized files), the output structures from
% CreateTMaskStruct_MSC.m, and a group-average comparison dconn
% (numvertices x numvertices). For each subject, an individual dconn is
% created from input CIFTIs by sampling the data from each session after 
% applying the temporal mask up until the target time in minutes is 
% reached. Data is sampled equally for each task starting with the first 
% odd-numbered session until either the target amount of data is reached or
% no more task data is available. Once no more task data is available, data
% from other tasks is sampled to compensate. Once the target amount of data
% is concatenated, pairwise temporal correlations are calculated between 
% all pairs of vertices. This generates a dconn which contains a seed map 
% for each vertex (numvertices x numvertices). This dconn is then 
% correlated row-wise with a group-average dconn to create a variant map 
% (numvertices x 1), which shows the similarity of the seed map of a given 
% vertex to the group average. The script can also output the dconns for 
% each individual and the masked, concatenated timeseries.
%
% INPUTS:
% -outdir: an output directory for the created files (a separate directory
% is created within this directory for each of the outputs)
% -MSCTaskdir: filepath where MSC task data tmasks are located
% -NumSampsFiles_path: path to structures that contain the number of sample
% points per task and session (see CreateTMaskStruct_MSC.m)
% -data_folder: folder containing FC processed surface data to load
% -template_path: path to template CIFTI file from which all output CIFTI
% files are written
% -subs: a cell array of subject IDs to calculate temporal correlations for
% -tasks: a cell array of tasks to calculate temporal correlations for
% -times: a vector of times (in minutes) to match the length of the
% concatenated data to
% -TR: The TR (sample rate) of the data (in seconds), used for matching the
% number of sample points to the desired time
%
% -CortexOnly: toggles whether the output CIFTI files and dconns should
% contain only cortical vertices on the surface (set to 1), or if they
% should also include subcortical voxels (set to 0)
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
data_folder = '/your/path/here/';  %% folder name corresponding to data to use
template_path = '/your/path/here/template_cifti.dtseries.nii';  %% path to .dtseries.nii file to use as a template

subs = {'MSC01','MSC02','MSC03','MSC04','MSC05','MSC06','MSC07','MSC09','MSC10'};
tasks = {'motor','mem','mixed'};

CortexOnly = 1;  %% Toggles whether to run correlations on cortex only
TaskResids = 0;  %% Toggles whether to use task residuals (= 1) or data with task activations
MakeDconn = 0;  %% Toggles whether to write a dconn
SaveTimeseries = 0;  %% Save concatenated timeseries for subject

times = [5:5:100];  %% Times (in minutes) to make spatial correlation maps
TR = 2.2;  %% The TR of the data (in seconds)


if CortexOnly      %% Select correct number of voxels for template
    voxnum = 59412;
else
    voxnum = 65625;
end

if ~exist([outdir 'Variant_Maps'], 'dir')  %% create directories that don't exist
    mkdir([outdir 'Variant_Maps']);
end
if MakeDconn && ~exist([outdir 'Dconns'], 'dir')
    mkdir([outdir 'Dconns']);
end
if SaveTimeseries && ~exist([outdir 'Masked_Timeseries'], 'dir')
    mkdir([outdir 'Masked_Timeseries']);
end

%% Loop over subjects to create variant maps for specific times for odd sessions

for i=1:numel(subs)
    
    for l = 1:numel(times)
    
        fprintf('Creating variant map for subject %s using %i minutes of data: %s\n', subs{i}, times(l), datestr(now));

        fname_desc = variant_file_desc(CortexOnly,times(l),TaskResids);  %% create filename description for naming outputs

        % Initialize concatenated session data
      	catData = [];
        
        % Calculate data sampling for odd numbered (test) sessions
        
        [sampspersessionmem, sampspersessionmixed, sampspersessionmotor, memmissing, mixedmissing, motormissing] = GetSampsPerSessionReliabilityTask(subs{i},times(l),TR,NumSampsFiles_path);
        
        if isempty([sampspersessionmem; sampspersessionmixed; sampspersessionmotor])   %% only create a variant map if enough data exists
            
            fprintf('Skipping variant map for subject %s, not enough data (%i minutes): %s\n', subs{i}, times(l), datestr(now));
            
        else
        
            for j=1:length(tasks)

                fprintf('Running spatial correlation to group for subject %s task %s: %s\n', subs{i}, tasks{j}, datestr(now));

                % Load tmasks
                if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')
                    MSCcondidir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_pass2_FDfilt/']);
                else
                    MSCcondidir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_pass2/']);
                end

                load ([MSCcondidir 'condindices.mat']);

                % remove tmasks for non-existent runs (due to low amount of
                % high-quality data)

                if ~isempty(memmissing) && strcmp(tasks{j},'mem')
                    TIndFin(memmissing) = [];   
                end

                if ~isempty(mixedmissing) && strcmp(tasks{j},'mixed')
                    TIndFin(mixedmissing) = [];
                end

                if ~isempty(motormissing)  && strcmp(tasks{j},'motor')
                    TIndFin(motormissing) = [];
                end 

                % set correct directory for data

                if TaskResids
                    if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')
                        MSCciftidir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_pass2_FDfilt/' data_folder]);
                    else
                        MSCciftidir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_pass2/' data_folder]);
                    end
                else
                    if strcmp(subs{i}, 'MSC03') || strcmp(subs{i}, 'MSC10')
                        MSCciftidir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_preGLM_pass2_FDfilt/' data_folder]);
                    else
                        MSCciftidir = strcat(MSCTaskdir, ['FCProc_' subs{i} '_' tasks{j} '_preGLM_pass2/' data_folder]);
                    end
                end

                % Load and concatentate session (vcid) data

                ciftifiles = dir(MSCciftidir);
                vcidlist = ciftifiles(contains({ciftifiles.name},'vc'));

                fprintf('%i sessions found for subject %s task %s: %s\n', size(vcidlist,1), subs{i}, tasks{j}, datestr(now));

                for k=1:length(vcidlist)

                    vcid = vcidlist(k).name;

                    fprintf('Loading timeseries for session %s for subject %s task %s: %s\n', vcid, subs{i}, tasks{j}, datestr(now));

                        data = ft_read_cifti_mod([MSCciftidir '/' vcid]);
                        data = data.data(1:voxnum,:);
                        if strcmp(tasks{j},'motor')
                            tmask = TIndFin(k).AllMotor;
                            sampspersession = sampspersessionmotor(k);
                        elseif strcmp(tasks{j},'mem')
                            tmask = TIndFin(k).AllMem;
                            sampspersession = sampspersessionmem(k);
                        else
                            tmask = TIndFin(k).AllGlass + TIndFin(k).AllSemantic;
                            sampspersession = sampspersessionmixed(k);
                        end

                        if sum(tmask) == 0      % Skip file if no task task data

                            data = [];
                            fprintf('Vcid %s for subject %s has no usable data for task %s\n', vcid, subs{i}, tasks{j});

                        else  %% select correct number of samples from each task

                            tmasknum = sampspersession;

                            if tmasknum > 0  
                                data = data(1:voxnum,logical(tmask));
                                if tmasknum < sum(tmask)
                                    data = data(:,1:tmasknum);    
                                end
                            else
                                data = [];
                            end
                        end

                        fprintf('Adding %i data points from session %s for subject %s to "test" data\n', size(data,2), vcid, subs{i});
                        catData = [catData data];
                end
            end
            
            %% Save out data

            if SaveTimeseries
                
                % creating template from MSC01 data
                timeseriestemplate = ft_read_cifti_mod(template_path);
                timeseriestemplate.data = [];
                
                %adding data to template
                timeseriestemplate.data = catData;
                
                fprintf('Data timseries is size %i by %i, %s\n', size(catData,1), size(catData,2), datestr(now));
                
                % write timeseries as cifti file
                timeseries_fname = [outdir 'Masked_Timeseries/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_MaskedTimeseries.dtseries.nii',subs{i},'AllTasks','Odd','WashU120',fname_desc)];
                
            end
            
            % writes Dconn output files and creates output file names for variant maps
            
            fprintf('Running Correlations: on data size %i by %i, %s\n', size(catData,1), size(catData,2), datestr(now));
            
            variant_fname = [outdir 'Variant_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_VariantMap.dtseries.nii',subs{i},'AllTasks','Odd','WashU120',fname_desc)];
            
            dconn_fname = [outdir 'Dconns/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s',subs{i},'AllTasks','Odd','WashU120',fname_desc)];
            
            make_dconn_wrapper_MSC(catData,template_path,variant_fname,CortexOnly,MakeDconn,dconn_fname)

        end
    end
end





function fname_desc = variant_file_desc(CortexOnly,time,TaskResids)

% create filename description

if CortexOnly && TaskResids
    
    fname_desc = sprintf('Resids-testhalf-cortexonly-%0.fminutes',time);

elseif CortexOnly
    
    fname_desc = sprintf('PreGLM-testhalf-cortexonly-%0.fminutes',time);
    
elseif TaskResids
    
    fname_desc = sprintf('Resids-testhalf-allvoxels-%0.fminutes',time);
    
else
    
    fname_desc = sprintf('PreGLM-testhalf-allvoxels-%0.fminutes',time);
    
end
end

