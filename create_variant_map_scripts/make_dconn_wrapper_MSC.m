function make_dconn_wrapper_MSC(input_data,template_path,variant_fname,CortexOnly,save_dconn,dconn_fname)
 
%% This function is a wrapper for the function that creates variant maps, used to provide the necessary data and filepaths for computing the spatial correlation
% 
% This function requires a CIFTI timeseries (numvertices x numtimepoints)
% to compute the temporal correlations between each vertex and all other
% vertices a CIFTI file. If desired, this function will save out a file
% containing all these temporal correlations (a dconn). It will also
% provide this dconn to another function to create a spatial correlation
% (variant) map between the input dconn and a group average dconn.
%
% INPUTS:
% -input_data: a CIFTI timeseries (numvertices x numtimepoints)
% -template_path: path to template CIFTI file from which all output CIFTI
% files are written
% -variant_fname: filepath and filename for the output variant file
% (automatically generated in previous script)
% -dconn_fname: filepath and filename for the output dconn file
% (automatically generated in previous script)
% -CortexOnly: toggles whether the output CIFTI files and dconns should
% contain only cortical vertices on the surface (set to 1), or if they
% should also include subcortical voxels (set to 0)
% -save_dconn: toggles whether to save out the created conn file (set to 1)
% or whether to only use it for creating a variant map (set to 0)
%
% OUTPUTS:
% -dconn_data: matrix containing the values for the temporal correlation
% between each vertex and every other vertex. Used to make the variant map
% in the createSptlcorr_MSC.m function
%
% Written by BK (01-2021)
%

%% Initialize Variables

atlas_path = '/your/path/here/';  %% path to atlas folder
groupavg_name = 'groupsubs_corrmat';  %% name of group dconn to use as comparison

%% Create dconn, save if desired, create variant map

dconn_data = paircorr_mod(input_data');  %% run a temporal correlation between each vertex and every other vertex for the input timeseries

if save_dconn  %% if saving dconn
    
    if CortexOnly  %% set the number of voxels for template
        voxnum = 59412;
    else
        voxnum = 65625;   
    end
    
    % load template:
    template = ft_read_cifti_mod(template_path);
    
    % modify as needed
    template.data = [];
    template.dimord = 'pos_pos';
    template.hdr.dim(7) = voxnum;
    template.hdr.dim(6) = template.hdr.dim(7);
    template.hdr.intent_code = 3001;
    template.hdr.intent_name = 'ConnDense';
    
    % substitute data
    template.data = dconn_dat;
    
    % write dconn
    disp('writing dconn...');
    ft_write_cifti_mod(dconn_fname,template);
    disp('done');
    
end


createSptlcorr_MSC(atlas_path,groupavg_name,CortexOnly,template_path,dconn_data,variant_fname)