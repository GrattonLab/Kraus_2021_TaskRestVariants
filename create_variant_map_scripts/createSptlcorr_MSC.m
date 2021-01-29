function createSptlcorr_MSC(atlas_path,groupavg_name,CortexOnly,template_path,dconn_data,variant_fname)

%% Creates a variant map from the spatial correlation between an individual's temporal correlations at each vertex and group average temporal correlations at the same vertex
%
% This function correlates the temporal correlations between each voxel and
% all other voxels (a seed map) between an individual's seed map and a
% group average seed map. An individual's seed map file (dconn) must be
% input along with the path and name of the group average dconn to compare
% it to. The function saves out a variant map (numvertices x 1) CIFTI,
% where each vertex containts a correlation value representing the
% similarity of the seed map of that vertex between an individual and the
% group average. Lower values indicate that there is less similarity at a
% given vertex between that individual's seed map and the group seed map.
%
% INPUTS:
%
% -atlas_path: the path to the directory that contains the group average
% dconn to use for creating the variant map
% -groupavg_name: the filename for the group average dconn to use for 
% creating the variant map
% -CortexOnly: toggles whether the output CIFTI files and dconns should
% contain only cortical vertices on the surface (set to 1), or if they
% should also include subcortical voxels (set to 0)
% -template_path: path to template CIFTI file from which all output CIFTI
% files are written
% -dconn_data: matrix containing the values for the temporal correlation
% between each vertex and every other vertex
% -variant_fname: filepath and filename for the output variant file
% (automatically generated in previous script)
%
% OUTPUTS:
%
% -Variant Map: A CIFTI file (numvertices x 1) where each vertex contains a
% correlation value representing the similarity of the seed map of that 
% vertex between an individual and the group average. Lower values indicate
% that there is less similarity at a given vertex between that individual's
% seed map and the group seed map
%
% Original code from BAS at WashU, edited by BK (01-2021)
%

    
% Remove NaNs (produced if vertices have no data)
dconn_data(isnan(dconn_data)) = 0;

% Apply the Fisher tranformation to all values in the input dconn
% row-by-row
for dconnrow = 1:size(dconn_data,1)
    dconn_data(dconnrow,:) = FisherTransform(dconn_data(dconnrow,:));
end

% Load group-average dconn
group = ft_read_cifti_mod([atlas_path '/' groupavg_name '.dconn.nii']);
if CortexOnly %% select the correct number of vertices
    group = group.data(1:59412,1:59412);
    
else
    group = group.data;
    
end

% Load template variants file
template = ft_read_cifti_mod(template_path);
template.data = [];

% Compare to group average
for i=1:size(group,1)
    template.data(i,1) = paircorr_mod(group(:,i),dconn_data(:,i));
end

% Write out the variants
ft_write_cifti_mod(variant_fname,template);


 

