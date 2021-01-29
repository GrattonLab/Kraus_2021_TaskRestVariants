function excluded_dat = variant_size_exclude(cifti_dat,minsize)

%% Excludes contiguous variants smaller than an input size threshold (in vertices)
%
% The input for this function is a cifti object containing a UniqueID
% variant map. Any variant smaller than the input minimum size (in
% vertices) is removed from the input cifti.data vector. Once all
% contiguous variants have been checked, the resulting vector with the
% excluded variants set to 0 is output.
%
% INPUTS:
%
% -cifti_dat: a cifti object (from a UniqueID map, see 
% threshold_variant_maps.m for additional documentation) to exclude 
% variants smaller than a certain size from
% -minsize: the minimum size threshold (in vertices) that a variant has to
% meet in order to be included in the output data
%
% OUTPUTS:
%
% -excluded_dat: a vector containing the variants from the UniqueID map
% that survived the size exclusion
%
% Written by BK (01-2021)
%

excluded_dat = zeros(size(cifti_dat.data));  %% empty vector for excluded data

allvars = unique(cifti_dat.data(cifti_dat.data>0));  %% get the unique ID number for each variant in each file

for var = 1:length(allvars)  %% loop through variants
    
    if length(find(cifti_dat.data == allvars(var))) >= minsize  %% if variant has at least as many vertices as the size threshold
        
        excluded_dat(cifti_dat.data == allvars(var)) = allvars(var);  %% add that variant's UniqueID value to its indices in the size excluded data
        
    end
end