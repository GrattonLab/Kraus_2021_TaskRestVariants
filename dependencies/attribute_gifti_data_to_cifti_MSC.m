function attribute_gifti_data_to_cifti_MSC(data_L, data_R, fname, atlasdir, workbenchdir, giftitemplatepath)

%% Write out GIFTI as CIFTI
%
% This function takes two GIFTI files as inputs (one for each hemisphere)
% and reconsitutes them into a single CIFTI file containing both
% hemispheres (after accounting for the medial wall). This CIFTI is then 
% saved out as a .dtseries.nii file
%
% INPUTS:
%
% -data_L: the filepath and filename for a GIFTI left hemisphere file
% -data_R: the filepath and filename for a GIFTI right hemisphere file
% -fname: the filepath and filename of a CIFTI file to create from the
% GIFTI files
% -atlasdir: the path to the directory that contains the 32k surface medial
% wall GIFTI files
% -giftitemplatepath: the filepath and filename to the GIFTI used as a
% template for converting GIFTIs to CIFTIs
% -workbenchdir: the filepath to connectome workbench
%
% -OUTPUTS:
%
% -passes concantenated GIFTI data to cifti_write_wHDR_MSC.m to write out
% the CIFTI file
%
%
% Originally from WashU, edited by BK (01-2021)
%

order_R = gifti(data_R);  %% load GIFTI data for right hemisphere
order_L = gifti(data_L);  %% load GIFTI data for leftt hemisphere
order_L = order_L.cdata;
order_R = order_R.cdata;

medial_R = gifti([atlasdir 'medial_wall.R.32k_fs_LR.func.gii']);  %% load mask for medial wall from 32k atlas for right hemisphere
medial_L = gifti([atlasdir 'medial_wall.L.32k_fs_LR.func.gii']);  %% load mask for medial wall from 32k atlas for left hemisphere
medial_R = medial_R.cdata;
medial_L = medial_L.cdata;

cifti_order = [order_L(~medial_L,:); order_R(~medial_R,:)];  %% exclude medial wall from GIFTI
cifti_order(size(cifti_order,1):66697,:) = 0;  %% set order of vertices in CIFTI accounting for medial wall in the beginning

cifti_write_wHDR_MSC(cifti_order, [], fname, workbenchdir, giftitemplatepath, 'dtseries');  %% write CIFTI from GIFTI data

