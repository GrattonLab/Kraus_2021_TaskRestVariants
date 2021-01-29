function threshold_variant_maps(thresholds,SNRexclusion,absolutethresholds,outputdir,spatialmapslist,workbenchdir,surfdir,SNRpath)

%% Create binarized and unique ID files (and text files with filepaths and subject IDs) for each input variant map
%
% This function reads a text file in the format 'pathtofile subjectid' with
% spaces as delimiters. Then, loads each variant map in the text file and
% thresholds it at the desired thresholds. Both a binarized map (1s and 0s)
% and a uniqueID map (each contiguous variant has a unique identifier) are
% created at each threshold input. Each of these filetypes is saved in a
% separate directory created within the output directory provided, along
% with an associated text file in the same format as the input text file
%
% INPUTS:
%
% -workbenchdir: the filepath to connectome workbench
% -surfdir: the filepath to the directory that has a path to a template for
% the surface of the brain
% -SNRpath: the filepath to the directory that contains the SNR CIFTIs for
% each subject. These files were created by mode 1000 normalizing the
% resting state data for each session within each subject. The mean BOLD
% value for each voxel in each session was calculated, and then this SNR
% value for each voxel was averaged across sessions. This data was then
% transformed to a CIFTI (numvertices x 1) and used for excluding low SNR
% regions (mean signal < 750 after mode 1000 normalization)
% -spatialmapslist: the filepath to the text file that contains a
% filepaths for each variant map and it's associated subject ID (space
% delimited) on each line of the file
% (e.g.,filepath/MSC01.dtseries.nii MSC01)
% -outputdir: the output directory for the thresholded variant files (note
% that subdirectories are created within this directory for the output
% files)
% -thresholds: a vector (or single value) of thresholds to define variants
% for
% -SNRexclusion: toggles whether to exclude low signal regions from
% consideration as variants (set to 1), otherwise allow all vertices to be
% defined as variants (set to 0)
% -absolutethresholds: toggles whether to use an absolute correlation (r)
% value to define variants (set to 1), otherwise thresholds are interpreted
% as the xth lowest percentage of correlation values (set to 0)
%
% -Lhemipath: the filename and extension for the left hemisphere surface
% midthickness GIFTI
% -Rhemipath: the filename and extension for the right hemisphere surface
% midthickness GIFTI
% -SNR_fstring: the filename and extension for the SNR map
%
% OUTPUTS: 
%
% -Binarized Map: A CIFTI file (numvertices x 1) with binarized (0s and 1s)
% values. The 1s correspond to vertices that were defined as variants given
% the inputs provided, and the 0s are vertices which were not defined as
% variants. This map is then used to create a UniqueID map using connectome
% workbench
% -UniqueID Map: A CIFTI file (numvertices x 1) where every contiguous
% parcel (variant) on the surface has been given a unique value (i.e., all
% the number 1s are part of a contiguous parcel that is distinct from all
% the number 2s). This allows easy identification of distinct parcels
% (variants) for further analysis
%
% Written by BK (01-2021)
%


%% Initialize Variables

Lhemipath = [surfdir 'Conte69.L.midthickness.32k_fs_LR.surf.gii'];  %% filename and extension for left hemisphere surface midthickness 
Rhemipath = [surfdir 'Conte69.R.midthickness.32k_fs_LR.surf.gii'];  %% filename and extension for right hemisphere surface midthickness 
SNR_fstring = 'SNRname.dtseries.nii';  %% filename and extension for SNR map

[variant_files, subjects] = textread(spatialmapslist,'%s%s'); % read in list of variant maps

% check and create output directories

if ~exist([outputdir 'Binarized_Maps'], 'dir')
    mkdir([outputdir 'Binarized_Maps']);
end
if ~exist([outputdir 'UniqueID_Maps'], 'dir')
    mkdir([outputdir 'UniqueID_Maps']);
end

%% Loop through each threshold and create a set of binarized and unique ID spatial maps for each

for thresh = 1:numel(thresholds)
    
    binarizedfnames = {};  % store filenames at each threshold
    uniqueIDfnames = {};
    
    threshold = thresholds(thresh);
    
    for file = 1:length(variant_files)  % loop through files
        
        subject = subjects{file};
        
        [sub, task, ses, grp, desc] = read_variant_filename(variant_files{file});  % read different parts of variant filename
        
        fprintf('Writing variants for subject %s, task %s, and session %s\n',sub, task, ses);
        
        cifti_file = ft_read_cifti_mod(variant_files{file});  % read each variant map
        
        if SNRexclusion   %% exclude individual SNR map
            
            SNRmap = ft_read_cifti_mod([SNRpath subject SNR_fstring]);
            SNRmap.data = SNRmap.data(1:59412,:);
            SNRexclude = find(SNRmap.data < 750);
            
            cifti_file.data(SNRexclude) = NaN;
            
        end
        
        if absolutethresholds   %% threshold variant map
            
            cifti_file_threshold = find(cifti_file.data < threshold);
            
            thresh_str = [num2str(threshold) 'r'];
            
        else
            
            cifti_file_threshold = find(cifti_file.data < prctile(cifti_file.data,threshold));
            
            thresh_str = [num2str(threshold) 'pct'];
            
        end
        
        % apply threshold to map
        
        cifti_file_thresh_dat = zeros(size(cifti_file.data));
        cifti_file_thresh_dat(cifti_file_threshold) = 1;
        
        cifti_file.data = cifti_file_thresh_dat;
        
        % create output filenames for each map
        
        outfilebinarized = [outputdir 'Binarized_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_BinarizedMap.dtseries.nii',sub,task,ses,grp,desc,thresh_str)];
        outfileuniqueID = [outputdir 'UniqueID_Maps/' sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_UniqueIDMap.dtseries.nii',sub,task,ses,grp,desc,thresh_str)];
        
        % write out binarized file
        
        ft_write_cifti_mod(outfilebinarized, cifti_file);
        
        % use binarized file to create the unique ID file with connectome workbench
        
        system([workbenchdir 'wb_command -cifti-find-clusters ' outfilebinarized ' 0 0 0 0 COLUMN ' outfileuniqueID ' -left-surface ' Lhemipath ' -right-surface ' Rhemipath]);
        
        % save each filename and subject ID for output text file
        
        binarizedfnames = [binarizedfnames;{outfilebinarized, subject}];
        uniqueIDfnames = [uniqueIDfnames;{outfileuniqueID subject}];
        
    end
    
    % write output file for each threshold
    
    write_textfiles_variants([outputdir 'Binarized_Maps/Binarized_Maps_thresh-' thresh_str '_SubjectList.txt'], binarizedfnames(:,1),binarizedfnames(:,2));
    write_textfiles_variants([outputdir 'UniqueID_Maps/UniqueID_Maps_thresh-' thresh_str '_SubjectList.txt'], uniqueIDfnames(:,1),uniqueIDfnames(:,2));
    
end





function [sub, task, ses, grp, desc] = read_variant_filename(var_path)

% Read parts of each filename from the variant filepath

sub = cell2mat(extractBetween(var_path,'sub-','_'));
task = cell2mat(extractBetween(var_path,'task-','_'));
ses = cell2mat(extractBetween(var_path,'ses-','_'));
grp = cell2mat(extractBetween(var_path,'grp-','_'));
desc = cell2mat(extractBetween(var_path,'desc-','_'));



function write_textfiles_variants(outputpath,filepaths,subjects)

% write text files with filepath and subject ID for each output map

fid = fopen(outputpath,'wt');
for file = 1:length(filepaths)
fprintf(fid,'%s ',filepaths{file});
fprintf(fid,'%s\n',subjects{file});
end
fclose(fid);

