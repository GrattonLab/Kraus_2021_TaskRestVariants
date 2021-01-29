
clear all

%% Uses seed maps to assign contiguous variants and vertices to functional networks for each state, in the same vertices in the other state, and in the same vertices across subjects
%
% This script assigns contiguous variants (and vertices) to functional
% networks based on the similarity between the variant seed map and the
% seed map for each group level network. A uniqueID map is used to identify
% the vertices that correspond to variants. The timeseries of these
% vertices are then correlated to the timeseries of all other vertices
% (greater than a specified geodesic distance away from the variant). The
% seed map for each vertex that composes a variant is then averaged to
% create a mean seed map for that variant. This seed map is binarized (1s 
% indicate a variant exists at a vertex, 0s indicate no variant is present)
% with the spatial locations of the top 5% of correlation values with the 
% seed map. The same is done for the group-average network templates, and
% the overlap between the spatial locations of each variant and all the
% network templates is evaluated. The network template with the highest
% spatial correlation is declared the "winner", and the variant is assigned
% to belonging to that network. If this variant either overlaps spatially
% with the template network it was assigned to (> 50% of its vertices), has
% a spatial correlation value with the template network in the lowest 5% of
% all template matches (.1428 in this sample), or there is a tie between
% the "winning" networks, then the variant is assigned to an "unknown"
% network. This process is also repeated for the individual vertices that
% compose each contiguous variant. Variant network matches are computed for
% variants within their original state (i.e., task or rest), for the same
% vertices between states (i.e., if a variant is observed during rest, the
% same vertices would also be assigned to a variant in the same subject's
% task data), and across subjects (e.g. the variant vertices at rest in
% MSC01 are also assigned to networks in the task data for MSC03). The
% results of these network assignments are then saved out into a structure
% for further analysis. The spatial maps containing the network assignments
% (colored according to the Power networks) can also be saved out as
% CIFTIs.
%
% INPUTS:
%
% -outputdir: the output directory for the structure and the output spatial
% and maps for each state and subject (note that subdirectories are created
% within this directory for the output files)
% -networktemplatepath: the path to the network template files for the
% group average templates and the consensus group network assignments
% -groupnettemplatefile: the filename and extension for the group average
% template network file, which contains the average group seed maps for the
% 14 Power networks
% -consensustemplatefile: the filename and extension for the group average
% network parcellation file, which contains the group network assignments 
% for each vertex to the 14 Power networks
% -geodesicdistfile: filename and extension for a geodesic distance file
% which is a matrix of the disance between every vertex and every other
% vertex. Used to exclude vertices within a minimum distance of the seed
% map being assigned (can be generated for each vertex via:
% wb_command -surface-geodesic-distance)
% -templatepath: the path to a template CIFTI to use for outputting network
% variant assignments on the surface
%
% -minsize: the minimum size threshold (in vertices) that a variant has to
% meet in order to be included in the output data
% -SizeExclude: toggles whether to exclude variants below a minimum size
% (in vertices) threshold (set to 1), or use all variant locations
% regardless of size (set to 0)
% -mindist: the minimum distance threshold (mm) that a vertex must be from
% any vertex in a seed map to be included in the spatial correlation with
% the group network templates
% -variantthreshold: a threshold of uniqueID maps to load for matching 
% variant seed maps to networks
% -spatialtemplatethreshold: the highest percent of spatial correlation
% values with each network to binarize for the spatial correlation between
% the group level template and each individual seed map
% -excludeconsensus: toggles whether to exclude variants that share > 50%
% of their vertices with the group-average network parcellation (set to 1),
% or to assign them to that network regardless of spatial overlap (set to 
% 0)
% -minmatchval: The lowest 5th percentile of the spatial correlation values
% between the group level template and each individual seed map, matches
% below this threshold are assigned to an "unknown" network
% -network_names: a cell array containing the abbreviations for each of the
% Power networks (corresponding to the numbers in the wbcolors variable)
% -wb_colors: a vector of numbers that correspond to the colors of the
% Power networks in connectome workbench
%
% -task_UniqueIDs: reads path to a text file containing the paths to each 
% of the uniqueID files for task data for all sessions. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
% -rest_UniqueIDs: reads path to a text file containing the paths to each 
% of the uniqueID files for rest data for all sessions. The format is 
% (pathtofile subjectID) and the file is space-delimited (see below for 
% more formatting details)
% -task_timeseries: reads path to a text file containing the paths to each 
% of the concatenated timeseries files for task data for all sessions. The 
% format is (pathtofile subjectID) and the file is space-delimited (see 
% below for more formatting details)
% -rest_timeseries: reads path to a text file containing the paths to each 
% of the concatenated timeseries files for rest data for all sessions. The 
% format is (pathtofile subjectID) and the file is space-delimited (see 
% below for more formatting details)
%
% OUTPUTS:
%
% -Network Assignment Maps (optional): A CIFTI file (numvertices x 1) with 
% values that correspond to the network assignment for each variant 
% (using the Power network colors in connectome workbench)
% -Variant_Match_Vertices_Struct: a structure that contains all of the
% network matches for each variant and vertex within state, between state,
% and across subjects. The outputs are:
%  -Subject: the subject ID the variant was defined in
%  -State: the state (task or rest) the variant was defined in
%  -Var_Num: the unique identifier for the parcel corresponding to each
%  variant
%  -Variant_Vertices: the indices of the CIFTI file that correspond to each
%  variant/vertex
%  -Network: the name of the network the variant is assigned to
%  -Network_Other_State: the name of the network the variant is assigned to
%  in the other state (e.g., the network corresponding to the vertices in
%  task if a network is identified in rest)
%  -Net_Assign_Perm: the name of the network the variant is assigned to in
%  another subject (e.g., the network assignment for the vertices of a
%  variant for MSC01 in rest using the timeseries data of MSC03 in task)
%  -Vertex_Assignments: the number corresponding to the network assignment
%  for each individual vertex that composes each variant
%  -Network_Other_State_Vertex: The number corresponding to the network
%  assignment for each individual vertex of a variant in the opposite state
%  (e.g., the network corresponding to the vertices in task if a network is
%  identified in rest)
%  -Net_Assign_Perm_Vertex: the number corresponding to the network
%  assignment for each individual vertex assigned in another subject (e.g., 
%  the network assignment for a vertex for MSC01 in rest using the 
%  timeseries data of MSC03 in task)
%  -Seed_Map_Other_State: the spatial correlation between the mean seed map
%  for a variant in one state versus the seed map for the same vertices in 
%  the opposite state (e.g., the correlation of the mean seed map for the 
%  same vertices in task versus rest)
%  -Seed_Map_Other_State_Perm: the spatial correlation between the mean 
%  seed map for a variant versus the seed map for the vertices in another
%  subject (e.g., the correlation of the mean seed map for a variant in 
%  MSC01 in rest with the same seed map for MSC03 in task)
%
% Written by BK (01-2021)
%


%% Initialize Variables

outputdir = '/your/path/here/';  %% output directory for spatial variant maps and variant structure
networktemplatepath = '/your/path/here/';  %% file path to network templates
groupnettemplatefile = 'WashU120_groupavg_networktemplates.dtseries.nii';  %% filename for group network template
consensustemplatefile = 'WashU120_groupavg_parcellation.dtseries.nii';  %% filename for consensus network template
geodesicdistfile = '/your/path/here/geo_distances.mat';  %% filename for file containing geodesic distances between all vertices
ciftitemplate = '/your/path/here/template_cifti.dtseries.nii';  %% Path where template cifti is located
SizeExclude = 1;  %% Toggles whether to exclude small variants
minsize = 50;  %% Minimum variant size for size exclusion
mindist = 30;  %% minimum geodesic distance (mm) that vertices must be apart in order to be used in seed map
variantthreshold = 5;  %% Threshold of variants to create maps from
spatialtemplatethreshold = 5;  %% percentile of highest correlations to select for spatial overlap between seed and group network templates
excludeconsensus = 1;  %% toggles whether to exclude variants that overlap > 50% spatially with the same network they are assigned to
minmatchval = 0.1428;      %% threshold below which to define variants as belonging to an "unknown" network (set by the lowest %% of Dice correlations with the template across all subjects)

network_names = {'DMN'	'Vis'	'FP'	'Unassigned'	'DAN'	'Unassigned2'	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'	'PMN'	'PON'};  %% The names for the Power functional networks
wb_colors = [1 2 3 5 7 8 9 10 11 12 13 14 15 16];    %% The colors that correspond to the Power functional networks in connectome workbench


% output structure containing network assignments for variants within
% states, between states, and across subjects

Variant_Match_Vertices_Struct = struct('Subject',{},'State',{},'Network',{},'Var_Num',{},'Vertex_Assignments',{},'Network_Other_State',{},'Net_Assign_Perm',{},'Network_Other_State_Vertex',{},'Net_Assign_Perm_Vertex',{},'Seed_Map_Other_State',{},'Seed_Map_Other_State_Perm',{},'Variant_Vertices',{});

% load geodesic distance file
cifti_geo = load([networktemplatepath geodesicdistfile]);
cifti_geo.distances = cifti_geo.distances(1:59412,1:59412);  %% select only surface vertices

%load group template-matrix of 16 networks x XX,XXX vertices
nettemplates = ft_read_cifti_mod([networktemplatepath groupnettemplatefile]);
nettemplates = nettemplates.data(1:59412,:)';

%threshold templates at xth percentile
temp = sort(nettemplates(:),'descend');
templateThresh = temp(round((spatialtemplatethreshold/100) * numel(temp)));
cifti_template_mat_full = nettemplates > templateThresh;
clear temp


% load consensus networks
if excludeconsensus
    
    %%% group-average network assignments (surface only), 59412 x 1
    consensus = ft_read_cifti_mod([networktemplatepath consensustemplatefile]);
   
end

% Create output directories and set output paths

if ~exist([outputdir 'Variant_Network_Assignment_Maps'], 'dir')
    mkdir([outputdir 'Variant_Network_Assignment_Maps']);
end
outputdirmaps = [outputdir 'Variant_Network_Assignment_Maps/'];

% read text files with paths to data

% Unique ID maps
% load uniqueID maps for analyses (see threshold_variant_maps.m for
% additional documentation on uniqueID maps) using a space-delimited text 
% file in the format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% the order of the data files for all of the subjects should be the same in
% all text files

[task_UniqueIDs, subjects, ~] = textread(['/your/path/here/MSC_alltasks_varmaps_' num2str(variantthreshold) '_uniqueID.txt'],'%s%s%s');

[rest_UniqueIDs, ~, ~] = textread(['/your/path/here/MSC_rest_varmaps_' num2str(variantthreshold) '_uniqueID.txt'],'%s%s%s');

% Timeseries files
% load concatenated surface timeseries for analyses 
% (numvertice x numtimepoints)using a space-delimited text file in the 
% format: pathtofile subID
% e.g. filepath/MSC01.dtseries.nii MSC01
% the order of the data files for all of the subjects should be the same in
% all text files

[task_timeseries,~,~] = textread('/your/path/here/MSC_alltask_timeseries.txt','%s%s%s');

[rest_timeseries,~,~] = textread('/your/path/here/MSC_rest_timeseries.txt','%s%s%s');


nfiles = length(rest_UniqueIDs);  %% total number of CIFTI files in text files

for x = 1:nfiles
    
    %% Load Files
    
	subject = subjects{x};
    
    disp(['Processing data for subject ' subject])
        
	cifti_rest_ID = ft_read_cifti_mod(rest_UniqueIDs{x});  % load task and rest data for each subject
	cifti_task_ID = ft_read_cifti_mod(task_UniqueIDs{x});
    
    template_rest = ft_read_cifti_mod(ciftitemplate);  % create templates for saving out spatial maps
    template_rest.data = zeros(size(template_rest.data));
    template_task = template_rest;

    if SizeExclude   % exclude variants smaller than size threshold
        
        cifti_rest_ID.data = variant_size_exclude(cifti_rest_ID,minsize);
        cifti_task_ID.data = variant_size_exclude(cifti_task_ID,minsize);
        
    end
    
	vars_rest = unique(cifti_rest_ID.data(cifti_rest_ID.data>0));  %% get unique variant IDs in each state
   	vars_task = unique(cifti_task_ID.data(cifti_task_ID.data>0));
    
    cifti_rest_timeseries = ft_read_cifti_mod(rest_timeseries{x});  %% load CIFTI timeseries for each state
    cifti_task_timeseries = ft_read_cifti_mod(task_timeseries{x});
    
	%% Match network variants in rest to each network
    
    %%% loop through variants, get avg seedmap, match to network templates %%%
    for var = 1:length(vars_rest)

        inds = find(cifti_rest_ID.data == vars_rest(var));  %% get CIFTI indices of each variant
        
        %%% Exclude vertices too close to seedmap
        
        excludedverts = [];
        
        for vrts = 1:length(inds)
            excludedverts = [excludedverts find(cifti_geo.distances(inds(vrts),:) < mindist)];
        end
        
        finalexclude = unique([inds' excludedverts]);  % vector of vertices to exclude from seedmaps
        
        %%% Assign whole variant to network
        
        varRmat = paircorr_mod(cifti_rest_timeseries.data(inds,:)',cifti_rest_timeseries.data'); % correlate each vertex within variant to all vertices in brain
        
        varRmatGeo = varRmat;
        
        varRmatGeo(:,finalexclude) = []; % exclude variant vertices within minimum geodesic distance of all other vertices
        
        varmean = mean(varRmatGeo,1); %average seedmap of variant
        
        realnetmatch = match_seed_to_template(varmean,cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds);

        if realnetmatch < 99   %% if variant not excluded for overlapping, add it to spatial map
            
            template_rest.data(inds) = realnetmatch;
            
        end
        
        %%% Assign variant to network in other state
        
        varRmat_other = paircorr_mod(cifti_task_timeseries.data(inds,:)',cifti_task_timeseries.data'); %correlate each vertex within variant to all vertices in brain
        
        varRmat_other(:,finalexclude) = [];
        
        varmean_other = mean(varRmat_other,1);  % average seedmap of variant
        
        spat_corr_other = corrcoef(varmean, varmean_other);  % correlate average seedmaps of variants across each state
        
        othernetmatch = match_seed_to_template(varmean_other,cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds);

        %%% Assign networks by vertex
        
        netmatchesverts = [];
        othermatchesverts = [];
        
        for verts = 1:size(varRmatGeo,1)
            
            netmatchesverts = [netmatchesverts match_vertex_to_template(varRmatGeo(verts,:),cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds(verts))];
            othermatchesverts = [othermatchesverts match_vertex_to_template(varRmat_other(verts,:),cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds(verts))];

        end
        
        %%% Do permutations across subjects
        
        NetAssignPerm = {};
        NetAssignPermVerts = {};
        SpatialCorrPerm = [];
        
        permsubs = [1:nfiles];  % get vector with number for each subject file
        permsubs(x) = [];  % leave the current subject out
        
        for y = 1:length(permsubs)
            
            NetAssignPermVertsSub = [];

            cifti_subperm_task_timeseries = ft_read_cifti_mod(task_timeseries{permsubs(y)});
            
            varRmat_perm = paircorr_mod(cifti_subperm_task_timeseries.data(inds,:)',cifti_subperm_task_timeseries.data'); % correlate each vertex within variant to all vertices in other subject's cortex
            
            varRmat_perm(:,finalexclude) = [];
            
            varmean_perm = mean(varRmat_perm,1); %average seedmap of variant
            
            spat_corr_other_perm = corrcoef(varmean, varmean_perm);  % correlate average seedmaps of variants across each state and subject
            
            SpatialCorrPerm = [SpatialCorrPerm; spat_corr_other_perm(1,2)];
            
            NetAssignPerm = [NetAssignPerm; match_seed_to_template(varmean_perm,cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds)];

            if excludeconsensus && NetAssignPerm{y} == 99  %% checks network assignment across subjects versus consensus
                NetAssignPerm{y} = 'Excluded';
            else
                NetAssignPerm{y} = network_names(NetAssignPerm{y});
            end
            
            for verts = 1:size(varRmat_perm,1)
                
                NetAssignPermVertsSub = [NetAssignPermVertsSub; match_vertex_to_template(varRmat_perm(verts,:),cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds(verts))];
                
            end
            
            NetAssignPermVerts = [NetAssignPermVerts NetAssignPermVertsSub];
            
        end
        
        % add values to a structure to save
        
        StructAdd = size(Variant_Match_Vertices_Struct,2)+1;
        
        Variant_Match_Vertices_Struct(StructAdd).Subject = subject;
        Variant_Match_Vertices_Struct(StructAdd).State = 'Rest';
        if excludeconsensus && realnetmatch == 99
            Variant_Match_Vertices_Struct(StructAdd).Network = 'Excluded';
        else
            Variant_Match_Vertices_Struct(StructAdd).Network = network_names(realnetmatch);
        end
        Variant_Match_Vertices_Struct(StructAdd).Var_Num = vars_rest(var);
        Variant_Match_Vertices_Struct(StructAdd).Vertex_Assignments = netmatchesverts;
        if excludeconsensus && othernetmatch == 99
            Variant_Match_Vertices_Struct(StructAdd).Network_Other_State = 'Excluded';
        else
            Variant_Match_Vertices_Struct(StructAdd).Network_Other_State = network_names(othernetmatch);
        end
        Variant_Match_Vertices_Struct(StructAdd).Net_Assign_Perm = NetAssignPerm;
        Variant_Match_Vertices_Struct(StructAdd).Network_Other_State_Vertex = othermatchesverts;
        Variant_Match_Vertices_Struct(StructAdd).Net_Assign_Perm_Vertex = NetAssignPermVerts;
        Variant_Match_Vertices_Struct(StructAdd).Seed_Map_Other_State = spat_corr_other(1,2);
        Variant_Match_Vertices_Struct(StructAdd).Seed_Map_Other_State_Perm = SpatialCorrPerm;
        Variant_Match_Vertices_Struct(StructAdd).Variant_Vertices = inds;
        
    end
    
    %% Do the same for all task variants for this subject
    
    %%% loop through variants, get avg seedmap, match to network templates %%%
    for var = 1:length(vars_task)

        inds = find(cifti_task_ID.data == vars_task(var));
        
        %%% Exclude vertices too close to seedmap
        
        excludedverts = [];
        
        for vrts = 1:length(inds)
            
            excludedverts = [excludedverts find(cifti_geo.distances(inds(vrts),:) < mindist)];
            
        end
        
        finalexclude = unique([inds' excludedverts]);  % vector of vertices to exclude from seedmaps
        
        %%% Assign whole variant to network
        
        varRmat = paircorr_mod(cifti_task_timeseries.data(inds,:)',cifti_task_timeseries.data'); %correlate each vertex within variant to all vertices in brain
        
        varRmatGeo = varRmat;
        
        varRmatGeo(:,finalexclude) = [];
        
        varmean = mean(varRmatGeo,1); %average seedmap of variant
        
        realnetmatch = match_seed_to_template(varmean,cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds);
        
        if realnetmatch < 99   %% if variant not excluded for overlapping, add it to spatial map
            
            template_task.data(inds) = realnetmatch;
            
        end
        
        %%% Assign variant to network in other state
        
        varRmat_other = paircorr_mod(cifti_rest_timeseries.data(inds,:)',cifti_rest_timeseries.data'); %correlate each vertex within variant to all vertices in brain
        
        varRmat_other(:,finalexclude) = [];
        
        varmean_other = mean(varRmat_other,1); %average seedmap of variant
        
        spat_corr_other = corrcoef(varmean, varmean_other);  % correlate average seedmaps of variants across each state
        
        othernetmatch = match_seed_to_template(varmean_other,cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds);
        
        %%% Assign networks by vertex
        
        netmatchesverts = [];
        othermatchesverts = [];
        
        for verts = 1:size(varRmatGeo,1)
            
            netmatchesverts = [netmatchesverts match_vertex_to_template(varRmatGeo(verts,:),cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds(verts))];
            othermatchesverts = [othermatchesverts match_vertex_to_template(varRmat_other(verts,:),cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds(verts))];

        end
        
        %%% Do permutations
        
        NetAssignPerm = {};
        NetAssignPermVerts = {};
        SpatialCorrPerm = [];
        
        
        permsubs = [1:nfiles];
        permsubs(x) = [];
        
        for y = 1:length(permsubs)

            NetAssignPermVertsSub = [];
            
            cifti_subperm_rest_timeseries = ft_read_cifti_mod(rest_timeseries{permsubs(y)});
            
            varRmat_perm = paircorr_mod(cifti_subperm_rest_timeseries.data(inds,:)',cifti_subperm_rest_timeseries.data'); %correlate each vertex within variant to all vertices in other subject's brain
            
            varRmat_perm(:,finalexclude) = [];
            
            varmean_perm = mean(varRmat_perm,1); %average seedmap of variant
            
            spat_corr_other_perm = corrcoef(varmean, varmean_perm);  % correlate average seedmaps of variants across each state and subject
            
            SpatialCorrPerm = [SpatialCorrPerm; spat_corr_other_perm(1,2)];
            
            NetAssignPerm = [NetAssignPerm; network_names(match_seed_to_template(varmean_perm,cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds))];

            if excludeconsensus && NetAssignPerm{y} == 99  %% checks network assignment across subjects versus consensus
                NetAssignPerm{y} = 'Excluded';
            else
                NetAssignPerm{y} = network_names(NetAssignPerm{y});
            end
            
            for verts = 1:size(varRmat_perm,1)
                
                NetAssignPermVertsSub = [NetAssignPermVertsSub; match_vertex_to_template(varRmat_perm(verts,:),cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds(verts))];
                
            end

            NetAssignPermVerts = [NetAssignPermVerts NetAssignPermVertsSub];     
            
        end
        
        % add values to a structure to save

        StructAdd = size(Variant_Match_Vertices_Struct,2)+1;
        
        Variant_Match_Vertices_Struct(StructAdd).Subject = subject;
        Variant_Match_Vertices_Struct(StructAdd).State = 'Task';
        if excludeconsensus && realnetmatch == 99
            Variant_Match_Vertices_Struct(StructAdd).Network = 'Excluded';
        else
            Variant_Match_Vertices_Struct(StructAdd).Network = network_names(realnetmatch);
        end
        Variant_Match_Vertices_Struct(StructAdd).Var_Num = vars_task(var);
        Variant_Match_Vertices_Struct(StructAdd).Vertex_Assignments = netmatchesverts;
        if excludeconsensus && othernetmatch == 99
            Variant_Match_Vertices_Struct(StructAdd).Network_Other_State = 'Excluded';
        else
            Variant_Match_Vertices_Struct(StructAdd).Network_Other_State = network_names(othernetmatch);
        end
        Variant_Match_Vertices_Struct(StructAdd).Net_Assign_Perm = NetAssignPerm;
        Variant_Match_Vertices_Struct(StructAdd).Network_Other_State_Vertex = othermatchesverts;
        Variant_Match_Vertices_Struct(StructAdd).Net_Assign_Perm_Vertex = NetAssignPermVerts;
        Variant_Match_Vertices_Struct(StructAdd).Seed_Map_Other_State = spat_corr_other(1,2);
        Variant_Match_Vertices_Struct(StructAdd).Seed_Map_Other_State_Perm = SpatialCorrPerm;
        Variant_Match_Vertices_Struct(StructAdd).Variant_Vertices = inds;
        
    end
    
    % write out variant spatial maps for each subject by state
    
    rest_fname = create_filename(subject,'Rest',variantthreshold);
    
    ft_write_cifti_mod([outputdirmaps rest_fname],template_rest);
    
    task_fname = create_filename(subject,'Task',variantthreshold);
    
    ft_write_cifti_mod([outputdirmaps task_fname],template_task);

    
end
    

if excludeconsensus
    
    save([outputdir 'Variants_Match_Vertices_Struct_Consensus.mat'], 'Variant_Match_Vertices_Struct');
    
else

    save([outputdir 'Variants_Match_Vertices_Struct.mat'], 'Variant_Match_Vertices_Struct');
    
end




function realwinner = match_seed_to_template(seedmap,cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds)

% matches a seedmap for a variant to the group network template based on the spatial
% overlap between the two templates. the xth highest % of correlations are
% selected for the spatial correlation

temp = sort(seedmap,'descend');
varthresh = temp(round((spatialtemplatethreshold/100) * numel(temp)));
varmean_thresh = seedmap > varthresh;

for net = 1:14
    
    %%% dice calculation b/w variant & each network %%%
    netTemplate = cifti_template_mat_full(net,:);
    
    netTemplate(:,finalexclude) = [];
    
    var_dice_to_template(net,:) = (sum((netTemplate & varmean_thresh),2)*2) ./ sum((netTemplate | varmean_thresh),2);
    
end


%%% pick highest eta value and assign variant to that network %%%
[maxval, winner] = max(var_dice_to_template,[],1);

if maxval > 0 && maxval > minmatchval && length(find(var_dice_to_template==maxval)) < 2   %% set variants with max dice correlation to the template lower than lowest 5% of matches, variants with multiple equally good matches, or variants with no good matches to an "unknown" network
    
    winner = wb_colors(winner);
    
else
    
    winner = 4;
    
end

if excludeconsensus
    
    %%% remove variant if over half its vertices overlap w consensus %%%
    if length(intersect(inds,find(consensus.data==winner))) > (length(inds)*0.5)
        
        realwinner = 99;
        
    else
        
        realwinner = winner;
        
    end
    
else
    
    realwinner = winner;
    
end

end


function realwinner = match_vertex_to_template(seedmap,cifti_template_mat_full,spatialtemplatethreshold,wb_colors,excludeconsensus,consensus,minmatchval,finalexclude,inds)

% matches a seedmap for a vertex to the group network template based on the spatial
% overlap between the two templates. the xth highest % of correlations are
% selected for the spatial correlation

temp = sort(seedmap,'descend');
varthresh = temp(round((spatialtemplatethreshold/100) * numel(temp)));
varmean_thresh = seedmap > varthresh;

for net = 1:14
    
    %%% dice calculation b/w variant & each network %%%
    netTemplate = cifti_template_mat_full(net,:);
    
    netTemplate(:,finalexclude) = [];
    
    var_dice_to_template(net,:) = (sum((netTemplate & varmean_thresh),2)*2) ./ sum((netTemplate | varmean_thresh),2);
    
end


%%% pick highest eta value and assign variant to that network %%%
[maxval, winner] = max(var_dice_to_template,[],1);

if maxval > 0 && maxval > minmatchval && length(find(var_dice_to_template==maxval)) < 2   %% set variants with max dice correlation to the template lower than lowest 5% of matches, variants with multiple equally good matches, or variants with no good matches to an "unknown" network
    
    winner = wb_colors(winner);
    
else
    
    winner = 4;
    
end

if excludeconsensus
    
    %%% remove variant if over half its vertices overlap w consensus %%%
    if consensus.data(inds) == winner
        
        realwinner = 99;
        
    else
        
        realwinner = winner;
        
    end
    
else
    
    realwinner = winner;
    
end

end


function filename = create_filename(sub,task,threshold)

% set output filename for each map

ses = 'All';
grp = 'WashU120';
desc = 'matched';
thresh_str = num2str(threshold);


filename = sprintf('sub-%s_task-%s_ses-%s_grp-%s_desc-%s_thresh-%s_NetworkAssignmentMap.dtseries.nii',sub,task,ses,grp,desc,thresh_str);

end


