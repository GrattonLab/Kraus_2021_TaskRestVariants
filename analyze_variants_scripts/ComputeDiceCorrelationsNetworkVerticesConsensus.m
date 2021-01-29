
clear all

%% Plot data from network matching structure for network variants, variant vertices, and seed map correlations
%
% This script loads in the output structure from AssignNetworksByVertex.m
% and summarizes the variant network assignment outcomes. This script plots
% pie charts to summarize the network assignments for variants and vertices
% and compares network assignments between states and across subjects. 
%
%
% INPUTS:
%
% -outdir: output directory for plots summarizing network assignment
% -structdir: the directory containing the structure with the network
% assignment data
% -threshold: the threshold used to binarize the variant maps used for
% network assignment
% -network_names: a cell array containing the abbreviations for each of the
% Power networks (corresponding to the numbers in the wbcolors variable)
% -rgb_colors: a 17x3 matrix containing all of the MATLAB colors for the
% Power networks (and an additional color for the consensus network)
% -Variants_Match_Vertices_Struct_Consensus: a structure containing network
% assignment outcomes (see AssignNetworksByVertex.m for additional
% documentation)
%
% OUTPUTS:
%
% -pie charts: plots summarizing the total number of variants and vertices
% assigned to each network in task and rest states
% -seed maps: plot summarizing the mean correlation between variant seed
% maps between states and across subjects
% -network assignment: plots summarizing the network assignments for
% variants and individual vertices between states and across subjects
% -network assignment and location: plot showing the likelihood of variant
% vertices both spatially overlapping between states and being assigned to
% the same network versus the same vertice across subjects
%
% Written by BK (01-2021)
%

%% Initialize Variables

outdir = '/your/path/here/';  %% output directory
structdir = '/your/path/here/';  %% directory containing structure

threshold = 5;      %% Threshold of matched data to load

% load struct with network matching data

load([structdir 'Variants_Match_Vertices_Struct_Consensus.mat'])

% set up network matching information for script

network_names = {'DMN'	'Vis'	'FP'	'Unassigned'	'DAN'	'Unassigned2'	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'	'PMN'	'PON'   'Consensus'};    %% The colors that correspond to the Power functional networks in connectome workbench
rgb_colors = [1 0 0;    %% The colors that correspond to the Power functional networks in MATLAB
              0 0 .6;
              1 1 0;
              .67 .67 .67;
              0 .8 0;
              .67 .67 .67;
              0 .6 .6;
              0 0 0;
              .3 0 .6;
              .2 1 1;
              1 .5 0;
              .6 .2 1;
              .2 1 .2;
              0 .2 .4;
              0 0 1;
              .8 .8 .6;
              .4 0 0];
          
%% Get unique subject IDs in structure

subs = {};

for x = 1:size(Variant_Match_Vertices_Struct,2)    
    
    if ismember(Variant_Match_Vertices_Struct(x).Subject,subs) == 0
        
        subs = [subs; Variant_Match_Vertices_Struct(x).Subject];
        
    end
end


%% Loop through variants in structure and calculate spatial overlap/network assignment metrics

% store total number of vertices and variants assigned to each network

totalverts = zeros(1,17);
restverts = zeros(1,17);
taskverts = zeros(1,17);

allnets = zeros(1,16);
restnets = zeros(1,16);
tasknets = zeros(1,16);

%%% Network variant map assignment

dicecorrsallvarsseed = [];
dicecorrsallvarspermseed = [];
dicecorrsrestseed = [];
dicecorrsrestpermseed = [];
dicecorrstaskseed = [];
dicecorrstaskpermseed = [];

%%% Vertex-wise assignment

dicecorrsallvars = [];
dicecorrsallvarsperm = [];
dicecorrsrest = [];
dicecorrsrestperm = [];
dicecorrstask = [];
dicecorrstaskperm = [];

%%% Vertex-wise assignment full match

dicecorrsallvarsfull = [];
dicecorrsallvarspermfull = [];
dicecorrsrestfull = [];
dicecorrsrestpermfull = [];
dicecorrstaskfull = [];
dicecorrstaskpermfull = [];

%%% Network variant seed map spatial correlations

spatcorrsallvarsseed = [];
spatcorrsallvarspermseed = [];
spatcorrsrestseed = [];
spatcorrsrestpermseed = [];
spatcorrstaskseed = [];
spatcorrstaskpermseed = [];


for subjects = 1:numel(subs)  %% loops through subjects found in structure

    % variables to store subject-level values

    dicecorrsallvarsseedsub = [];
    dicecorrsallvarspermseedsub = [];
    dicecorrsrestseedsub = [];
    dicecorrsrestpermseedsub = [];
    dicecorrstaskseedsub = [];
    dicecorrstaskpermseedsub = [];
    
    dicecorrsallvarssub = [];
    dicecorrsallvarspermsub = [];
    dicecorrsrestsub = [];
    dicecorrsrestpermsub = [];
    dicecorrstasksub = [];
    dicecorrstaskpermsub = [];
    
    dicecorrsallvarsfullsub = [];
    dicecorrsallvarspermfullsub = [];
    dicecorrsrestfullsub = [];
    dicecorrsrestpermfullsub = [];
    dicecorrstaskfullsub = [];
    dicecorrstaskpermfullsub = [];
    
    spatcorrsallvarsseedsub = [];
    spatcorrsallvarspermseedsub = [];
    spatcorrsrestseedsub = [];
    spatcorrsrestpermseedsub = [];
    spatcorrstaskseedsub = [];
    spatcorrstaskpermseedsub = [];

    for vars = 1:size(Variant_Match_Vertices_Struct,2)  %% for each variant
        
        %% Create binary vectors for dice correlations and gather other variables for plotting
        
        %% Get variant count for each state
        
        if strcmp(Variant_Match_Vertices_Struct(vars).Subject, subs(subjects))  %% if variant is from current subject

            if ~strcmp(Variant_Match_Vertices_Struct(vars).Network, 'Excluded')   %% if variant was not excluded for overlapping with the consensus template
                
                allnets(strcmp(Variant_Match_Vertices_Struct(vars).Network,network_names)) = allnets(strcmp(Variant_Match_Vertices_Struct(vars).Network,network_names))+1;  %% add variant to network assignment count for both states combined
                
                spatcorrsallvarsseedsub = [spatcorrsallvarsseedsub Variant_Match_Vertices_Struct(vars).Seed_Map_Other_State];  %% get spatial correlation value for seed maps between states for each variant
                
                spatcorrsallvarspermseedsub = [spatcorrsallvarspermseedsub; Variant_Match_Vertices_Struct(vars).Seed_Map_Other_State_Perm'];  %% get spatial correlation value for seed maps across subjects for each variant
                
                if strcmp(Variant_Match_Vertices_Struct(vars).State, 'Rest')  %% assign values for variants in rest states
                    
                    restnets(strcmp(Variant_Match_Vertices_Struct(vars).Network,network_names)) = restnets(strcmp(Variant_Match_Vertices_Struct(vars).Network,network_names))+1;  %% add variant to network assignment count for rest state
                    
                    spatcorrsrestseedsub = [spatcorrsrestseedsub Variant_Match_Vertices_Struct(vars).Seed_Map_Other_State];  %% get spatial correlation value for seed maps between states for each variant
                    
                    spatcorrsrestpermseedsub = [spatcorrsrestpermseedsub; Variant_Match_Vertices_Struct(vars).Seed_Map_Other_State_Perm'];  %% get spatial correlation value for seed maps across subjects for each variant
                    
                elseif strcmp(Variant_Match_Vertices_Struct(vars).State, 'Task')  %% assign values for variants in task states
                    
                    tasknets(strcmp(Variant_Match_Vertices_Struct(vars).Network,network_names)) = tasknets(strcmp(Variant_Match_Vertices_Struct(vars).Network,network_names))+1;  %% add variant to network assignment count for task state
                
                    spatcorrstaskseedsub = [spatcorrstaskseedsub Variant_Match_Vertices_Struct(vars).Seed_Map_Other_State];  %% get spatial correlation value for seed maps between states for each variant
                    
                    spatcorrstaskpermseedsub = [spatcorrstaskpermseedsub; Variant_Match_Vertices_Struct(vars).Seed_Map_Other_State_Perm'];  %% get spatial correlation value for seed maps across subjects for each variant
                    
                end
            end
            
            %% Get vertex count for each state

            for count = 1:length(Variant_Match_Vertices_Struct(vars).Vertex_Assignments)   %% for each vertex in each variant
                
                if Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count) == 99   %% if vertex excluded for overlapping with consensus network

                    %totalverts(17) = totalverts(17)+1;  %% add 1 to consensus count
                    
                    if strcmp(Variant_Match_Vertices_Struct(vars).State, 'Rest')  %% assign values for variants in rest states
                        
                        %restverts(17) = restverts(17)+1;  %% add 1 to consensus count

                    elseif strcmp(Variant_Match_Vertices_Struct(vars).State, 'Task')  %% assign values for variants in task states
                        
                        %taskverts(17) = taskverts(17)+1;  %% add 1 to consensus count

                    end
                else
                    
                    totalverts(Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count)) = totalverts(Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count))+1;  %% else, add 1 to count for that network
                
                    if strcmp(Variant_Match_Vertices_Struct(vars).State, 'Rest')  %% assign values for variants in rest states
                        
                        restverts(Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count)) = restverts(Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count))+1;  %% else, add 1 to count for that network

                    elseif strcmp(Variant_Match_Vertices_Struct(vars).State, 'Task')  %% assign values for variants in task states
                        
                        taskverts(Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count)) = taskverts(Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count))+1;  %% else, add 1 to count for that network

                    end
                end
            end
            
            %% Check variant network (and location) matches between states and across subjects
            
         	betweenstatevariantvec = between_state_variant_assignment(Variant_Match_Vertices_Struct,vars);  %% find binary correlation for variant network assignment across states
            
            dicecorrsallvarsseedsub = [dicecorrsallvarsseedsub; betweenstatevariantvec];  %% adds correlations for each subject to matrix
            
            betweenstatevertexvec = between_state_vertex_assignment(Variant_Match_Vertices_Struct,vars);  %% find binary correlation for vertex network assignment across states
            
         	dicecorrsallvarssub = [dicecorrsallvarssub; betweenstatevertexvec];  %% adds correlations for each subject to matrix
            
            acrosssubjectvariantvec = across_subject_variant_assignment(Variant_Match_Vertices_Struct,vars);  %% check if the same vertices for each network variant were assigned to the same network across subjects
            
            dicecorrsallvarspermseedsub = [dicecorrsallvarspermseedsub; acrosssubjectvariantvec];  %% adds correlations for each subject to matrix
            
          	acrosssubjectvertexvec = across_subject_vertex_assignment(Variant_Match_Vertices_Struct,vars);  %% check if each vertex was assigned to the same network in each subject
            
            dicecorrsallvarspermsub = [dicecorrsallvarspermsub; acrosssubjectvertexvec];  %% adds correlations for each subject to matrix
            
            betweenstatefullvec = between_state_vertex_location_assignment(Variant_Match_Vertices_Struct,vars,subs,subjects);  %% match vertex location and network assignment between states
            
            dicecorrsallvarsfullsub = [dicecorrsallvarsfullsub; betweenstatefullvec];  %% adds correlations for each subject to matrix
            
            acrosssubjectsfullvec = across_subject_vertex_location_assignment(Variant_Match_Vertices_Struct,vars,subs,subjects);  %% match vertex location and network assignment across subjects  

            dicecorrsallvarspermfullsub = [dicecorrsallvarspermfullsub; acrosssubjectsfullvec];  %% adds correlations for each subject to matrix
            
            if strcmp(Variant_Match_Vertices_Struct(vars).State, 'Rest')  %% assign values for variants in rest states
                
                dicecorrsrestseedsub = [dicecorrsrestseedsub; betweenstatevariantvec];  %% add between state variant dice correlation values to rest
                dicecorrsrestsub = [dicecorrsrestsub; betweenstatevertexvec];  %% add vertex between state dice correlation values to rest
                dicecorrsrestpermsub = [dicecorrsrestpermsub; acrosssubjectvertexvec];  %% add vertex across subject dice correlation values to rest
                dicecorrsrestpermseedsub = [dicecorrsrestpermseedsub; acrosssubjectvariantvec];  %% add variant across subject dice correlation values to rest
                dicecorrsrestfullsub = [dicecorrsrestfullsub; betweenstatefullvec];  %% add vertex between state location and network correlation values to rest
                dicecorrsrestpermfullsub = [dicecorrsrestpermfullsub; acrosssubjectsfullvec];  %% add vertex across subject location and network correlation values to rest
                
            elseif strcmp(Variant_Match_Vertices_Struct(vars).State, 'Task')  %% assign values for variants in task states
                
                dicecorrstaskseedsub = [dicecorrstaskseedsub; betweenstatevariantvec];  %% add between state variant dice correlation values to task
                dicecorrstasksub = [dicecorrstasksub; betweenstatevertexvec];  %% add vertex between state dice correlation values to task
                dicecorrstaskpermsub = [dicecorrstaskpermsub; acrosssubjectvertexvec];  %% add vertex across subject dice correlation values to task
                dicecorrstaskpermseedsub = [dicecorrstaskpermseedsub; acrosssubjectvariantvec];  %% add variant across subject dice correlation values to task
                dicecorrstaskfullsub = [dicecorrstaskfullsub; betweenstatefullvec];  %% add vertex between state location and network correlation values to task
                dicecorrstaskpermfullsub = [dicecorrstaskpermfullsub; acrosssubjectsfullvec];  %% add vertex across subject location and network correlation values to task
                
            end
            
        end
    end
    
    %% Get summary statistics for each subject
    
    % seed map correlations between states
    
    spatcorrsallvarsseed = [spatcorrsallvarsseed; mean(spatcorrsallvarsseedsub)];  %% calculate the mean of the spatial correlations for the seedmaps between states within each subject
    spatcorrsrestseed = [spatcorrsrestseed; mean(spatcorrsrestseedsub)];
    spatcorrstaskseed = [spatcorrstaskseed; mean(spatcorrstaskseedsub)];

    % seed map correlations across subjects
    
   	spatcorrsallvarspermseed = [spatcorrsallvarspermseed; mean(spatcorrsallvarspermseedsub,1)];  %% calculate the mean of the spatial correlations for the seedmaps between states across subjects
    spatcorrsrestpermseed = [spatcorrsrestpermseed; mean(spatcorrsrestpermseedsub,1)];
    spatcorrstaskpermseed = [spatcorrstaskpermseed; mean(spatcorrstaskpermseedsub,1)];
    
    % between state variant vertex assignment

    dicecorrsallvars = [dicecorrsallvars; dice_coefficient_mod(dicecorrsallvarssub(:,1),dicecorrsallvarssub(:,2))];  % dice correlation between network assignment of vertices in current state and other state for current subject
    dicecorrsrest = [dicecorrsrest; dice_coefficient_mod(dicecorrsrestsub(:,1),dicecorrsrestsub(:,2))];
    dicecorrstask = [dicecorrstask; dice_coefficient_mod(dicecorrstasksub(:,1),dicecorrstasksub(:,2))];
    
    % between state network variant asssignment
    
    dicecorrsallvarsseed = [dicecorrsallvarsseed; dice_coefficient_mod(dicecorrsallvarsseedsub(:,1),dicecorrsallvarsseedsub(:,2))];  % dice correlation between network assignment of variants in current state and other state for current subject
    dicecorrsrestseed = [dicecorrsrestseed; dice_coefficient_mod(dicecorrsrestseedsub(:,1),dicecorrsrestseedsub(:,2))];
    dicecorrstaskseed = [dicecorrstaskseed; dice_coefficient_mod(dicecorrstaskseedsub(:,1),dicecorrstaskseedsub(:,2))];
    
    % between state variant vertex location and network asssignment
    
    dicecorrsallvarsfull = [dicecorrsallvarsfull; dice_coefficient_mod(dicecorrsallvarsfullsub(:,1),dicecorrsallvarsfullsub(:,2))];  % dice correlation between network assignment of vertices and their location in current state and other state for current subject
    dicecorrsrestfull = [dicecorrsrestfull; dice_coefficient_mod(dicecorrsrestfullsub(:,1),dicecorrsrestfullsub(:,2))];
    dicecorrstaskfull = [dicecorrstaskfull; dice_coefficient_mod(dicecorrstaskfullsub(:,1),dicecorrstaskfullsub(:,2))];
    
    % calculate network variant assignments across subjects
    
    dicecorrstempseedsub = [];  % set up temporary variables for looping through subjects
    dicecorrstemprestseedsub = [];
    dicecorrstemptaskseedsub = [];
    
    for perms = 1:numel(Variant_Match_Vertices_Struct(vars).Net_Assign_Perm)
        
        dicecorrstempseedsub = [dicecorrstempseedsub dice_coefficient_mod(dicecorrsallvarsseedsub(:,1),dicecorrsallvarspermseedsub(:,perms))];  % for each subject and all other subjects, compare the network assignments for the spatial locations of each variant in the current subject
        dicecorrstemprestseedsub = [dicecorrstemprestseedsub dice_coefficient_mod(dicecorrsrestseedsub(:,1),dicecorrsrestpermseedsub(:,perms))];
        dicecorrstemptaskseedsub = [dicecorrstemptaskseedsub dice_coefficient_mod(dicecorrstaskseedsub(:,1),dicecorrstaskpermseedsub(:,perms))];
        
    end
    
    dicecorrsallvarspermseed = [dicecorrsallvarspermseed; dicecorrstempseedsub];  %% add each subject's data to output
    dicecorrsrestpermseed = [dicecorrsrestpermseed; dicecorrstemprestseedsub];
    dicecorrstaskpermseed = [dicecorrstaskpermseed; dicecorrstemptaskseedsub];
    
    % calculate variant vertex assignments across subjects
    
    dicecorrsallvarspermtemp = [];  % set up temporary variables for looping through subjects
    dicecorrsrestpermtemp = [];
    dicecorrstaskpermtemp = [];
    
    for perms = 1:size(dicecorrsallvarspermsub,2)
        
        dicecorrsallvarspermtemp = [dicecorrsallvarspermtemp dice_coefficient_mod(dicecorrsallvarssub(:,1),dicecorrsallvarspermsub(:,perms))];  % for each subject and all other subjects, compare the network assignments for the each vertex composing a variant in the current subject
        dicecorrsrestpermtemp = [dicecorrsrestpermtemp dice_coefficient_mod(dicecorrsrestsub(:,1),dicecorrsrestpermsub(:,perms))];
        dicecorrstaskpermtemp = [dicecorrstaskpermtemp dice_coefficient_mod(dicecorrstasksub(:,1),dicecorrstaskpermsub(:,perms))];
        
    end
    
    dicecorrsallvarsperm = [dicecorrsallvarsperm; dicecorrsallvarspermtemp];  %% add each subject's data to output
    dicecorrsrestperm = [dicecorrsrestperm; dicecorrsrestpermtemp];
    dicecorrstaskperm = [dicecorrstaskperm; dicecorrstaskpermtemp];
    
    % calculate variant vertex location and assignments across subjects
    
    dicecorrsallvarspermfulltemp = [];  % set up temporary variables for looping through subjects
    dicecorrsrestpermfulltemp = [];
    dicecorrstaskpermfulltemp = [];
    
    for perms = 1:size(dicecorrsallvarspermfullsub,2)
        
        dicecorrsallvarspermfulltemp = [dicecorrsallvarspermfulltemp dice_coefficient_mod(dicecorrsallvarsfullsub(:,1),dicecorrsallvarspermfullsub(:,perms))];  % for each subject and all other subjects, compare the network assignments and spatial locations for the each vertex composing a variant in the current subject
        dicecorrsrestpermfulltemp = [dicecorrsrestpermfulltemp dice_coefficient_mod(dicecorrsrestfullsub(:,1),dicecorrsrestpermfullsub(:,perms))];
        dicecorrstaskpermfulltemp = [dicecorrstaskpermfulltemp dice_coefficient_mod(dicecorrstaskfullsub(:,1),dicecorrstaskpermfullsub(:,perms))];
        
    end
    
    dicecorrsallvarspermfull = [dicecorrsallvarspermfull; dicecorrsallvarspermfulltemp];  %% add each subject's data to output
    dicecorrsrestpermfull = [dicecorrsrestpermfull; dicecorrsrestpermfulltemp];
    dicecorrstaskpermfull = [dicecorrstaskpermfull; dicecorrstaskpermfulltemp];
    
end


%% Plot Data

% Scatterplot for variant assignment across states and subjects

jitterAmount = 0.1;
jitterValuesX = 2*(rand(size(dicecorrsallvarsseed))-0.5)*jitterAmount;   % +/-jitterAmount max

figure;
scatterpos = [1:3];
scatter(scatterpos+jitterValuesX(9), [dicecorrsallvarsseed(9) dicecorrsrestseed(9) dicecorrstaskseed(9)], 150, 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
xlim([0 4]);
ylim([0 1]);
hold on
scatter(scatterpos+jitterValuesX(8), [dicecorrsallvarsseed(8) dicecorrsrestseed(8) dicecorrstaskseed(8)], 150, 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
hold on
scatter(scatterpos+jitterValuesX(7), [dicecorrsallvarsseed(7) dicecorrsrestseed(7) dicecorrstaskseed(7)], 150, 'filled', 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(6), [dicecorrsallvarsseed(6) dicecorrsrestseed(6) dicecorrstaskseed(6)], 150, 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter(scatterpos+jitterValuesX(5), [dicecorrsallvarsseed(5) dicecorrsrestseed(5) dicecorrstaskseed(5)], 150, 'filled', 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(4), [dicecorrsallvarsseed(4) dicecorrsrestseed(4) dicecorrstaskseed(4)], 150, 'filled', 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter(scatterpos+jitterValuesX(3), [dicecorrsallvarsseed(3) dicecorrsrestseed(3) dicecorrstaskseed(3)], 150, 'filled', 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter(scatterpos+jitterValuesX(2), [dicecorrsallvarsseed(2) dicecorrsrestseed(2) dicecorrstaskseed(2)], 150, 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter(scatterpos+jitterValuesX(1), [dicecorrsallvarsseed(1) dicecorrsrestseed(1) dicecorrstaskseed(1)], 150, 'filled', 'MarkerFaceColor', [0, 0, 0]);
hold on
scatter([repmat(scatterpos(1),1,size(dicecorrsallvarspermseed,1)*size(dicecorrsallvarspermseed,2)) repmat(scatterpos(2),1,size(dicecorrsrestpermseed,1)*size(dicecorrsrestpermseed,2)) repmat(scatterpos(3),1,size(dicecorrstaskpermseed,1)*size(dicecorrstaskpermseed,2))], [reshape(dicecorrsallvarspermseed,1,size(dicecorrsallvarspermseed,1)*size(dicecorrsallvarspermseed,2)) reshape(dicecorrsrestpermseed,1,size(dicecorrsrestpermseed,1)*size(dicecorrsrestpermseed,2)) reshape(dicecorrstaskpermseed,1,size(dicecorrstaskpermseed,1)*size(dicecorrstaskpermseed,2))], 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);


set(gca,'xtick',[scatterpos(1) scatterpos(2) scatterpos(3)])
set(gca,'xticklabel',{'All Variants', 'Rest Variants', 'Task Variants'}, 'FontSize',18)
ylabel('Dice Correlation');
title('Network Assignment Variants', 'FontSize',18)

m = findobj(gca,'Type','scatter');
[l, hobj, ~, ~] = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
l.FontSize = 18;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',18);
ax = gca;
ax.FontSize = 18;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.35, 0.35, 0.5]);

saveas(gcf,[outdir 'Variant_Dice_Correlations_' num2str(threshold) '_Percent.tiff'],'tiff');

close gcf



% Scatterplot for variant vertex assignment across states and subjects

jitterAmount = 0.1;
jitterValuesX = 2*(rand(size(dicecorrsallvars))-0.5)*jitterAmount;   % +/-jitterAmount max

figure;
scatterpos = [1:3];
scatter(scatterpos+jitterValuesX(9), [dicecorrsallvars(9) dicecorrsrest(9) dicecorrstask(9)], 150, 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
xlim([0 4]);
ylim([0 1]);
hold on
scatter(scatterpos+jitterValuesX(8), [dicecorrsallvars(8) dicecorrsrest(8) dicecorrstask(8)], 150, 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
hold on
scatter(scatterpos+jitterValuesX(7), [dicecorrsallvars(7) dicecorrsrest(7) dicecorrstask(7)], 150, 'filled', 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(6), [dicecorrsallvars(6) dicecorrsrest(6) dicecorrstask(6)], 150, 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter(scatterpos+jitterValuesX(5), [dicecorrsallvars(5) dicecorrsrest(5) dicecorrstask(5)], 150, 'filled', 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(4), [dicecorrsallvars(4) dicecorrsrest(4) dicecorrstask(4)], 150, 'filled', 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter(scatterpos+jitterValuesX(3), [dicecorrsallvars(3) dicecorrsrest(3) dicecorrstask(3)], 150, 'filled', 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter(scatterpos+jitterValuesX(2), [dicecorrsallvars(2) dicecorrsrest(2) dicecorrstask(2)], 150, 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter(scatterpos+jitterValuesX(1), [dicecorrsallvars(1) dicecorrsrest(1) dicecorrstask(1)], 150, 'filled', 'MarkerFaceColor', [0, 0, 0]);
hold on
scatter([repmat(scatterpos(1),1,size(dicecorrsallvarsperm,1)*size(dicecorrsallvarsperm,2)) repmat(scatterpos(2),1,size(dicecorrsrestperm,1)*size(dicecorrsrestperm,2)) repmat(scatterpos(3),1,size(dicecorrstaskperm,1)*size(dicecorrstaskperm,2))], [reshape(dicecorrsallvarsperm,1,size(dicecorrsallvarsperm,1)*size(dicecorrsallvarsperm,2)) reshape(dicecorrsrestperm,1,size(dicecorrsrestperm,1)*size(dicecorrsrestperm,2)) reshape(dicecorrstaskperm,1,size(dicecorrstaskperm,1)*size(dicecorrstaskperm,2))], 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);


set(gca,'xtick',[scatterpos(1) scatterpos(2) scatterpos(3)])
set(gca,'xticklabel',{'All Variants', 'Rest Variants', 'Task Variants'}, 'FontSize',18)
ylabel('Dice Correlation');
title('Network Assignment Vertices', 'FontSize',18)

m = findobj(gca,'Type','scatter');
[l, hobj, ~, ~] = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
l.FontSize = 18;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',18);
ax = gca;
ax.FontSize = 18;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.35, 0.35, 0.5]);

saveas(gcf,[outdir 'Vertex_Dice_Correlations_' num2str(threshold) '_Percent.tiff'],'tiff');

close gcf



% Scatterplot for variant vertex assignment and location across states and subjects

jitterAmount = 0.1;
jitterValuesX = 2*(rand(size(dicecorrsallvarsfull))-0.5)*jitterAmount;   % +/-jitterAmount max

figure;
scatterpos = [1:3];
scatter(scatterpos+jitterValuesX(9), [dicecorrsallvarsfull(9) dicecorrsrestfull(9) dicecorrstaskfull(9)], 150, 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
xlim([0 4]);
ylim([0 1]);
hold on
scatter(scatterpos+jitterValuesX(8), [dicecorrsallvarsfull(8) dicecorrsrestfull(8) dicecorrstaskfull(8)], 150, 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
hold on
scatter(scatterpos+jitterValuesX(7), [dicecorrsallvarsfull(7) dicecorrsrestfull(7) dicecorrstaskfull(7)], 150, 'filled', 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(6), [dicecorrsallvarsfull(6) dicecorrsrestfull(6) dicecorrstaskfull(6)], 150, 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter(scatterpos+jitterValuesX(5), [dicecorrsallvarsfull(5) dicecorrsrestfull(5) dicecorrstaskfull(5)], 150, 'filled', 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(4), [dicecorrsallvarsfull(4) dicecorrsrestfull(4) dicecorrstaskfull(4)], 150, 'filled', 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter(scatterpos+jitterValuesX(3), [dicecorrsallvarsfull(3) dicecorrsrestfull(3) dicecorrstaskfull(3)], 150, 'filled', 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter(scatterpos+jitterValuesX(2), [dicecorrsallvarsfull(2) dicecorrsrestfull(2) dicecorrstaskfull(2)], 150, 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter(scatterpos+jitterValuesX(1), [dicecorrsallvarsfull(1) dicecorrsrestfull(1) dicecorrstaskfull(1)], 150, 'filled', 'MarkerFaceColor', [0, 0, 0]);
hold on
scatter([repmat(scatterpos(1),1,size(dicecorrsallvarspermfull,1)*size(dicecorrsallvarspermfull,2)) repmat(scatterpos(2),1,size(dicecorrsrestpermfull,1)*size(dicecorrsrestpermfull,2)) repmat(scatterpos(3),1,size(dicecorrstaskpermfull,1)*size(dicecorrstaskpermfull,2))], [reshape(dicecorrsallvarspermfull,1,size(dicecorrsallvarspermfull,1)*size(dicecorrsallvarspermfull,2)) reshape(dicecorrsrestpermfull,1,size(dicecorrsrestpermfull,1)*size(dicecorrsrestpermfull,2)) reshape(dicecorrstaskpermfull,1,size(dicecorrstaskpermfull,1)*size(dicecorrstaskpermfull,2))], 150, 'jitter','on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [0, 0, 0]);


set(gca,'xtick',[scatterpos(1) scatterpos(2) scatterpos(3)])
set(gca,'xticklabel',{'All Variants', 'Rest Variants', 'Task Variants'}, 'FontSize',18)
ylabel('Dice Correlation');
title('Location and Network Assignment Vertices', 'FontSize',18)

m = findobj(gca,'Type','scatter');
[l, hobj, ~, ~] = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
l.FontSize = 18;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',18);
ax = gca;
ax.FontSize = 18;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.35, 0.35, 0.5]);

saveas(gcf,[outdir 'NetworkAndLocation_Dice_Correlations_' num2str(threshold) '_Percent.tiff'],'tiff');

close gcf


% Scatterplot for variant seed map correlations

jitterAmount = 0.2;
jitterValuesX = 2*(rand(size(spatcorrsallvarsseed))-0.5)*jitterAmount;   % +/-jitterAmount max

figure;
scatterpos = [1:3];
scatter(scatterpos+jitterValuesX(9), [spatcorrsallvarsseed(9) spatcorrsrestseed(9) spatcorrstaskseed(9)], 150, 'filled', 'MarkerFaceColor', [1, 0.5, 0]);
xlim([0 4]);
ylim([-.4 1]);
hold on
scatter(scatterpos+jitterValuesX(8), [spatcorrsallvarsseed(8) spatcorrsrestseed(8) spatcorrstaskseed(8)], 150, 'filled', 'MarkerFaceColor', [0, 0.6, 0.6]);
hold on
scatter(scatterpos+jitterValuesX(7), [spatcorrsallvarsseed(7) spatcorrsrestseed(7) spatcorrstaskseed(7)], 150, 'filled', 'MarkerFaceColor', [1, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(6), [spatcorrsallvarsseed(6) spatcorrsrestseed(6) spatcorrstaskseed(6)], 150, 'filled', 'MarkerFaceColor', [0.2, 1, 1]);
hold on
scatter(scatterpos+jitterValuesX(5), [spatcorrsallvarsseed(5) spatcorrsrestseed(5) spatcorrstaskseed(5)], 150, 'filled', 'MarkerFaceColor', [0, 0, 1]);
hold on
scatter(scatterpos+jitterValuesX(4), [spatcorrsallvarsseed(4) spatcorrsrestseed(4) spatcorrstaskseed(4)], 150, 'filled', 'MarkerFaceColor', [1, 0, 0]);
hold on
scatter(scatterpos+jitterValuesX(3), [spatcorrsallvarsseed(3) spatcorrsrestseed(3) spatcorrstaskseed(3)], 150, 'filled', 'MarkerFaceColor', [0, 1, 0]);
hold on
scatter(scatterpos+jitterValuesX(2), [spatcorrsallvarsseed(2) spatcorrsrestseed(2) spatcorrstaskseed(2)], 150, 'filled', 'MarkerFaceColor', [0.9, 0.9, 0]);
hold on
scatter(scatterpos+jitterValuesX(1), [spatcorrsallvarsseed(1) spatcorrsrestseed(1) spatcorrstaskseed(1)], 150, 'filled', 'MarkerFaceColor', [0, 0, 0]);
hold on
scatter([repmat(scatterpos(1),1,size(spatcorrsallvarspermseed,1)*size(spatcorrsallvarspermseed,2)) repmat(scatterpos(2),1,size(spatcorrsrestpermseed,1)*size(spatcorrsrestpermseed,2)) repmat(scatterpos(3),1,size(spatcorrstaskpermseed,1)*size(spatcorrstaskpermseed,2))], [reshape(spatcorrsallvarspermseed,1,size(spatcorrsallvarspermseed,1)*size(spatcorrsallvarspermseed,2)) reshape(spatcorrsrestpermseed,1,size(spatcorrsrestpermseed,1)*size(spatcorrsrestpermseed,2)) reshape(spatcorrstaskpermseed,1,size(spatcorrstaskpermseed,1)*size(spatcorrstaskpermseed,2))], 150, 'jitter','on', 'jitterAmount', 0.2, 'MarkerEdgeColor', [0, 0, 0]);



set(gca,'xtick',[scatterpos(1) scatterpos(2) scatterpos(3)])
set(gca,'xticklabel',{'All Variants', 'Rest Variants', 'Task Variants'}, 'FontSize',18)
ylabel('Pearson Correlation (r)');
title('Variant Seed Map Correlations')

m = findobj(gca,'Type','scatter');
[l, hobj, hout, mout] = legend(m(2:10), 'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10', 'Location', 'NorthEast');
l.FontSize = 18;
s = findobj(hobj,'type','patch');
set(s,'MarkerSize',10);
x = findobj(hobj,'type','text');
set(x,'FontSize',18);
ax = gca;
ax.FontSize = 18;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.4, 0.5]);

saveas(gcf,[outdir 'Variant_Seedmap_Correlations_' num2str(threshold) '_Percent.tiff'],'tiff');

close gcf


%%% Pie Charts


% All variants pie chart vertices

networksexist = find(totalverts > 0);
network_legend = network_names(networksexist);
color_legend = rgb_colors(networksexist,:);

filename = ['AllVariants_PieChart_' num2str(threshold) 'pct_AllSubs_Vertices.jpg'];


h = pie(totalverts(networksexist));
delete(findobj(h,'Type','text'))

for k = 1:length(h)/2
    
    set(h(k*2-1), 'FaceColor', color_legend(k,:));
    
end

legend(network_legend,'Orientation','vertical','Location','eastoutside', 'FontSize',20)

print(gcf,[outdir filename],'-dpng','-r300');

close gcf



% Rest variants vertices pie chart

networksexist = find(restverts > 0);
network_legend = network_names(networksexist);
color_legend = rgb_colors(networksexist,:);

filename = ['RestVariants_PieChart_' num2str(threshold) 'pct_AllSubs_Vertices.jpg'];


h = pie(restverts(networksexist));
delete(findobj(h,'Type','text'))

for k = 1:length(h)/2
    
    set(h(k*2-1), 'FaceColor', color_legend(k,:));
    
end

legend(network_legend,'Orientation','vertical','Location','eastoutside', 'FontSize',20)

print(gcf,[outdir filename],'-dpng','-r300');

close gcf



% Task variants vertices pie chart

networksexist = find(taskverts > 0);
network_legend = network_names(networksexist);
color_legend = rgb_colors(networksexist,:);

filename = ['TaskVariants_PieChart_' num2str(threshold) 'pct_AllSubs_Vertices.jpg'];


h = pie(taskverts(networksexist));
delete(findobj(h,'Type','text'))

for k = 1:length(h)/2
    
    set(h(k*2-1), 'FaceColor', color_legend(k,:));
    
end

legend(network_legend,'Orientation','vertical','Location','eastoutside', 'FontSize',20)

print(gcf,[outdir filename],'-dpng','-r300');

close gcf


% All variants pie chart

networksexist = find(allnets > 0);
network_legend = network_names(networksexist);
color_legend = rgb_colors(networksexist,:);

filename = ['AllVariants_PieChart_' num2str(threshold) 'pct_AllSubs_Variants.jpg'];


h = pie(allnets(networksexist));
delete(findobj(h,'Type','text'))

for k = 1:length(h)/2
    
    set(h(k*2-1), 'FaceColor', color_legend(k,:));
    
end

legend(network_legend,'Orientation','vertical','Location','eastoutside', 'FontSize',20)

print(gcf,[outdir filename],'-dpng','-r300');

close gcf


% Rest variants pie chart 

networksexist = find(restnets > 0);
network_legend = network_names(networksexist);
color_legend = rgb_colors(networksexist,:);

filename = ['RestVariants_PieChart_' num2str(threshold) 'pct_AllSubs_Variants.jpg'];


h = pie(restnets(networksexist));
delete(findobj(h,'Type','text'))

for k = 1:length(h)/2
    
    set(h(k*2-1), 'FaceColor', color_legend(k,:));
    
end

legend(network_legend,'Orientation','vertical','Location','eastoutside', 'FontSize',20)

print(gcf,[outdir filename],'-dpng','-r300');

close gcf


% Task variants pie chart

networksexist = find(tasknets > 0);
network_legend = network_names(networksexist);
color_legend = rgb_colors(networksexist,:);

filename = ['TaskVariants_PieChart_' num2str(threshold) 'pct_AllSubs_Variants.jpg'];


h = pie(tasknets(networksexist));
delete(findobj(h,'Type','text'))

for k = 1:length(h)/2
    
    set(h(k*2-1), 'FaceColor', color_legend(k,:));
    
end

legend(network_legend,'Orientation','vertical','Location','eastoutside', 'FontSize',20)

print(gcf,[outdir filename],'-dpng','-r300');

close gcf


%% Functions for calculating network assignment statistics

function input_vec = between_state_variant_assignment(Variant_Match_Vertices_Struct,vars)

%% Checks the match of a network variant assignment in one state to the assignment of the same vertices in the other state

input_vec = [];

if strcmp(Variant_Match_Vertices_Struct(vars).Network, 'Excluded')  %% if the variant was excluded for overlapping with the consensus network
    
    input_vec = [input_vec; 0 0];    %% do not include it in the dice correlation for network assignment
    
elseif Variant_Match_Vertices_Struct(vars).Network_Other_State ~= 99 && strcmp(Variant_Match_Vertices_Struct(vars).Network, Variant_Match_Vertices_Struct(vars).Network_Other_State)  %% if the variants are assigned to the same network
    
    input_vec = [input_vec; 1 1];  %% include it as a match for the dice correlation
    
else
    
    input_vec = [input_vec; 1 0];  %% else, set it as a non-match
    
end

end


function input_vec = across_subject_variant_assignment(Variant_Match_Vertices_Struct,vars)

%% Checks the match of a network variant assignment in one subject to the assignment of the same vertices in the other state of another subject


input_vec = [];

for perms = 1:numel(Variant_Match_Vertices_Struct(vars).Net_Assign_Perm)   %% for each variant in another subject
    
    if strcmp(Variant_Match_Vertices_Struct(vars).Network, 'Excluded')  %% if the variant in current subject was excluded for overlapping with the consensus network
        
        input_vec = [input_vec 0];  %% do not include it in the dice correlation for network assignment (both vectors will be 0s)
        
    elseif strcmp(Variant_Match_Vertices_Struct(vars).Network,Variant_Match_Vertices_Struct(vars).Net_Assign_Perm{perms})  %% if the variants are assigned to the same network in both subjects
        
        input_vec = [input_vec 1];  %% include it as a match for the dice correlation
        
    else
        
        input_vec = [input_vec 0];  %% else, set it as a non-match
        
    end

end

end


function input_vec = between_state_vertex_assignment(Variant_Match_Vertices_Struct,vars)

%% Checks the match of a network vertex in one state to the same vertex in the other state

input_vec = [];

for count = 1:length(Variant_Match_Vertices_Struct(vars).Vertex_Assignments)   %% for each vertex in each variant

    if Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count) == 99   %% if the vertex was excluded for overlapping with the consensus network
        
        input_vec = [input_vec; 0 0];  %% do not include it in the dice correlation for network assignment
        
    elseif Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count) == Variant_Match_Vertices_Struct(vars).Network_Other_State_Vertex(count)  %% if the vertices are assigned to the same network
        
        input_vec = [input_vec; 1 1];  %% include it as a match for the dice correlation
        
    else
        
        input_vec = [input_vec; 1 0];  %% else, set it as a non-match
        
    end

end

end


function input_vec = across_subject_vertex_assignment(Variant_Match_Vertices_Struct,vars)

%% Checks the match of a vertex network assignment in one subject to the same vertex in the other state of another subject

input_vec = [];

for nsubs = 1:length(Variant_Match_Vertices_Struct(vars).Net_Assign_Perm_Vertex)   %% for all other subjects
    
    tempdice = [];
    
    for count = 1:length(Variant_Match_Vertices_Struct(vars).Net_Assign_Perm_Vertex{nsubs})  %% for each vertex in another subject
        
        if Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count) == 99   %% if the vertex in current subject was excluded for overlapping with the consensus network
            
            tempdice = [tempdice; 0];  %% do not include it in the dice correlation for network assignment (both vectors will be 0s)
            
        elseif Variant_Match_Vertices_Struct(vars).Vertex_Assignments(count) == Variant_Match_Vertices_Struct(vars).Net_Assign_Perm_Vertex{nsubs}(count)  %% if the vertices are assigned to the same network in both subjects
            
            tempdice = [tempdice; 1];  %% include it as a match for the dice correlation
            
        else
            
            tempdice = [tempdice; 0];  %% else, set it as a non-match
            
        end
        
    end
    
    input_vec = [input_vec tempdice];  %% add values to the output vectors for correlation
    
end

end


function input_vec = between_state_vertex_location_assignment(Variant_Match_Vertices_Struct,vars,subs,subjects)

%% Checks if a vertex matches in both spatial location and network assignment across states

input_vec = [];

for verts = 1:length(Variant_Match_Vertices_Struct(vars).Vertex_Assignments)  %% for each vertex in a given variant
    
    if Variant_Match_Vertices_Struct(vars).Vertex_Assignments(verts) == 99  %% if the vertex was excluded for overlapping with the consensus network
        
        tmpdice1 = 0;  %% ignore this vertex in the correlation (set both values in each vector to 0)
        tmpdice2 = 0;
        
        input_vec = [input_vec; tmpdice1 tmpdice2];  %% add values to the output vectors for correlation
        
    else  %% if the vertex was not excluded
        
        tmpdice1 = 1;  %% set the value of the original vertex to 1 (as it is the "true" assignment)
        tmpdice2 = 0;  %% set the value of the vertex in the opposite state to 0 (non-matching) by default
        
        for varsloop = 1:size(Variant_Match_Vertices_Struct,2)  %% loop through all other variants to check if any other vertices are in the same location with the same network assignment
            
            if strcmp(Variant_Match_Vertices_Struct(varsloop).Subject, subs(subjects)) && ~strcmp(Variant_Match_Vertices_Struct(vars).State, Variant_Match_Vertices_Struct(varsloop).State) && tmpdice2 == 0 %% finds each variant for the current subject in the opposite state if no match has been found yet
                
                for vertsmatch = 1:length(Variant_Match_Vertices_Struct(varsloop).Variant_Vertices)  %% loops through each vertex in each variant
                    
                    if Variant_Match_Vertices_Struct(vars).Variant_Vertices(verts) == Variant_Match_Vertices_Struct(varsloop).Variant_Vertices(vertsmatch) && Variant_Match_Vertices_Struct(vars).Vertex_Assignments(verts) == Variant_Match_Vertices_Struct(varsloop).Vertex_Assignments(vertsmatch)  %% if any vertices have the same index on the CIFTI spatial map and are assigned to the same network
                        
                        tmpdice2 = 1;  %% include it as a match for the dice correlation
                        
                    end
                end
            end
        end
        
        input_vec = [input_vec; tmpdice1 tmpdice2];  %% add values to the output vectors for correlation
        
    end
end

end


function input_vec = across_subject_vertex_location_assignment(Variant_Match_Vertices_Struct,vars,subs,subjects)

%% Checks if a vertex matches in both spatial location and network assignment across states

input_vec = [];

substemp = subs;  %% creates a temporary variable with all subject names
substemp(subjects) = [];  %% leaves current subject out of variable

for verts = 1:length(Variant_Match_Vertices_Struct(vars).Vertex_Assignments)  %% for each vertex in the current variant
    
    dicetemp = [];  %% set matrix to contain comparisons across subjects
    
    if Variant_Match_Vertices_Struct(vars).Vertex_Assignments(verts) == 99  %% if the vertex in current subject was excluded for overlapping with the consensus network
        
        tmpdice = 0;  %% ignore this vertex in the correlation (set both values in each vector to 0)
        
        dicetemp = [dicetemp repmat(tmpdice,1,numel(substemp))];  %% do the same for the comparisons to all other subjects
        
    else  %% else, if vertex did not overlap with consensus
        
        for substmp = 1:numel(substemp)  %% for every other subject
            
            tmpdice = 0;  %% set the value of the vertex in the opposite state to 0 (non-matching) by default
            
            for varsloop = 1:size(Variant_Match_Vertices_Struct,2)  %% loop through all other variants to check if any other vertices are in the same location with the same network assignment
                
                if strcmp(Variant_Match_Vertices_Struct(varsloop).Subject, substemp(substmp)) && ~strcmp(Variant_Match_Vertices_Struct(vars).State, Variant_Match_Vertices_Struct(varsloop).State) && tmpdice == 0  %% finds each variant for the current (other) subject in the opposite state if no match has been found yet
                    
                    for vertsmatch = 1:length(Variant_Match_Vertices_Struct(varsloop).Variant_Vertices)  %% loops through each vertex in each variant
                        
                        if Variant_Match_Vertices_Struct(vars).Variant_Vertices(verts) == Variant_Match_Vertices_Struct(varsloop).Variant_Vertices(vertsmatch) && Variant_Match_Vertices_Struct(vars).Vertex_Assignments(verts) == Variant_Match_Vertices_Struct(varsloop).Vertex_Assignments(vertsmatch)  %% if any vertices in this subject in the opposite state have the same index on the CIFTI spatial map and are assigned to the same network
                            
                            tmpdice = 1;  %% include it as a match for the dice correlation
                            
                        end
                    end
                end
            end
            
            dicetemp = [dicetemp tmpdice];  %% add values to the output vectors for correlation (columns are subjects, rows are variants)
            
        end
    end
    
    input_vec = [input_vec; dicetemp];  %% add values to the output vectors for correlation
    
end

end

