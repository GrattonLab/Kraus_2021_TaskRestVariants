function verts=generate_rotated_parcellation(giftidat,iterations,rotations,hem,outfile,neighborfilepath,atlasdir,surfdir)

%% Rotates parcels (variants) within each hemisphere for each subject
%
% This function randomly rotates parcels (vairants) in a GIFTI (hemisphere)
% the amount of times that is specified by iterations. The rotation in each
% plane is provided to the function for each iteration as well. For each
% iteration, each parcel is first rotated according to the rotations 
% provided, then if necessary randomly rotated until it does not spatially
% overlap with the medial wall or with another parcel. A GIFTI containing
% all of these rotations is then saved out to the specified filepath and
% the vertices matrix (vertices x rotations) is also output into the
% workspace.
%
% INPUTS:
%
% -giftidat: a gifti data vector (numvertices x 1) that contains unique IDs
% for each parcel to be rotated (each vertex in a contiguous parcel should
% have the same value in the vector)
% -iterations: the number of rotations to perform within each hemisphere
% for each subject and state
% -rotations: a structure containing rotation values for the x, y, and z
% planes, there should be at least as many rotations as iterations
% requested
% -hem: indicates which hemisphere GIFTI is being rotated ('L' or 'R')
% -outfile: the output filepath and filename for the GIFTI file containing
% the rotations
% -neighborfilepath: path to node neighbor file generated from
% caret '-surface-topology-neighbors' (closest "neighbors" to each vertex)
% -atlasdir: the path to the directory that contains the 32k surface medial
% wall GIFTI files
% -surfdir: the path to the directory that contains the 32k surface medial
% wall GIFTI files
%
% OUTPUTS:
%
% -GIFTI rotation files: GIFTI spatial maps (numvertices x iterations) for
% which each timepoint in the file represents a spatial rotation. The GIFTI
% rotation files are used to make the CIFTI rotation files
% -GIFTI rotation matrix: matrix representing GIFTI spatial maps 
% (numvertices x iterations) is output into the workspace
%
% Originally from BAS at WashU. Edited by BK (01-2021)
%


%% Initialize Variables

% Set paths and constants
bufsize=16384;
% Read in node neighbor file generated from caret
% -surface-topology-neighbors (closest "neighbors" to each vertex)
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread(neighborfilepath,'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;
nonanneighbors = neighbors; nonanneighbors(isnan(neighbors)) = 380;
% Load normal wall mask
maskname = [atlasdir 'medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
% Load inflated sphere surface
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
spherecoords = sphere.vertices;

%% Randomly rotate parcels and ensure that none end up on the medial wall or overlap with other variants

parcelIDs = unique(giftidat(giftidat>0));  %% unique identifiers for each parcel
numparcels = length(parcelIDs);  %% number of parcels

disp(['Running ' num2str(iterations) ' random rotations of ' num2str(numparcels) ' parcels']);

rotmap = zeros(32492,iterations);  %% empty GIFTI matrix for each iteration of rotations
test = zeros(32492,3);

for iternum = 1:iterations
    
    xrot = rotations.rotations.xrot(iternum);  %% get the rotation value for each iteration
    yrot = rotations.rotations.yrot(iternum);
    zrot = rotations.rotations.zrot(iternum);
    
    rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];  %% translate the rotation value into coordinates
    rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
    rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
    
    parcelinds = [];
    
    for parcelnum = 1:numparcels
        
        parcelID = parcelIDs(parcelnum);
        
        realparcelindices = find(giftidat==parcelID);

        string{parcelnum} = ['Rotation ' num2str(iternum) ': parcel ' num2str(parcelnum)];
        if parcelnum==1; fprintf('%s',string{parcelnum}); else fprintf([repmat('\b',1,length(string{parcelnum-1})) '%s'],string{parcelnum}); end
        
        while true      %% Loop to randomly rotate parcels until they end up not on the medial wall and not overlapping with other variants

            indexcoords = spherecoords(realparcelindices,:)';  %% create the rotation by multiplying the rotation coordinates by the coords of the original parcel indices
         	xrotcoords = rotmat_x * indexcoords;
           	xyrotcoords = rotmat_y * xrotcoords;
           	xyzrotcoords = rotmat_z * xyrotcoords;

           	clear rotparcelindices
            for n = 1:length(realparcelindices)  %% generate rotated parcel indices for each real parcel index
                test(:,1) = xyzrotcoords(1,n); test(:,2) = xyzrotcoords(2,n); test(:,3) = xyzrotcoords(3,n);  %% get rotated parcel coordinates
                diff_test = sum(abs(spherecoords - test),2);  %% find the closest match to the differences between the 2 sets of coordinates
                [~, rotparcelindices(n)] = min(diff_test);  %% add index to the rotated parcel vector
            end
            
            % Stop loop if parcel is not on medial wall and does not
            % overlap with existing parcels, else generate random numbers
            % and re-rotate parcel
            
            if (numel(intersect(rotparcelindices,find(mask == 1))) == 0) && (isempty(parcelinds) || (numel(intersect(rotparcelindices,parcelinds)) == 0))
                
                parcelinds = [parcelinds rotparcelindices];  %% add indices to output if they do not overlap with medial wall or another rotation
                
                break
                
            else  %% else if an overlap occurs

                rng('shuffle')  %% shuffle the random number generator
                
                xrot = min(rotations.rotations.xrot) + (max(rotations.rotations.xrot)-min(rotations.rotations.xrot))*rand(1,1);  %% recalculate the rotation for this variant
                yrot = min(rotations.rotations.yrot) + (max(rotations.rotations.yrot)-min(rotations.rotations.yrot))*rand(1,1);
                zrot = min(rotations.rotations.zrot) + (max(rotations.rotations.zrot)-min(rotations.rotations.zrot))*rand(1,1);
                
                rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];  %% translate the rotation value into coordinates
                rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
                rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
            end
        end
        
        % add vertices to each rotated parcel until it becomes the same
        % size as the original parcel
        
        doneadding = 0;
        while doneadding ==0
            rotatedparcel = zeros(size(mask));
            rotatedparcel(rotparcelindices) = 1;
            rotneighvals = rotatedparcel(neighbors(13:end,2:7));
            rotneighvals_top = rotatedparcel(neighbors(1:12,2:6));
            rotverts_toadd = [((sum(rotneighvals_top,2) > 3) .* (rotatedparcel(1:12)==0)) ; ((sum(rotneighvals,2) > 4) .* (rotatedparcel(13:end)==0))];
            if nnz(rotverts_toadd) == 0
                doneadding = 1;
            else
                rotatedparcel = rotatedparcel + rotverts_toadd;
                rotparcelindices = find(rotatedparcel);
            end
        end
        
        % check that each rotated parcel is the same size as the original
        % parcel, add or subtract vertices until this is the case
        
        samesize = 0;
        while samesize==0
            
            deltaverts = length(rotparcelindices) - length(realparcelindices);
            rotatedparcel = zeros(size(rotatedparcel)); rotatedparcel(rotparcelindices) = 1; rotatedparcel(logical(mask)) = 2;
            if sign(deltaverts) == 1
                
                borderverts = find((rotatedparcel==1) .* any(rotatedparcel(nonanneighbors(:,2:end))==0,2));
                
                if length(borderverts) >= deltaverts
                    rotparcelindices = setdiff(rotparcelindices,borderverts(1:deltaverts));
                    samesize = 1;
                else
                    rotparcelindices = setdiff(rotparcelindices,borderverts);
                end
            elseif sign(deltaverts) == -1
                
                borderverts = find((rotatedparcel==0) .* any(rotatedparcel(nonanneighbors(:,2:end))==1,2));
                
                if length(borderverts) >= abs(deltaverts)
                    rotparcelindices = union(rotparcelindices,borderverts(1:abs(deltaverts)));
                    samesize = 1;
                else
                    rotparcelindices = union(rotparcelindices,borderverts);
                end
            else
                samesize = 1;
            end
        end
        
        rotmap(rotparcelindices,iternum) = parcelID;  %% write out final variant onto the spatial rotation map
        
    end
    fprintf(repmat('\b',1,length(string{parcelnum})));
    
end
    
disp('finished all rotations');


save(gifti(single(rotmap)),outfile);  %% save GIFTI rotation map

verts = rotmap;  %% output GIFTI from function


