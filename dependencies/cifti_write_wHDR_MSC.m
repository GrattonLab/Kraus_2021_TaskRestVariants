function cifti_write_wHDR_MSC(data,ciftitemplatefile,outnamestem,workbenchdir,giftitemplatepath,varargin)  

%% Write out input GIFTI data from separate hemispheres into a CIFTI file containing both hemispheres
%
% This function takes an input concatenated GIFTI dataset 
% (numvertices x timepoints), changes the header, and converts it into 
% CIFTI format.
%
% INPUTS:
%
% -data: two concatenated GIFTIs (one from each hemisphere) after
% accounting for the medial wall
% -ciftitemplatefile/giftitemplatepath: the template file to use for
% writing each GIFTI to a CIFTI (can either use a CIFTI or GIFTI as the
% template)
% -outnamestem: the output filepath and filename for the CIFTI file
% -giftitemplatepath: the filepath and filename to the GIFTI used as a
% template for converting GIFTIs to CIFTIs
% -workbenchdir: the filepath to connectome workbench
% -varargin: the extension for the output file (dtseries, dconn, or pconn)
%
% OUTPUTS:
%
% -Writes a CIFTI file (numvertices x timepoints) at the specified filepath
%
% Originally from WashU, edited by BK (01-2021)
%

%% Write out CIFTI file

% set output filetype and extension (assumes dtseries if no input)

if ~isempty(varargin)
    filetype = varargin{1};
else
    filetype = 'dtseries';
end
extension = ['.' filetype '.nii'];

if length(outnamestem) > length(extension) && strcmp(outnamestem(end-(length(extension)-1):end),extension)
    outnamestem = outnamestem(1:end-length(extension));
end

% either create a GIFTI template from a CIFTI, or use an existing GIFTI
% template

deletestuff= 0;
if isempty(ciftitemplatefile)
    ciftitemplatefile = giftitemplatepath;
elseif strcmp(ciftitemplatefile(end-(length(extension)-1):end),extension)
    system([workbenchdir 'wb_command -cifti-convert -to-gifti-ext ' ciftitemplatefile ' WritingTemp.func.gii']);
    ciftitemplatefile = 'WritingTemp.func.gii';
    deletestuff = 1;
end

% read template file for header information
bufsize = 524288;
ciftiheadertext = textread(ciftitemplatefile,'%s','delimiter','\r','bufsize',bufsize);
if deletestuff
    delete('WritingTemp.func*')
end

% write header for file and write out final file in correct format (dconn,
% pconn, or dtseries)

switch filetype
    
    case 'dconn'
        
        for row = 1:length(ciftiheadertext)
            thisline = ciftiheadertext{row};
            if length(thisline) >= 4 && strcmp(thisline(1:4),'Dim1')
                thisline = ['Dim1="' num2str(size(data,2)) '"'];
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline) >= 16 && strcmp(thisline(1:16),'ExternalFileName')
                thisline = ['ExternalFileName="' outnamestem '.func.dat"'];
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif row>1 && strcmp(ciftiheadertext{row-1},'<Name>Provenance</Name>')
                thisline = '<Value>created with Matlab cifti_write function</Value>';
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline)>=70 && strcmp(thisline(1:70),'<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_BRAIN_MODELS"')
                thisline = '<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_BRAIN_MODELS" AppliesToMatrixDimension="0,1">';
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline)>=69 && strcmp(thisline(1:69),'<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_TIME_POINTS"')
            else
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            end           
            
        end
        system([workbenchdir 'wb_command -cifti-convert -from-gifti-ext ' outnamestem '.func.gii ' outnamestem '.dconn.nii']);
        
     case 'pconn'
        
        for row = 1:length(ciftiheadertext)
            thisline = ciftiheadertext{row};
            if length(thisline) >= 4 && strcmp(thisline(1:4),'Dim1')
                thisline = ['Dim1="' num2str(size(data,2)) '"'];
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline) >= 16 && strcmp(thisline(1:16),'ExternalFileName')
                thisline = ['ExternalFileName="' outnamestem '.func.dat"'];
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif row>1 && strcmp(ciftiheadertext{row-1},'<Name>Provenance</Name>')
                thisline = '<Value>created with Matlab cifti_write function</Value>';
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline)>=70 && strcmp(thisline(1:70),'<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_BRAIN_MODELS"')
                thisline = '<MatrixIndicesMap AppliesToMatrixDimension="0,1" IndicesMapToDataType="CIFTI_INDEX_TYPE_PARCELS">';
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            elseif length(thisline)>=69 && strcmp(thisline(1:69),'<MatrixIndicesMap IndicesMapToDataType="CIFTI_INDEX_TYPE_TIME_POINTS"')
            else
                dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            end           
            
        end
        system([workbenchdir 'wb_command -cifti-convert -from-gifti-ext ' outnamestem '.func.gii ' outnamestem '.pconn.nii']);

    case 'dtseries'
        for row = 1:length(ciftiheadertext)
            thisline = ciftiheadertext{row};
            if length(thisline) >= 4 && strcmp(thisline(1:4),'Dim1')
                thisline = ['Dim1="' num2str(size(data,2)) '"'];
            elseif length(thisline) >= 16 && strcmp(thisline(1:16),'ExternalFileName')
                thisline = ['ExternalFileName="' outnamestem '.func.dat"'];
            elseif row>1 && strcmp(ciftiheadertext{row-1},'<Name>Provenance</Name>')
                thisline = '<Value>created with Matlab cifti_write function</Value>';
            end
            
            dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
            
        end
        system([workbenchdir 'wb_command -cifti-convert -from-gifti-ext ' outnamestem '.func.gii ' outnamestem '.dtseries.nii']);

end

delete([outnamestem '.func.*'])

