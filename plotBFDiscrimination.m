function plotBFDiscrimination (figextension)

    % This function prints each profile and its corresponding levels
    % This is designed to be called from the MATLAB console after run the
    % BedformsATM_ScaleBasedDiscrimination
    % User should chose the .MAT files to print

    % INPUTS
    %     1. figextension: Format to save the figure.
    % 
    %%
    switch figextension
        case 'jpg'
        case 'fig'
        case 'PDF'
        case 'tiff'
        otherwise
            display( 'The figure format inserted is not suported.' )
            display( 'You can print in ''jpg'', ''tiff'', ''PDF'' or ''fig'' format.');
            return;
    end
    
    DocumentsPath = pwd;
    projectfolder = fullfile(DocumentsPath, 'Projects');
    
    files = uipickfiles('FilterSpec',projectfolder);
    [~,nfiles] = size(files);
    for fil = 1:nfiles
        load(files{fil}, 'bedformData','toprint','sampleName','sampleUnit' ); 
        plotBFS (     bedformData.signalAbcise, bedformData.n,...
                      bedformData.n13, bedformData.n23,...
                      bedformData.n33, bedformData.name,...
                      figextension, toprint, sampleName, sampleUnit );
    end
