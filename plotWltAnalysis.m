function plotWltAnalysis (figextension)
    % This function prints each profile and its corresponding wavelet
    % spectrum and global significance
    % This is designed to be called from the MATLAB console after run the
    % BedformsATM_WaveletAnalysis
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
      
        load(files{fil}, 'j', 'motherWavelet', 'wltParameter',...
        'TitleFig' , 'signalDataX', 'signalData', 'wltPower', 'fourierPeriod',...
        'coneOfInfluence', 'wltSignif', 'globalSignif', 'globalWs', 'lag1',...
        'printpath', 'peaksCross', 'xPeakspos', 'yPeakspos','signalName',...
        'signalUnit', 'sampleName', 'sampleUnit');
        
        titleFigA  = [TitleFig ' ' num2str(j)];
        
        if(fil == 1)
           
            titleFont = 16;
            labelFont = 16;
            axisFont = 14;
            [~,~,~] = mkdir(printpath);
            if strcmp(figextension,'fig')
                h = figure;
            else
                h = figure('Visible','off');
            end
            
            plot_a = subplot('position',[0.13 0.78 0.83 0.15]);
            plot_b = subplot('position',[0.13 0.42 0.83 0.27]);
            plot_c = subplot('position',[0.13 0.10 0.83 0.20]);
            
        end
        
        plotFiguresWlt ( motherWavelet, wltParameter,...
        titleFigA , signalDataX, signalData, wltPower, fourierPeriod,...
        coneOfInfluence, wltSignif, globalSignif, globalWs, lag1,...
        1, figextension, printpath, plot_a, plot_b, plot_c,...
        signalName, signalUnit, sampleName, sampleUnit,...
        peaksCross, xPeakspos, yPeakspos);
 
        
    end

    close(h);