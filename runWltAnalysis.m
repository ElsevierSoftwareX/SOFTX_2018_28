function runWltAnalysis(isTest, numberOfFiles, projectfolder, SaveFigure, ...
         figextension, signalName, signalUnit, sampleName, sampleUnit, ...
         TitleFig, n_2D, n_3D, motherWavelet, wltParameter, ...
         deltaFreq, signifLevel, h, plot_a, plot_b, plot_c)
%% .
    % This function execute the wavelet analysis for a profile.
    
    % This function performs the wavelet analysis, print the resulting
    % figure depending of the user and save the .MAT files that stores the
    % results.

%%
% Settings to save figures

outputPath = fullfile(projectfolder, 'Wavelet_Output');
[~,~,~] = mkdir(outputPath);
printpath = fullfile(outputPath, 'Figures');


% Create a figure to save outputs
if SaveFigure
    [~,~,~] = mkdir(printpath);

    % Invisible fig for speed up the runtime
    h = figure('Visible','off');
    
    plot_a = subplot('position',[0.13 0.78 0.83 0.15]);
    plot_b = subplot('position',[0.13 0.42 0.83 0.27]);
    plot_c = subplot('position',[0.13 0.10 0.83 0.20]);

end

%% 
% Run program

allPath = fullfile(outputPath, 'Data_Output');
[~,~,~] = mkdir (allPath);

if ~isTest
    wb = waitbar(0,'1','Name','Progress...');
end

for j = 1: numberOfFiles
    
    if ~isTest
        waitbar(j/numberOfFiles,wb,sprintf(['Working on signal ' num2str(j) '...']));
    end
    
    n_2D_yJ = n_2D{j}; % Data in 2D
    n_3D_yJ = n_3D{j};  % Data in 3D
    
    signalDataX = n_2D_yJ (:,2);  
    signalData  = transpose (n_2D_yJ (:,3));
    
    sampleFreq  = signalDataX(2)-signalDataX(1);
    
    titleFigA   = [TitleFig ' ' num2str(j)];
    filename    = [strcat(strtrim(titleFigA), '_',  motherWavelet,  num2str(wltParameter))];
    %% Wavelet
    [wltPower, fourierPeriod, coneOfInfluence] =...
     waveletTransform (motherWavelet,wltParameter, signalData,...
                       sampleFreq   ,deltaFreq );
    %% Statistical analysis
    lowerScale = 0.001;
    upperScale = 10000;

    [wltSignif, globalSignif,globalWs,lag1,scaleAveSignif,scaleAvg] =...
    statisticsWlt(motherWavelet, wltParameter, signalData,sampleFreq,...
    deltaFreq, signifLevel, wltPower,lowerScale, upperScale);

    %% Plot figures
    [peaksCross, xPeakspos,yPeakspos] = getWltpeaks(globalWs,globalSignif,fourierPeriod);
    
    if SaveFigure
        figure(h);
        set(h,'Visible','off');
        plotFiguresWlt ( motherWavelet, wltParameter,...
        titleFigA , signalDataX, signalData, wltPower, fourierPeriod,...
        coneOfInfluence, wltSignif, globalSignif, globalWs, lag1,...
        SaveFigure, figextension, printpath, plot_a, plot_b, plot_c,...
        signalName, signalUnit, sampleName, sampleUnit,...
        peaksCross, xPeakspos, yPeakspos);
    
    end
    
    %% Save all variables
    save ([allPath filesep filename '.mat'],'j', 'signalDataX', 'signalData', 'sampleFreq',...
    'motherWavelet', 'wltParameter', 'deltaFreq','signifLevel', 'n_3D_yJ',...
    ... % from waveletTransform
    'wltPower', 'fourierPeriod','coneOfInfluence',...
    ... % from statisticsWlt
    'wltSignif', 'globalSignif', 'globalWs','scaleAveSignif','scaleAvg',...
    'lag1',...
    ... % To print Data
    'TitleFig', 'printpath', 'peaksCross', 'xPeakspos', 'yPeakspos','signalName',...
    'signalUnit', 'sampleName', 'sampleUnit',...
    '-mat');

end;

if ~isTest
    delete(wb);
end

if SaveFigure 
    close(h);
end