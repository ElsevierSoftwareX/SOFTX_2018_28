function varargout = BedformsATM_App3_ScaleBasedDiscrimination(varargin)
% BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION MATLAB code for BedformsATM_App3_ScaleBasedDiscrimination.fig
%      BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION, by itself, creates a new BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION or raises the existing
%      singleton*.
%
%      H = BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION returns the handle to a new BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION or the handle to
%      the existing singleton*.
%
%      BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION.M with the given input arguments.
%
%      BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION('Property','Value',...) creates a new BEDFORMSATM_APP3_SCALEBASEDDISCRIMINATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BedformsATM_App3_ScaleBasedDiscrimination_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BedformsATM_App3_ScaleBasedDiscrimination_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BedformsATM_App3_ScaleBasedDiscrimination

% Last Modified by GUIDE v2.5 19-Nov-2016 13:14:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BedformsATM_App3_ScaleBasedDiscrimination_OpeningFcn, ...
                   'gui_OutputFcn',  @BedformsATM_App3_ScaleBasedDiscrimination_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BedformsATM_App3_ScaleBasedDiscrimination is made visible.
function BedformsATM_App3_ScaleBasedDiscrimination_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BedformsATM_App3_ScaleBasedDiscrimination (see VARARGIN)

axes(handles.axes1);
imshow('Scales.png');

handles.RunCorrelation = 0;
handles.RunStatistic   = 0;

handles.SaveFigure     = 0;
handles.figextension   = 'PDF';

% Choose default command line output for BedformsATM_App3_ScaleBasedDiscrimination
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = BedformsATM_App3_ScaleBasedDiscrimination_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in SelectProject.
function SelectProject_Callback(hObject, eventdata, handles)
% hObject    handle to SelectProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get projects path for uipickfiles
DocumentsPath = pwd;

% Get .MAT files
handles.projectfolder = fullfile(DocumentsPath, 'Projects');
handles.Files = uipickfiles('FilterSpec',handles.projectfolder);
[~,handles.j] = size(handles.Files);

if ( handles.j==0 )
    msgbox('No files had been loaded.');
    return;
end

load (handles.Files{1},'signalData','signifLevel','deltaFreq',...
                       'motherWavelet','printpath','sampleFreq');

handles.ableCorrelation = signalData(1,1)<0;

% Get the path to save the output
aux = find(printpath==filesep);
handles.printpath = [ printpath(1:aux(end-1)) 'BFPDiscrimination' filesep];

[~,~,~] = mkdir(handles.printpath);

% Read data to show in the main window
set (handles.FreqIncrText   ,'String',num2str(deltaFreq     ));
set (handles.SignifLevelText,'String',num2str(signifLevel   ));
set (handles.NumProfilesText,'String',num2str(handles.j));
set (handles.WavText        ,'String',motherWavelet);

% Save the largest scale of each profile to set the target
handles.largestscale(handles.j) = 0;

% Count the number of profiles with largest scale less than 100 (bars)
cont = 0;

for file = 1:handles.j 
    
    load (handles.Files{file},'globalWs','globalSignif',...
          'coneOfInfluence','fourierPeriod','j'); 
    
    maxPeriod = max(coneOfInfluence);
    
    fun = log10(globalWs);
    
    m = mean(fun);
    sd = std(fun);
    fun = (fun-m)/sd;
    
    param = 0.0001*(max(fun)-min(fun));
    [maxp,~] = peakdet( fun, param );
    
    [x,~] = size(maxp);
    
    target = fourierPeriod(maxp(1,1));
    
    for i=1:x
        % Just if the peak pass the confidence line and is inside the cone
        % of influence
        if ( globalWs( maxp(i,1) ) >= globalSignif( maxp(i,1) ) ) && fourierPeriod(maxp(i,1)) <= maxPeriod
            target = fourierPeriod(maxp(i,1));
            %nxtPeriod < 200 && 
            if i+1<=x && fourierPeriod(maxp(i+1,1)) <= maxPeriod
                target = fourierPeriod(maxp(i+1,1)); 
            end%
        end
        
        if target<100
            cont = cont + 1;
            target = maxPeriod;
        end
    end
    
    handles.largestscale(j) = target;
end

% Largest Scale
m = mean(handles.largestscale);
handles.BFPScale2 = m;
log_m = log10(m);

if cont > 0.95 * handles.j
    handles.scale3 = 'Dunes';
    set(handles.Scale3Text, 'string', 'Large Dunes');
    x = [190 + 60*log_m, 190 + 60*log_m];
    w = 190 + 60*log_m;
else
    handles.scale3 = 'Bars';
    set(handles.Scale3Text, 'string', 'Bars');
    x = [190 + 60*log_m, 190 + 60*log_m];
    w = 190 + 60*log_m;
end

y = [  326   140];
line(x,y,'LineWidth',2.0,'Color','r');

if sampleFreq<=0.12
    handles.scale1 = 'Ripples';
    set(handles.Scale1Text, 'string', 'Ripples');
    x = [  100   100];
    v = 100;
    y = [  326   140];
    line(x+0.5,y,'LineWidth',2.0,'Color','r');
    handles.BFPScale1 = 0.60;
else
    handles.scale1 = 'Dunes';
    set(handles.Scale1Text, 'string', 'Small Dunes');
    x = [  159   159];
    v = 159;
    y = [  326   140];
    handles.BFPScale1 = 5;
    line(x+0.5,y,'LineWidth',2.0,'Color','r');
end

y = [ 200 200 ];
x = [   w   v ];

line(x+0.5,y,'LineWidth',2.0,'Color','r','Marker','d');

set(handles.MinWL, 'string', num2str( handles.BFPScale1, '%.2f' ) );
set(handles.MaxWL, 'string', num2str( handles.BFPScale2, '%.2f' ) );

if strcmp(handles.scale3,'Dunes') && strcmp(handles.scale1,'Dunes')
    set(handles.Scale2Text, 'string', 'Middle Dunes');
elseif strcmp(handles.scale3,'Dunes') && strcmp(handles.scale1,'Ripples')
    set(handles.Scale2Text, 'string', 'Small and Middle Dunes');
elseif strcmp(handles.scale3,'Bars') && strcmp(handles.scale1,'Ripples')
    set(handles.Scale2Text, 'string', 'Dunes');
else 
    set(handles.Scale2Text, 'string', 'Middle and Large Dunes');
end

% Update handles structure
guidata(hObject, handles);

function MinWL_Callback(hObject, eventdata, handles)
% hObject    handle to MinWL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    value = str2double(get (hObject, 'string') );
    if ~isnan(value)
        handles.BFPScale1 = value;
    end
    set(hObject, 'string', num2str( handles.BFPScale1, '%.2f' ) );
guidata(hObject,handles);


function MaxWL_Callback(hObject, eventdata, handles)
% hObject    handle to MaxWL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    value = str2double(get (hObject, 'string') );
    if ~isnan(value)
        handles.BFPScale2 = value;
    end
    set(hObject, 'string', num2str( handles.BFPScale2, '%.2f' ) );
guidata(hObject,handles);

% --- Executes on button press in SaveFig.
function SaveFig_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.SaveFigure = get(hObject,'Value'); % Saves the figure in the fileFolder   

if handles.SaveFigure
    str = {'fig','jpg','tiff','PDF'};
    [choice,ok] = listdlg('PromptString','Select a format:',...
                                'SelectionMode','single',...
                                'ListString',str);
    if ~ok
        handles.figextension = '-';
        set(hObject,'Value',0);
        return;
    end
    handles.figextension = str{choice};
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in CorrAnalysis.
function CorrAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to CorrAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.ableCorrelation
    handles.RunCorrelation = get(hObject,'Value');
else
    msgbox('The correlation analysis is only available for water depth bedform fields.');
    set(hObject,'Value',0);
end

% Update handles structure
guidata(hObject, handles);


%% Not available for this version
% --- Executes on button press in StatAnalysis.
function StatAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to StatAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.RunStatistic = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in RunTheProgram.
function RunTheProgram_Callback(hObject, eventdata, handles)
% hObject    handle to RunTheProgram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%close();
%% Start discrimination
load (handles.Files{1}, 'motherWavelet', 'wltParameter', 'deltaFreq' ,...
                        'TitleFig'     , 'signalName'  , 'sampleName',...
                        'signalUnit'   , 'sampleUnit'  , 'signifLevel',...
                        'signalDataX');

bedformProfile(handles.j).name = 'fake value to preallocate the structure';
forAnalysis       = cell(1,handles.j);

% Structure to save info
% .name         <> Name of the project + ordinal number
% .path         <> Path to Name_BedformsDiscrimination.mat File
% .signalAbcise <> X coordinates
% .n            <> Profile points
% .n13          <> 1st level
% .n23          <> 2nd level
% .n33          <> 3rd level

% Arrays to plot each level

inodes = length(signalDataX);
jnodes = handles.j;

FileBed_dataMX   = zeros(inodes,jnodes); 
FileBed_dataMY   = zeros(inodes,jnodes); 
FileBed_dataMZ   = zeros(inodes,jnodes); 
FileBed_dataMZ13 = zeros(inodes,jnodes);
FileBed_dataMZ23 = zeros(inodes,jnodes);
FileBed_dataMZ33 = zeros(inodes,jnodes);


% Scales to perform the significance analysis
lowerScale = 0.001;
upperScale = 10000;

% ...\MyProject\BFPDiscrimination\
toprint = handles.printpath;
[~,~,~] = mkdir([handles.printpath 'Output']);

sParamPower = [-12:-1 0:2:60];

%% Getting the smoothed version of the signal

wb = waitbar(0,'1','Name','Progress...');

for file = 1:handles.j
    
    waitbar(file/handles.j,wb,sprintf(['Working on signal ' num2str(file) '...']));
    load (handles.Files{file},'signalDataX','signalData','sampleFreq','n_3D_yJ',...
                                'globalWs','fourierPeriod','globalSignif','j');

    %%
    % Profile info
    bedformProfile(j).name         = [TitleFig ' ' num2str(j)];
    bedformProfile(j).signalAbcise = signalDataX;
    bedformProfile(j).n            = signalData;
    
    % Smooth versions of the signal
    smoothData = zeros(length(signalData),length(sParamPower));
    
    m  = mean(signalData);
    sd =  std(signalData);
    NsignalData = (signalData-m)/sd;
    
    % Actual value 12.5% of max amplitude and using normalized data
    param = 0.125 * ( max( NsignalData )-min(NsignalData));
    [~, mintab] = peakdet(NsignalData, param);

    troughSignal = signalData(mintab(:,1))';
    
    for nSParam = 1:length(sParamPower)

        Z = SMOOTHN(signalData,2^(sParamPower(nSParam)))'; % Base of the S Parameter
        delta = 0;
        troughSmooth = [];
        m  = mean(Z);
        sd =  std(Z);
        NsignalData = (Z-m)/sd;
        param = 0.125 * ( max( NsignalData )-min(NsignalData));
        [~, mintab] = peakdet(NsignalData, param);
        
        if isempty(mintab)==0
            troughSmooth = Z(mintab(:,1));
        end;
        
        if length(troughSignal)==length(troughSmooth)
            delta = mean(abs(troughSmooth-troughSignal));
        end;
        Z = Z-delta; 
        smoothData(:,nSParam) = Z;

    end;

    %% Getting the 3rd level of the BFP 
    target = handles.largestscale(j);
    
    level2Index = 1;
    level3Index = 1;    
    signal3 = [];
    
    for rrg = 21:1:length(sParamPower)
        % Wavelet transform
            [wltPower, fourierPeriod, ~] =...
            waveletTransform (motherWavelet,wltParameter, smoothData(:,rrg)',...
                               sampleFreq   ,deltaFreq );

        % Wavelet significance
            [~, ~,globalWs,~,~,~] =...
            statisticsWlt(motherWavelet, wltParameter, smoothData(:,rrg)',sampleFreq,...
            deltaFreq, signifLevel, wltPower,lowerScale, upperScale);
        
        % Estimate the closer peak to the target
        ScaledSig = log10(globalWs);
        m    = mean(ScaledSig);
        sd   = std (ScaledSig);
        ScaledSig = (ScaledSig-m)/sd;
        
        param = 0.01*(max(ScaledSig)-min(ScaledSig));
        
        [maxtab,~] = peakdet( ScaledSig, param );

        if ( ~isempty(maxtab) ) 
            maxp = fourierPeriod(maxtab(:,1));
        end
        x = length(maxp);
        if x>0
            if x<=5
                signal3(level3Index,1) = rrg;
                [~,cv] = searchclosest(maxp,target);
                signal3(level3Index,2) = abs(target-cv);
                level3Index = level3Index+1;
            end
        end;
    end;    
    
    % Getting the signal at level 3
    position3 = max(find(signal3(:,2)==min(signal3(:,2))));
    level3Position = signal3(position3,1);
    % 3rd Level of th BFP
    bedformProfile(j).n33 = smoothData(:,level3Position);

    % *********************************************************************
    
    %% Getting the 2nd level of the BFP 
    forRipples = zeros(length(signalData),20);
    % Estimate the closer peak to the ripple threshold
    target = handles.BFPScale1;
    signal2 = [];
    
    for rrg = 1:1:20
        
        forRipples(:,rrg) = signalData'-smoothData(:,rrg);

        % Wavelet transform
            [wltPower, fourierPeriod, ~] =...
            waveletTransform (motherWavelet,wltParameter, forRipples(:,rrg)',...
                               sampleFreq   ,deltaFreq );
        % Wavelet significance
            [~, ~,globalWs,~,~,~] =...
            statisticsWlt(motherWavelet, wltParameter, forRipples(:,rrg)',sampleFreq,...
            deltaFreq, signifLevel, wltPower,lowerScale, upperScale);
  
        ScaledSig = log10(globalWs);
        m = mean(ScaledSig);
        sd = std(ScaledSig);
        
        ScaledSig = (ScaledSig-m)/sd;
        
        param = 0.05*(max(ScaledSig)-min(ScaledSig));
        
        [maxtab,~] = peakdet( ScaledSig, param );
        
        if ( ~isempty(maxtab) ) 
            maxp = fourierPeriod(maxtab(:,1));
        end
        
        if isempty(maxp)==0 
            signal2(level2Index,1) = rrg;
            signal2(level2Index,2) = abs( target-maxp(1) );
            level2Index = level2Index+1;
        end;

    end;
    
    % Getting the signal at level 2
    position2 = min( find(signal2(:,2)==min(signal2(:,2))) );
    level2Position = signal2(position2,1);

    bedformProfile(j).n23trended = smoothData(:,level2Position);

    bedformProfile(j).n23 = bedformProfile(j).n23trended-...
                                 bedformProfile(j).n33;
    
    bedformProfile(j).n13 = bedformProfile(j).n'-...
                                 bedformProfile(j).n23trended;

    % *********************************************************************
    
    %%
    bedformProfile(j).n33 = bedformProfile(j).n33+...
                                 mean(bedformProfile(j).n23)+...
                                 mean(bedformProfile(j).n13);

    bedformProfile(j).n23 = bedformProfile(j).n23-mean(bedformProfile(j).n23);
    bedformProfile(j).n13 = bedformProfile(j).n13-mean(bedformProfile(j).n13);
    
    % Threshold to minimize deep-through deformations
    
    [~,filteredSignalData] = modifiedZscoreForOutliers(bedformProfile(j).n13,1.25);
    bedformProfile(j).n13 = filteredSignalData;
                          
    FileBed_dataMX  (:,j)  = bedformProfile(j).signalAbcise;
    FileBed_dataMY  (:,j)  = n_3D_yJ(:,2);
    FileBed_dataMZ  (:,j)  = bedformProfile(j).n;
    FileBed_dataMZ13(:,j)  = bedformProfile(j).n13;
    FileBed_dataMZ23(:,j)  = bedformProfile(j).n23;
    FileBed_dataMZ33(:,j)  = bedformProfile(j).n33;
        
end;

delete(wb);

%% Smooth each Level
FileBed_dataMZ33 = SMOOTHN(FileBed_dataMZ33,2^7);
FileBed_dataMZ23 = SMOOTHN(FileBed_dataMZ23,2^5);

smooth_error = FileBed_dataMZ - FileBed_dataMZ33 - FileBed_dataMZ23 - FileBed_dataMZ13;

% FileBed_dataMZ33 = FileBed_dataMZ33;
FileBed_dataMZ23 = FileBed_dataMZ23 + smooth_error;

for j = 1:handles.j
    bedformProfile(j).n23 = FileBed_dataMZ23(:,j);
    bedformProfile(j).n33 = FileBed_dataMZ33(:,j);
    
    % *********************************************************************
    
    %% Creating the *.mat file for print and statistical analysis. 
       
    nameOfMATFile = [handles.printpath 'Output' filesep bedformProfile(j).name,'-BedformsDiscrimination'];
    BedformsDiscrimination = [bedformProfile(j).signalAbcise...
                              bedformProfile(j).n'...
                              bedformProfile(j).n33...
                              bedformProfile(j).n23...
                              bedformProfile(j).n13];
    bedformData = bedformProfile(j);
    
    save ([nameOfMATFile '.mat'], 'BedformsDiscrimination', 'bedformData',...
        'toprint','sampleName','sampleUnit', 'motherWavelet', 'wltParameter',...
        'deltaFreq', 'signifLevel' ); 
    forAnalysis{1,j} = [nameOfMATFile '.mat'];
                          
end


if (handles.SaveFigure)
    wb = waitbar(0,'1','Name','Progress...');
    for j = 1:handles.j

        waitbar(j/handles.j,wb,sprintf(['Printing signal ' num2str(j) '...']));

        %% Plotting the signal at different levels
        
        plotBFS ( bedformProfile(j).signalAbcise, bedformProfile(j).n,...
                      bedformProfile(j).n13, bedformProfile(j).n23,...
                      bedformProfile(j).n33, bedformProfile(j).name,...
                      handles.figextension, handles.printpath, sampleName, sampleUnit );
    end
    delete(wb);
end

% For Statistical Analysis
BFPScale1 = handles.BFPScale1;
BFPScale2 = handles.BFPScale2;

% forAnalysis stores the path to all BedformsDiscrimination .MAT Files
save ( [handles.printpath TitleFig 'Discrimination_ForAnalysis.mat'], 'forAnalysis',...
       'BFPScale1','BFPScale2','FileBed_dataMX','FileBed_dataMY','FileBed_dataMZ',...
       'FileBed_dataMZ13','FileBed_dataMZ23','FileBed_dataMZ33','motherWavelet',...
       'wltParameter', 'deltaFreq', 'signifLevel');
   
% Not available for this version 
% if ( handles.RunStatistic )
%     runStatisticAnalysis  ( forAnalysis, BFPScale1, BFPScale2 );
% end

if ( handles.RunCorrelation )
    runCorrelationAnalysis( forAnalysis );
end

%% Surf each level
if (handles.j>10)
    SurfBedformLevel( FileBed_dataMX, FileBed_dataMY, FileBed_dataMZ, FileBed_dataMZ13,1,handles.printpath );
    SurfBedformLevel( FileBed_dataMX, FileBed_dataMY, FileBed_dataMZ, FileBed_dataMZ23,2,handles.printpath );
    SurfBedformLevel( FileBed_dataMX, FileBed_dataMY, FileBed_dataMZ, FileBed_dataMZ33,3,handles.printpath );
end
display('Successful run!!!');
