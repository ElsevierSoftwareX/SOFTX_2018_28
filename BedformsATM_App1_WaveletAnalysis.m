function varargout = BedformsATM_App1_WaveletAnalysis(varargin)
% BEDFORMSATM_APP1_WAVELETANALYSIS MATLAB code for BedformsATM_App1_WaveletAnalysis.fig
%      BEDFORMSATM_APP1_WAVELETANALYSIS, by itself, creates a new BEDFORMSATM_APP1_WAVELETANALYSIS or raises the existing
%      singleton*.
%
%      H = BEDFORMSATM_APP1_WAVELETANALYSIS returns the handle to a new BEDFORMSATM_APP1_WAVELETANALYSIS or the handle to
%      the existing singleton*.
%
%      BEDFORMSATM_APP1_WAVELETANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEDFORMSATM_APP1_WAVELETANALYSIS.M with the given input arguments.
%
%      BEDFORMSATM_APP1_WAVELETANALYSIS('Property','Value',...) creates a new BEDFORMSATM_APP1_WAVELETANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BedformsATM_App1_WaveletAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BedformsATM_App1_WaveletAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BedformsATM_App1_WaveletAnalysis

% Last Modified by GUIDE v2.5 19-Nov-2016 13:05:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BedformsATM_App1_WaveletAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @BedformsATM_App1_WaveletAnalysis_OutputFcn, ...
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

% --- Executes just before BedformsATM_App1_WaveletAnalysis is made visible.
function BedformsATM_App1_WaveletAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BedformsATM_App1_WaveletAnalysis (see VARARGIN)

% Get MyDocuments path
handles.DocumentsPath = pwd;

% Set default wavelet parameters
handles.signifLevel = 0.95;
handles.deltaFreq   = 0.05;

handles.motherWavelet = 'MORLET';
handles.wltParameter  = 6;

% Set default output parameters
handles.TitleFigA    = '-';
handles.SignalName   = '-';
handles.SampleName   = '-';
handles.SignalUnit   = '-';
handles.SampleUnit   = '-';
handles.SaveFigure   =  0 ;
handles.figextension = 'x';

% Create arrays with tags to enable/disable panels

handles.DataFilesPanel = [ handles.SelectData ...
                          ,handles.ImportData ...
                          ,handles.Profiles];

handles.OutputPanel    = [ handles.FigATitle ...
                          ,handles.SignalName_popup ...
                          ,handles.SampleName_popup ...
                          ,handles.SignalUnit_popup ...
                          ,handles.SampleUnit_popup];
                       
set ( handles.DataFilesPanel, 'Enable', 'off' );
set ( handles.OutputPanel   , 'Enable', 'off' );
set ( handles.SaveFig       , 'Enable', 'off' );
set ( handles.RunTheProgram , 'Enable', 'off' );


set (handles.step1, 'FontWeight', 'bold','FontSize',10,'ForegroundColor', 'Blue');
set (handles.TabWavelet, 'BackgroundColor', [240/255, 240/255, 230/255],'ForegroundColor','Black');
set (handles.TabInfo   , 'BackgroundColor', [ 50/255, 200/255, 0      ],'ForegroundColor','White');

axes( handles.LogoPUCP );
imshow( 'Logo_PUCP.png' );

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = BedformsATM_App1_WaveletAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Tabs
% ------------------------------------------------------------------------

% --- Executes on button press in TabInfo.
function TabInfo_Callback(hObject, eventdata, handles)
% hObject    handle to TabInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set (handles.InfoPanel, 'Visible', 'on');
    
set (handles.TabWavelet, 'BackgroundColor', [240/255, 240/255, 230/255],'ForegroundColor','Black');
set (handles.TabInfo   , 'BackgroundColor', [ 50/255, 200/255, 0      ],'ForegroundColor','White');

set (handles.step1, 'Visible', 'off');
set (handles.step2, 'Visible', 'off');
set (handles.step3, 'Visible', 'off');
set (handles.step4, 'Visible', 'off');
set (handles.step5, 'Visible', 'off');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in TabWavelet.
function TabWavelet_Callback(hObject, eventdata, handles)
% hObject    handle to TabWavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set (handles.InfoPanel, 'Visible', 'off');

set (handles.TabInfo   , 'BackgroundColor', [240/255, 240/255, 230/255],'ForegroundColor','Black');
set (handles.TabWavelet, 'BackgroundColor', [ 50/255, 200/255, 0      ],'ForegroundColor','White');

set (handles.step1, 'Visible', 'on');
set (handles.step2, 'Visible', 'on');
set (handles.step3, 'Visible', 'on');
set (handles.step4, 'Visible', 'on');
set (handles.step5, 'Visible', 'on');

% Update handles structure
guidata(hObject, handles);

% Load Data
% ------------------------------------------------------------------------

% --- Executes on button press in CreateProject.
function CreateProject_Callback(hObject, eventdata, handles)
% hObject    handle to CreateProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get project's name from user input

strText = sprintf(['Enter a name for the project.\n' ...
                   'Use only alphanumeric characters']);
strTitle = 'Create a Project';

% The function insert_name returns an alphanumerc string from user input
handles.NameProject = insert_name({strText}, strTitle,{''});

if (isempty(handles.NameProject)) return; end

% Create a new folder, set alert if it already exists
handles.ProjectFolder = fullfile(handles.DocumentsPath, 'Projects', handles.NameProject);

% Check if the folder alredy exist
while ( exist(handles.ProjectFolder,'dir') )
        
    % Decide whether to clear the existing project or enter another name
    choice = questdlg(sprintf(['This project alredy exists.\n' ...
                               'If you continue all files will be cleared.']), ...
                               'The project name is already used', ...
                               'Continue','Enter other name','Continue');
    switch choice
        case 'Continue'
            % Clear existing project
            rmdir(handles.ProjectFolder, 's');

        case 'Enter other name'
            % Get another name
            handles.NameProject = insert_name({strText}, strTitle,{''});
            
            if (isempty(handles.NameProject)) return; end

            handles.ProjectFolder = fullfile(handles.DocumentsPath, 'Projects',...
                                     handles.NameProject);

     end

end

mkdir(handles.ProjectFolder);

% Enable Data Files 
set ( handles.SelectData, 'Enable', 'on' );

set (handles.step1, 'FontWeight', 'normal', 'FontSize', 8.0,'ForegroundColor','Black');
set (handles.step1, 'String', [ get(handles.step1,'String') ' (DONE)' ]);
set (handles.step2, 'FontWeight', 'bold', 'FontSize', 10.0, 'ForegroundColor', 'Blue');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SelectData.
function SelectData_Callback(hObject, eventdata, handles)
% hObject    handle to SelectData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call BedformsATM_WaveletAnalysis_ImportFiles
BedformsATM_App1_WaveletAnalysis_ImportFiles( {handles.ProjectFolder} );

set(handles.Profiles, 'String', 'Press the Preview Button');

set ( handles.DataFilesPanel, 'Enable', 'on' );

set (handles.step2, 'FontWeight', 'normal', 'FontSize', 8.0,'ForegroundColor','Black');
set (handles.step2, 'String', [ get(handles.step2,'String') ' (DONE)' ]);
set (handles.step3, 'FontWeight', 'bold', 'FontSize', 10.0, 'ForegroundColor', 'Blue');

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in ImportData.
function ImportData_Callback(hObject, eventdata, handles)
% hObject    handle to ImportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load signal data and coordinates data
load( fullfile(handles.ProjectFolder ,'input','inputfile.mat'), 'outFile','Dformat' );
load( fullfile(handles.ProjectFolder ,'input','XYZdata.mat'), 'XYZ' )

% handles.n          stores profiles in 2D
% handles.XYZcoor    stores XYZ coordinates 3D
handles.n = outFile;
handles.XYZcoor    = XYZ;
clear outFile;
clear XYZ;
handles.displayXY = Dformat;

[~,handles.numberOfFiles] = size(handles.n);

% prev stores a cellarray to show profiles in the pop-up menu
prev = cell(handles.numberOfFiles+1,1);
prev{1} = 'Select preview profile';

for nProfile = 1 : handles.numberOfFiles
    prev{nProfile+1} = [ 'Profile - ' num2str(nProfile) ]; 
end

set ( handles.Profiles, 'String', prev );

set ( handles.SaveFig      ,'Enable', 'on');
set ( handles.OutputPanel  ,'Enable', 'on');
set ( handles.RunTheProgram,'Enable', 'on');

set (handles.step3, 'FontWeight', 'normal', 'FontSize', 8.0,'ForegroundColor','Black');
set (handles.step3, 'String', [ get(handles.step3,'String') ' (DONE)' ]);
set (handles.step4, 'FontWeight', 'bold', 'FontSize', 10.0, 'ForegroundColor', 'Blue');
set (handles.step5, 'FontWeight', 'bold', 'FontSize', 10.0, 'ForegroundColor', 'Blue');

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in Profiles.
function Profiles_Callback(hObject, eventdata, handles)
% hObject    handle to Profiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get number of profile to plot
profile = get(hObject,'Value')-1;

% Text size for the plot
titleFont = 16;
labelFont = 16;
axisFont = 10;

% If plant view will be plot, load data to do it
% Window's size depends of handles.displayXY
if handles.displayXY
    coordinateX = handles.XYZcoor {profile}(:,1);
    coordinateY = handles.XYZcoor {profile}(:,2);
    figure('WindowStyle', 'Modal', 'Position', [350 100 800 600 ]);
else
    figure('WindowStyle', 'Modal', 'Position', [350 300 800 300 ]);
end

signalDataX = handles.n{profile}(:,2);
signalData  = handles.n{profile}(:,3);

% Range of the plots
plotXRange = [min(signalDataX), max(signalDataX)];

if  handles.displayXY 
    % Create an axes
    subplot('position',[0.10 0.60 0.83 0.30]);
end
plot(signalDataX, signalData, 'b','linewidth',0.8);
ylabel(strcat('Depth', {' '},'(m)'),'fontsize',labelFont);
xlabel(strcat('X', {' '},'(m)'),'fontsize',labelFont);
title ([ 'Profile - ' num2str(profile) ] ,'fontsize',titleFont,'FontWeight','bold');
set(gca, 'XLim', plotXRange(:));
set(gca,'fontsize',axisFont);
grid on;

if  handles.displayXY
    % Create another axes
    subplot('position',[0.10 0.10 0.83 0.30]);
    plot(coordinateX, coordinateY, 'black','linewidth',2.0);
    ylabel(strcat('Y', {' '},'(m)'),'fontsize',labelFont);
    xlabel(strcat('X', {' '},'(m)'),'fontsize',labelFont);
    title ([ 'Plant View - ' num2str(profile) ] ,'fontsize',titleFont,'FontWeight','bold');
    set(gca, 'fontsize', axisFont);
    grid on;
end

% Reset pop-up menu
set(hObject,'Value',1);

% Define parameters
% ------------------------------------------------------------------------
function SignifLevel_Callback(hObject, eventdata, handles)
% hObject    handle to SignifLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'Value');
switch contents
    case 1
        handles.signifLevel = 0.95;
    case 2
        handles.signifLevel = 0.80;
    case 3
        handles.signifLevel = 0.60;
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on slider movement.
function FreqIncrement_Callback(hObject, eventdata, handles)
% hObject    handle to FreqIncrement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% minimum value to get acceptable accuracy
handles.deltaFreq = round( get(hObject,'Value') * 100 ) / 100; 

set(hObject, 'Value', handles.deltaFreq );
set(handles.FreqIncrementText, 'String', num2str(handles.deltaFreq) );

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in MorletFn.
function MorletFn_Callback(hObject, eventdata, handles)
% hObject    handle to MorletFn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~get(hObject,'Value') 
    msgbox( 'At least one wavelet function must be selected.' );
    set (hObject,'Value',1);
    return;
end

handles.motherWavelet = 'MORLET';     % Use capitals 

% Get Ko
set(handles.Ko, 'Enable', 'on');
handles.wltParameter = get(handles.Ko,'Value')+5;

set(handles.DOG, 'Value' , 0    );
set(handles.DOGs_popup, 'Enable', 'off');

% Update handles structure
guidata(hObject, handles);

function Ko_Callback(hObject, eventdata, handles)
% hObject    handle to Ko (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.wltParameter = get(hObject,'Value')+5;
handles.wltParameter 
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in DOG.
function DOG_Callback(hObject, eventdata, handles)
% hObject    handle to DOG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~get(hObject,'Value') 
    msgbox( 'At least one wavelet function must be selected.' );
    set (hObject,'Value',1);
    return;
end

set(handles.DOGs_popup , 'Enable', 'on');
param = get(handles.DOGs_popup,'Value');

switch param
    case 1
        handles.motherWavelet = 'MEXICANHAT'; % Use capitals
        handles.wltParameter = 2;
    case 2
        handles.motherWavelet = 'DOG'; % Use capitals
        handles.wltParameter = 6;
end

set(handles.MorletFn   , 'Value' , 0    );
set(handles.Ko         , 'Enable', 'off');

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in DOGs_popup.
function DOGs_popup_Callback(hObject, eventdata, handles)
% hObject    handle to DOGs_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

choice = get(hObject,'value');
switch choice
    case 1
        handles.motherWavelet = 'MEXICANHAT';
        handles.wltParameter = 2;
    case 2
        handles.motherWavelet = 'DOG';
        handles.wltParameter = 6;
end
% Update handles structure
guidata(hObject, handles);

% Define output settings
% ------------------------------------------------------------------------

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

% --- Executes on selection change in FigATitle.
function FigATitle_Callback(hObject, eventdata, handles)
% hObject    handle to FigATitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.TitleFigA = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in SignalName_popup.
function SignalName_popup_Callback(hObject, eventdata, handles)
% hObject    handle to SignalName_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
handles.SignalName = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in SampleName_popup.
function SampleName_popup_Callback(hObject, eventdata, handles)
% hObject    handle to SampleName_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
handles.SampleName = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in SignalUnit_popup.
function SignalUnit_popup_Callback(hObject, eventdata, handles)
% hObject    handle to SignalUnit_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
handles.SignalUnit = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in SampleUnit_popup.
function SampleUnit_popup_Callback(hObject, eventdata, handles)
% hObject    handle to SampleUnit_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = get(hObject,'String');
handles.SampleUnit = contents{get(hObject,'Value')};

% Update handles structure
guidata(hObject, handles);

% Run the program
% ------------------------------------------------------------------------

% --- Executes on button press in timetest.
function timetest_Callback(hObject, eventdata, handles)
    
    timetab = {'Without saving figure'       , 0, 0 ;...
               'Saving figure in fig format' , 0, '';...
    	       'Saving figure in jpg format' , 0, '';...
               'Saving figure in tiff format', 0, '';...
               'Saving figure in PDF format' , 0, ''};
    

    wb = waitbar(0,'1','Name','Running Time Test');
    %Time test whitout save
    tic;
    waitbar(1/5,wb,sprintf(['Please wait, don''t open the project folder.']));
    runWltAnalysis(1, 1, handles.ProjectFolder,...
         0, handles.figextension, handles.SignalName,...
         handles.SignalUnit, handles.SampleName, handles.SampleUnit, ...
         handles.TitleFigA, handles.n, handles.XYZcoor, handles.motherWavelet,...
         handles.wltParameter, handles.deltaFreq, handles.signifLevel)
    timetab{1,2} = toc;
    tottime = timetab{1,2}*1.1*handles.numberOfFiles;
    if tottime > 3600
        timetab{1,3} = [ num2str(round(tottime/3600)) ' hours'];
    elseif tottime > 60
        timetab{1,3} = [ num2str(round(tottime/60)) ' minutes'];
    else 
        timetab{1,3} = [ num2str(round(tottime)) ' seconds'];
    end
    formats = {'fig','jpg','tiff','PDF'};
    for i=1:4
        waitbar((i+1)/5,wb,sprintf(['Please wait, don''t open the project folder.']));
        tic;
         runWltAnalysis(1, 1, handles.ProjectFolder,...
              1, formats{i}, handles.SignalName,...
              handles.SignalUnit  , handles.SampleName, handles.SampleUnit, ...
              handles.TitleFigA   , handles.n, handles.XYZcoor, handles.motherWavelet,...
              handles.wltParameter, handles.deltaFreq , handles.signifLevel)
        timetab{i+1,2} = toc;
        tottime = timetab{i+1,2}*1.1*handles.numberOfFiles;
        if tottime > 3600
            timetab{i+1,3} = [ num2str(round(tottime/3600)) ' hours'];
        elseif tottime > 60
            timetab{i+1,3} = [ num2str(round(tottime/60)) ' minutes'];
        else 
            timetab{i+1,3} = [ num2str(round(tottime)) ' seconds'];
        end
    end
    
    rmdir(fullfile(handles.ProjectFolder, 'Wavelet_Output'), 's');
    delete(wb);
    f = figure('Position',[400 300 550 350], 'name', 'Time Test',...
        'WindowStyle', 'modal','NumberTitle','off');
    
    message = sprintf(['Saving the figures can greatly increase the runtime of the program\n(See runtime estimates in the table below).\n'...
                       'We recommend you to print the images after running this application. After doing so, you can use the "plotWltAnalysis" function that can be called from the MATLAB console.\n'...
                       'This function needs a string parameter (e.g., ''fig'', ''jpg'', ''tiff'' or ''PDF''), depending on the desired output format and prompts you to select the profiles to print (i.e. the .MAT files in the folder Data_Output within the project folder).\n'...
                       'The following table shows an estimate of the runtime:']);
%     
    uicontrol(f,'Units','normalized','style', 'text', 'string', message, 'position', [0.05 0.50 0.90 0.45]);
    columnname  = {'Type', 'Runtime per profile (s) ', 'Total estimated runtime'};
    columnwidth = {180,130,170};
    columnformat = {'char', 'char', 'numeric'}; 
    
    uitable('Units','normalized', 'Data', timetab,... 
            'ColumnName', columnname,'Position', [0.05 0.05 0.91 0.40],...
            'ColumnFormat', columnformat,'ColumnWidth',columnwidth,...
            'RowName',[],'FontSize',10); 
     
% --- Executes on button press in RunTheProgram.
function RunTheProgram_Callback(hObject, eventdata, handles)
% hObject    handle to RunTheProgram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Verify output parameters if user will save the outputdata
if handles.SaveFigure
    if (handles.TitleFigA  == '-') 
        msgbox( 'Title Fig has not been defined.' );
        return;
    elseif(handles.SignalName == '-')
        msgbox( 'Signal Name has not been defined.' );
        return;
    elseif  (handles.SampleName == '-')
        msgbox( 'Sample Name has not been defined.' );
        return;
    elseif (handles.SignalUnit == '-')
        msgbox( 'Signal Unit has not been defined.' );
        return;
    elseif (handles.SampleUnit == '-')
        msgbox( 'Sample Unit has not been defined.' );
        return;
    elseif (handles.numberOfFiles <= 0)
        msgbox( 'No data to process.' );
        return;
    end
end
    
if ~handles.SaveFigure
    choice =  questdlg( 'Are you sure you do not want to save the figure(s)?', ...
                        'Do you want to continue?', ...
                        'Yes', 'No', 'No');
    if strcmp(choice,'No')
        return;
    end
    if (handles.TitleFigA  == '-') 
        handles.TitleFigA = handles.NameProject;
    end
    if(handles.SignalName == '-')
        handles.SignalName = 'Elevation';
    end
    if  (handles.SampleName == '-')
        handles.SampleName = 'X';
    end
    if (handles.SignalUnit == '-')
        handles.SignalUnit = 'm';
    end
    if (handles.SampleUnit == '-')
        handles.SampleUnit = 'm';
    end
end

close;

runWltAnalysis(0, handles.numberOfFiles, handles.ProjectFolder,...
         handles.SaveFigure, handles.figextension, handles.SignalName,...
         handles.SignalUnit, handles.SampleName, handles.SampleUnit, ...
         handles.TitleFigA, handles.n, handles.XYZcoor,...
         handles.motherWavelet, handles.wltParameter, handles.deltaFreq, handles.signifLevel)
display('Successful run!!!');

% Other functions
% ------------------------------------------------------------------------

function [Name] = insert_name(strText, strTitle, preText)
    new_name= inputdlg(strText,strTitle, [1 50], preText);
                                   
    if (isempty(new_name))
        Name ='';
        return
    end    

    Name = new_name{:};
    if (~is_alpha_num(Name)) 
        uiwait(msgbox(sprintf( ['The entered  name has an invalid character.\n' ...
                                'Please introduce a valid name.']),...
                                'Invalid project name','modal'));

        Name = insert_name(strText, strTitle, preText);
    end
    
    
function alphanum = is_alpha_num (s)

alphanum = (all( (s>='0' & s<='9') | ...
                 (s>='a' & s<='z') | ...
                 (s>='A' & s<='Z') ));
