function varargout = BedformsATM_App1_WaveletAnalysis_ImportFiles(varargin)
% BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES MATLAB code for BedformsATM_App1_WaveletAnalysis_ImportFiles.fig
%      BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES, by itself, creates a new BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES or raises the existing
%      singleton*.
%
%      H = BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES returns the handle to a new BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES or the handle to
%      the existing singleton*.
%
%      BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES.M with the given input arguments.
%
%      BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES('Property','Value',...) creates a new BEDFORMSATM_APP1_WAVELETANALYSIS_IMPORTFILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BedformsATM_App1_WaveletAnalysis_ImportFiles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BedformsATM_App1_WaveletAnalysis_ImportFiles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BedformsATM_App1_WaveletAnalysis_ImportFiles

% Last Modified by GUIDE v2.5 19-Nov-2016 13:08:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BedformsATM_App1_WaveletAnalysis_ImportFiles_OpeningFcn, ...
                   'gui_OutputFcn',  @BedformsATM_App1_WaveletAnalysis_ImportFiles_OutputFcn, ...
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


% --- Executes just before BedformsATM_App1_WaveletAnalysis_ImportFiles is made visible.
function BedformsATM_App1_WaveletAnalysis_ImportFiles_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BedformsATM_App1_WaveletAnalysis_ImportFiles (see VARARGIN)

% handles.ProjectFolder is the path to save input data
% The path to the project's folder is recived as input

handles.ProjectFolder = varargin{1}{1};

% handles.inputformat   Specified the format of the input data
%                   0 = N, X, Y Format (used as default)
%                   1 = X, Y, Z Linear Format
%                   2 = X, Y, Z Curve Format
%                   3 = X, Y    Format

% As default the input format is N, X, Y

set(handles.NXY,'Value',1);
handles.inputformat = 0;

% handles.verticalformat   Specified the format of the vertical coordinate
%                      0 = Meters above the sea level (MASL) or relative coordinate
%                      1 = Water Depth format

% As default the vertical format is MASL
set(handles.WaterDepth,'Value',1);
handles.verticalformat = 0;

% handles.fileformat    Specified the file format
%                   1 = .mat file
%                   2 = text file

% As default the file format is TextFile
set(handles.textFile, 'Value', 1);
handles.fileformat = 2;


% As default the ordinal format is not enabled 
handles.ordinalFormat = 0;

% Choose default command line output for BedformsATM_App1_WaveletAnalysis_ImportFiles
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

axes(handles.axes1);
imshow('RecPlot.png');

axes(handles.axes2);
imshow('CurPlot.png');

% UIWAIT makes BedformsATM_App1_WaveletAnalysis_ImportFiles wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BedformsATM_App1_WaveletAnalysis_ImportFiles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;

% Ordinal format
% ------------------------------------------------------------------------

% --- Executes on button press in OrdinalFormat.
function OrdinalFormat_Callback(hObject, eventdata, handles)
    
% hObject    handle to OrdinalFormat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ordinalFormat = get(hObject,'Value');

% Update handles structure
guidata(hObject, handles);

% Input format
% ------------------------------------------------------------------------

% --- Executes on button press in XYZLinear.
function XYZLinear_Callback(hObject, eventdata, handles)
% hObject    handle to XYZLinear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~get(hObject,'Value') 
    msgbox( 'At least one input format must be selected.' );
    set (hObject,'Value',1);
    return;
end

handles.inputformat = 1;

set(handles.XYZCurve  , 'Value', 0);
set(handles.XY        , 'Value', 0);
set(handles.NXY       , 'Value', 0);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in XYZCurve.
function XYZCurve_Callback(hObject, eventdata, handles)
% hObject    handle to XYZCurve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of XYZCurve

if ~get(hObject,'Value') 
    msgbox( 'At least one input format must be selected.' );
    set (hObject,'Value',1);
    return;
end

handles.inputformat = 2;

set(handles.XYZLinear , 'Value', 0);
set(handles.XY        , 'Value', 0);
set(handles.NXY       , 'Value', 0);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in XY.
function XY_Callback(hObject, eventdata, handles)
% hObject    handle to XY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of XY

if ~get(hObject,'Value') 
    msgbox( 'At least one input format must be selected.' );
    set (hObject,'Value',1);
    return;
end 

handles.inputformat = 3;

set(handles.XYZLinear , 'Value', 0);
set(handles.XYZCurve  , 'Value', 0);
set(handles.NXY       , 'Value', 0);

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in NXY.
function NXY_Callback(hObject, eventdata, handles)
% hObject    handle to NXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NXY

if ~get(hObject,'Value') 
    msgbox( 'At least one input format must be selected.' );
    set (hObject,'Value',1);
    return;
end

handles.inputformat = 0;

set(handles.XYZLinear , 'Value', 0);
set(handles.XYZCurve  , 'Value', 0);
set(handles.XY        , 'Value', 0);

% Update handles structure
guidata(hObject, handles);

% Vertical format
% ------------------------------------------------------------------------

% --- Executes on button press in MASL.
function MASL_Callback(hObject, eventdata, handles)
% hObject    handle to MASL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MASL

if ~get(hObject,'Value') 
    msgbox( 'At least one vertical coordinate format must be selected.' );
    set (hObject,'Value',1);
    return;
end

handles.verticalformat = 0;

set(handles.WaterDepth , 'Value', 0);


% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in WaterDepth.
function WaterDepth_Callback(hObject, eventdata, handles)
% hObject    handle to WaterDepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of WaterDepth

if ~get(hObject,'Value') 
    msgbox( 'At least one vertical coordinate format must be selected.' );
    set (hObject,'Value',1);
    return;
end

handles.verticalformat = 1;

set(handles.MASL , 'Value', 0);

% Update handles structure
guidata(hObject, handles);

% File format
% ------------------------------------------------------------------------

% --- Executes on button press in matFile.
function matFile_Callback(hObject, eventdata, handles)
% hObject    handle to matFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~get(hObject,'Value') 
    msgbox( 'At least one file format must be selected.' );
    set (hObject,'Value',1);
    return;
end

handles.fileformat = 1;
handles.ordinalFormat = 0;

set(handles.OrdinalFormat, 'Value', 0, 'Enable', 'off');
set(handles.textFile  , 'Value', 0);

% Update handles structure
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of matFile

% --- Executes on button press in textFile.
function textFile_Callback(hObject, eventdata, handles)
% hObject    handle to textFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~get(hObject,'Value') 
    msgbox( 'At least one file format must be selected.' );
    set (hObject,'Value',1);
    return;
end

handles.fileformat = 2;

set(handles.OrdinalFormat,  'Enable', 'on');

set(handles.matFile  , 'Value', 0);

% Update handles structure
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of textFile

% Select Data
% ------------------------------------------------------------------------

% --- Executes on button press in SelectFiles.
function SelectFiles_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.signalFile = uipickfiles;

% In case of a MAT file loadedSignalFile will store a cellarray
% containing all input profiles. 
% In case of text files handles.signalFiles will store the paths to files
switch handles.fileformat    
    case 1
        loaded = load( handles.signalFile{1} );
        var    = fieldnames(loaded);
        loadedSignalFile = loaded.( var{1} );
        % handles.numberOfFiles is the quantity of selected files
        [~,handles.numberOfFiles] = size(loadedSignalFile);
    case 2
        % handles.numberOfFiles is the size of the cell array
        [~,handles.numberOfFiles] = size(handles.signalFile);
end

% For XY and NXY format ask for the distance between profiles
if (handles.inputformat==0 || handles.inputformat==3)
    number = inputdlg(sprintf('Plese enter the distance between consecutive profiles:'),...
                 'Distance', [1 50]);
    distprof = str2double(number{:});
end

% Create input folder
path = fullfile(handles.ProjectFolder, 'input');
[~,~,~] = mkdir(path);

% Input data will be saved as a cell array
outFile = cell(1,handles.numberOfFiles);

% Save plant view data for X,Y,Z Curve format prev and hovmoller
XYZ = cell (1,handles.numberOfFiles);

datasize = zeros(1,handles.numberOfFiles);
sampleFreq = zeros(1,handles.numberOfFiles);

for NFiles = 1: handles.numberOfFiles
    
    % Load file
    switch handles.fileformat
        case 1
            profileNum = NFiles;            
            allData = loadedSignalFile{profileNum}; % Load one cell
        case 2
            idData = fopen(handles.signalFile{NFiles},'r');
            
            % Check if the ordinal format is enablead
            % If it is enabled read profileNum as the first number
            % Otherwise profileNum is NFiles ( loop variable )
            if handles.ordinalFormat == get(hObject,'Value'); 
                profileNum = cell2mat(textscan(idData,'%d trash'));
            else
                profileNum = NFiles;
            end

            % Read from File Text
            % Just when handles.inputformat is 3 (X Y format) the array is [nData x 2] size
            switch handles.inputformat
                case 3
                    allData=cell2mat(textscan(idData,'%f %f'));
                otherwise
                    allData=cell2mat(textscan(idData,'%f %f %f'));
            end
            % Close the opened file
            fclose(idData);
    end
    
    % Number of points in the profile
    [nData,~] = size(allData);
    
    % signalData stores the processed data according format N, X, Y
    signalData = zeros(nData,3);
    
    % Ordinal numbers
    signalData(:,1) = 1:nData;
    
    % Save the size of the input
    datasize(profileNum) = nData;
    
    % sampleFreq is the average distance between 2 consecutive points
    Freq = zeros(1,nData-1);
    XYZ{profileNum} = zeros(nData,3);
    switch handles.inputformat
        case 1 % X, Y, Z Linear Format
            signalData(:,3) = allData (:,3);
            for V = 2: nData
                Freq(V-1) = sqrt((allData(V,1)-allData(V-1,1))*(allData(V,1)-allData(V-1,1)) ...
                               + (allData(V,2)-allData(V-1,2))*(allData(V,2)-allData(V-1,2)) );
            end
            XYZ{profileNum} = allData;
        case 2 % X, Y, Z Curve Format
            signalData(:,3) = allData (:,3);
            for V = 2: nData
                Freq(V-1) = sqrt((allData(V,1)-allData(V-1,1))*(allData(V,1)-allData(V-1,1)) ...
                               + (allData(V,2)-allData(V-1,2))*(allData(V,2)-allData(V-1,2)) );
            end
            XYZ{profileNum} = allData;
        case 3 % X, Y Format
            signalData(:,3) = allData (:,2);
            for V = 2: nData
                Freq(V-1) = allData(V,1) - allData(V-1,1);
            end 
            XYZ{profileNum} (:,1) = allData(:,1);
            XYZ{profileNum} (:,2) = profileNum*distprof;
            XYZ{profileNum} (:,3) = allData(:,2);
        case 0 % N, X, Y Format
            signalData(:,3) = allData (:,3);
            for V = 2: nData
                Freq(V-1) = allData(V,2) - allData(V-1,2); 
            end
            XYZ{profileNum} (:,1) = allData(:,2);
            XYZ{profileNum} (:,2) = profileNum*distprof;
            XYZ{profileNum} (:,3) = allData(:,3);
    end
    
    if ~isEquidistant(Freq)
        msgbox(['Profile ' num2str(NFiles) ' have different sample frequencies']);
        return;
    end
    
    sampleFreq(profileNum)  = mean(Freq);
    signalData(:,2) = (0:nData-1)*sampleFreq(profileNum);
    
    if handles.verticalformat == 1
        signalData(:,3) = abs(signalData(:,3));
    end
    outFile{profileNum} = signalData;
end


% Check the size of all loaded files, all must be equal
if ~all(datasize==datasize(1))
    msgbox( ['Some profiles have different number of points. All profiles must have'...
             ' the same quantity of points.']);
    return;
end

% Check if the frequency is constant in all profiles
if ~isEquidistant(sampleFreq)
    msgbox('Profiles have different sample frequencies');
    return;
end

freq = mean(sampleFreq);
Dformat = (handles.inputformat == 1 || handles.inputformat==2);


% Save data as a mat file
save (fullfile(path, 'inputfile.mat'), 'outFile', 'Dformat', '-mat');
save (fullfile(path, 'XYZdata.mat'), 'XYZ', '-mat');

handles.output = 1;

% Update handles structure
guidata(hObject, handles);

close;

% Display a message with the total number of loaded files and the averged
% sample frequency
msgbox( [int2str( handles.numberOfFiles ) ,' Files had been imported.'...
        'The sample period is ' num2str(freq) '.' ] );



%Auxiliar functions

function isEquid = isEquidistant ( x )
    y = sort(x); 
    m = mean(y);
    sd = std(y);
    v = length(find(y<m-3*sd)) + length(find(y>m+3*sd));
    Q1 = median(y(find(y<median(y)))); 
    Q3 = median(y(find(y>median(y)))); 
    IQR = Q3-Q1; 
    outliers = length(find(y<Q1-3*IQR)) + length(find(y>Q3+3*IQR));
    isEquid = (( outliers< 0.05*length(x) ) || v<0.05*length(x)) & max(x)<4*m;
