function varargout = BedformsATM_App4_3DAnalysis(varargin)
% BEDFORMSATM_APP4_3DANALYSIS MATLAB code for BedformsATM_App4_3DAnalysis.fig
%      BEDFORMSATM_APP4_3DANALYSIS, by itself, creates a new BEDFORMSATM_APP4_3DANALYSIS or raises the existing
%      singleton*.
%
%      H = BEDFORMSATM_APP4_3DANALYSIS returns the handle to a new BEDFORMSATM_APP4_3DANALYSIS or the handle to
%      the existing singleton*.
%
%      BEDFORMSATM_APP4_3DANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEDFORMSATM_APP4_3DANALYSIS.M with the given input arguments.
%
%      BEDFORMSATM_APP4_3DANALYSIS('Property','Value',...) creates a new BEDFORMSATM_APP4_3DANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BedformsATM_App4_3DAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BedformsATM_App4_3DAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BedformsATM_App4_3DAnalysis

% Last Modified by GUIDE v2.5 19-Nov-2016 13:16:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BedformsATM_App4_3DAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @BedformsATM_App4_3DAnalysis_OutputFcn, ...
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


% --- Executes just before BedformsATM_App4_3DAnalysis is made visible.
function BedformsATM_App4_3DAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BedformsATM_App4_3DAnalysis (see VARARGIN)

handles.SaveFigure     = 0;
handles.figextension   = 'PDF';

% Both must be the same
handles.NameProject = 0;

handles.Select_Input = 4; % All hierachies by default

handles.Plot_Domain = [];
handles.Plot_Window = [];

set(findall(handles.Panel_H,'Enable','On'),'Enable','Off');

% Choose default command line output for BedformsATM_App4_3DAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = BedformsATM_App4_3DAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in SelectWavelet.
function SelectDiscrimination_Callback(hObject, eventdata, handles)
% hObject    handle to SelectWavelet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(hObject, 'BackGroundColor', [1 0.5 0]);
handles.NameProject.Discrimination = 1;

% Get projects path for uipickfiles
DocumentsPath = pwd;

% Get .MAT for_Analysis from Discrimination
handles.projectfolder = fullfile(DocumentsPath, 'Projects');
loaded_path = uipickfiles('FilterSpec',handles.projectfolder, ...
                           'NumFiles', 1, 'Output', 'char');
      
if ( loaded_path==0 )
    msgbox('No file had been loaded.');
    return;
end

load (loaded_path);

% Transformation of coordinates

deltaX = norm ( [FileBed_dataMX(1,1), FileBed_dataMY(1,1)] - [FileBed_dataMX(2,1), FileBed_dataMY(2,1)] );
deltaY = norm ( [FileBed_dataMX(1,1), FileBed_dataMY(1,1)] - [FileBed_dataMX(1,2), FileBed_dataMY(1,2)] );

[Xsize,Ysize] = size (FileBed_dataMX);

handles.FileBed.X = [0:Xsize-1]'*deltaX;
handles.FileBed.X = repmat(handles.FileBed.X, 1, Ysize);

handles.FileBed.Y = [0:Ysize-1]*deltaY;
handles.FileBed.Y = repmat(handles.FileBed.Y, Xsize, 1);

handles.FileBed.Z = FileBed_dataMZ;
handles.FileBed.Z13 = FileBed_dataMZ13;
handles.FileBed.Z23 = FileBed_dataMZ23;
handles.FileBed.Z33 = FileBed_dataMZ33;

handles.WindowData.Base_L(1) = BFPScale1;
handles.WindowData.Base_L(3) = BFPScale2;
handles.WindowData.Base_L(4) = BFPScale2;


pos_H23_analysis = [1 2] * round (Ysize/3);
pos_H23_analysis = [1 pos_H23_analysis Ysize];

BFPScale = 0;

for k =1:4
    Z = FileBed_dataMZ23(:, pos_H23_analysis(k));
    
    [wltPower, fourierPeriod, ~] =...
     waveletTransform (motherWavelet,wltParameter, Z', deltaX   ,deltaFreq );
    %% Statistical analysis
    lowerScale = 0.001;
    upperScale = 10000;

    [~, globalSignif,globalWs,~,~,~] =...
    statisticsWlt(motherWavelet, wltParameter, Z',deltaX,...
    deltaFreq, signifLevel, wltPower,lowerScale, upperScale);

    fun = log10(globalWs);
    
    m = mean(fun);
    sd = std(fun);
    fun = (fun-m)/sd;
    
    param = 0.0001*(max(fun)-min(fun));
    [maxp,~] = peakdet( fun, param );
    
    [x,~] = size(maxp);
    target = BFPScale2;
    for i=1:x
        % Just if the peak pass the confidence line and is inside the cone
        % of influence
        if ( globalWs( maxp(i,1) ) >= globalSignif( maxp(i,1) ) ) && fourierPeriod(maxp(i,1)) <= BFPScale2
            target = fourierPeriod(maxp(i,1));
        end
    end
    BFPScale = BFPScale + target;
end

handles.WindowData.Base_L(2) = BFPScale/4;

% Get the name of the project
aux = find(loaded_path==filesep);
pos = strfind ( loaded_path, [filesep 'Projects' filesep] );
pos = find(aux == pos) + 1;
handles.NameProject = loaded_path((aux(pos)+1):(aux(pos+1)-1));

[~,~,~] = mkdir ( [handles.projectfolder filesep handles.NameProject filesep '3DAnalysis'] );

handles.output_path = [handles.projectfolder filesep handles.NameProject filesep '3DAnalysis' filesep];

set(findall(handles.Panel_H,'Enable','Off'),'Enable','On');
set(findall(handles.Panel_Input,'Enable','On'),'Enable','Inactive');
set(handles.Input0, 'Value', 1);

set(hObject, 'BackGroundColor', 'g');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Input0.
function Input0_Callback(hObject, eventdata, handles)
% hObject    handle to Input0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set (handles.Input0, 'Value', 1);
    set (handles.Input1, 'Value', 0);
    set (handles.Input2, 'Value', 0);
    set (handles.Input3, 'Value', 0);

    handles.Select_Input = 4;
    
guidata(hObject, handles);

% --- Executes on button press in Input1.
function Input1_Callback(hObject, eventdata, handles)
% hObject    handle to Input0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set (handles.Input0, 'Value', 0);
    set (handles.Input1, 'Value', 1);
    set (handles.Input2, 'Value', 0);
    set (handles.Input3, 'Value', 0);

    handles.Select_Input = 1;
    
guidata(hObject, handles);

% --- Executes on button press in Input2.
function Input2_Callback(hObject, eventdata, handles)
% hObject    handle to Input0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set (handles.Input0, 'Value', 0);
    set (handles.Input1, 'Value', 0);
    set (handles.Input2, 'Value', 1);
    set (handles.Input3, 'Value', 0);

    handles.Select_Input = 2;
    
guidata(hObject, handles);

% --- Executes on button press in Input3.
function Input3_Callback(hObject, eventdata, handles)
% hObject    handle to Input0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    set (handles.Input0, 'Value', 0);
    set (handles.Input1, 'Value', 0);
    set (handles.Input2, 'Value', 0);
    set (handles.Input3, 'Value', 1);

    handles.Select_Input = 3;
    
guidata(hObject, handles);

% --- Executes on button press in LoadData.
function LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.Axes_Preview, 'Visible', 'On');
set (handles.Spe_Domain, 'Value', 0);
set (handles.All_Domain, 'Value', 0);

axes(handles.Axes_Preview); hold off;

handles.Plot_Domain = [];
handles.Plot_Window = [];

if     (handles.Select_Input == 1) 
    H = handles.FileBed.Z13;
elseif (handles.Select_Input == 2) 
    H = handles.FileBed.Z23;
elseif (handles.Select_Input == 3) 
    H = handles.FileBed.Z33;
else 
    H = handles.FileBed.Z;
end
colormap gray;
surf(handles.FileBed.X,handles.FileBed.Y,H,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud','SpecularStrength',0);
hold on;
axis off;

daspect([5 5 1]);
axis tight;
view(3);
camlight left;
view(2);

Xlim = get(handles.Axes_Preview, 'XLim');
MinX = Xlim(1);
MaxX = Xlim(2);

Ylim = get(handles.Axes_Preview, 'YLim');
MinY = Ylim(1);
MaxY = Ylim(2);

Zlim = get(handles.Axes_Preview, 'ZLim');
handles.MaxZ = Zlim(2);

handles.Size_X = MaxX - MinX;
handles.Size_Y = MaxY - MinY;

handles.proportion_sides = handles.Size_Y/handles.Size_X;

handles.main_Lx = min (handles.WindowData.Base_L(handles.Select_Input), min(handles.Size_X, handles.Size_Y)/4) ;
handles.main_Ly = handles.main_Lx;

set( handles.SL_WL, 'Min'   , min ( min(handles.Size_X, handles.Size_Y)/4, handles.main_Lx )  );
set( handles.SL_WL, 'Max'   , min ( min(handles.Size_X, handles.Size_Y)/3, handles.main_Lx*3 ) );
set( handles.SL_WL, 'Value' , handles.main_Lx );
guidata(hObject, handles);

% --- Executes on button press in All_Domain.
function All_Domain_Callback(hObject, eventdata, handles)
% hObject    handle to All_Domain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set (hObject, 'Value', 1);
set (handles.Spe_Domain, 'Value', 0);

% Domain is based on window size
% Window size could vary between Lx and 3Lx
% Lx : lenght of the main entity in the analizad data

handles.Lx = handles.main_Lx;
handles.Ly = handles.main_Ly;

handles.Xo = handles.Lx/2;
handles.Yo = handles.Ly/2;

handles.Dx = handles.Size_X - handles.Lx;
handles.Dy = handles.Size_Y - handles.Ly;

guidata(hObject, handles);
refresh_Sliders(hObject,1);

% --- Executes on button press in Spe_Domain.
function Spe_Domain_Callback(hObject, eventdata, handles)
% hObject    handle to Spe_Domain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set (hObject, 'Value', 1);
set (handles.All_Domain, 'Value', 0);

% Agregar una restriccion
handles.Lx = handles.main_Lx;
handles.Ly = handles.main_Ly;

handles.Xo = handles.Lx/2;
handles.Yo = handles.Ly/2;

handles.Dx = min (2*handles.Lx, handles.Size_X - handles.Lx);
handles.Dy = min (2*handles.Ly, handles.Size_Y - handles.Ly);

guidata(hObject, handles);
refresh_Sliders(hObject, 0);

function refresh_Sliders(hObject, type_analysis)
    
    handles  = guidata(hObject);
    
    % type_analysis = 1 : All domain, variable window
    % type_analysis = 0 : Custom domain, static window
    
    set( handles.SL_WL, 'Enable', 'Off' );
    
    set( handles.SL_Xo, 'Enable', 'Off' );
    set( handles.SL_Yo, 'Enable', 'Off' );
    
    set( handles.SL_DL, 'Enable', 'Off' );
    
    if ( type_analysis )
        set( handles.SL_WL, 'Value', handles.Lx );
        set( handles.SL_WL, 'Enable', 'On' );
    else
        
        minXo = handles.Lx/2;
        maxXo = handles.Size_X - handles.main_Lx/2 - handles.Dx;
        
        minYo = handles.Ly/2;
        maxYo = handles.Size_Y - handles.main_Ly/2 - handles.Dy;
        
        minDx = 2*handles.main_Lx;
        
        maxDx = min ( 4*handles.main_Lx, handles.Size_X - handles.main_Lx/2 - handles.Xo);
        
        maxDy = min ( 4*handles.main_Ly, handles.Size_Y - handles.main_Ly/2 - handles.Yo);
        maxDx = min (maxDx, maxDy);
        
        set( handles.SL_Xo, 'Min', minXo );
        set( handles.SL_Xo, 'Max', maxXo );
        set( handles.SL_Xo, 'Value', handles.Xo );
        
        set( handles.SL_Yo, 'Min', minYo );
        set( handles.SL_Yo, 'Max', maxYo );
        set( handles.SL_Yo, 'Value', handles.Yo );
        
        set( handles.SL_DL, 'Min', minDx );
        set( handles.SL_DL, 'Max', maxDx );
        set( handles.SL_DL, 'Value', handles.Dx );
        
        set( handles.SL_Xo, 'Enable', 'On' );
        set( handles.SL_Yo, 'Enable', 'On' );
        set( handles.SL_DL, 'Enable', 'On' );
    end
    
guidata(hObject, handles);
refresh_Text(hObject);
    
function refresh_Text ( hObject )
    handles  = guidata(hObject);
    
    set (handles.T_Lx, 'String', sprintf ('%.2f', handles.Lx));
    set (handles.T_Ly, 'String', sprintf ('%.2f', handles.Ly));
    set (handles.T_Xo, 'String', sprintf ('%.2f', handles.Xo));
    set (handles.T_Yo, 'String', sprintf ('%.2f', handles.Yo));
    set (handles.T_DLx, 'String', sprintf ('%.2f', handles.Dx));
    set (handles.T_DLy, 'String', sprintf ('%.2f', handles.Dy));
    guidata (hObject, handles);
    
surf_analysis(hObject);
    
function surf_analysis(hObject)

    handles  = guidata(hObject);
        
    X = [handles.Xo, handles.Xo + handles.Lx; ...
         handles.Xo, handles.Xo + handles.Lx];
    
    Y = [handles.Yo             , handles.Yo; ...
         handles.Yo + handles.Ly, handles.Yo + handles.Ly];

    Z = ones(2) * handles.MaxZ + 1;
    if isempty(handles.Plot_Window)
        handles.Plot_Window = surf ( X,Y,Z, 'FaceColor','r', ...
        'EdgeColor','none','FaceLighting','gouraud','FaceAlpha',0.8);
    else
        handles.Plot_Window.XData = X;
        handles.Plot_Window.YData = Y;
        handles.Plot_Window.ZData = Z;
    end
    
    X = [handles.Xo, handles.Xo + handles.Dx; ...
         handles.Xo, handles.Xo + handles.Dx];
    
    Y = [handles.Yo             , handles.Yo; ...
         handles.Yo + handles.Dy, handles.Yo + handles.Dy];
    
    Z = ones(2) * handles.MaxZ;
    if isempty(handles.Plot_Domain)
        handles.Plot_Domain = surf ( X,Y,Z, 'FaceColor','b', ...
        'EdgeColor','none','FaceLighting','gouraud','FaceAlpha',0.4);
    else
        handles.Plot_Domain.XData = X;
        handles.Plot_Domain.YData = Y;
        handles.Plot_Domain.ZData = Z;
    end
    axis tight;
guidata (hObject, handles);

% --- Executes on slider movement.
function SL_Xo_Callback(hObject, eventdata, handles)
% hObject    handle to SL_Xo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.Xo = get (hObject, 'Value');
    guidata (hObject, handles);
refresh_Sliders(hObject, 0);
    
    % --- Executes on slider movement.
function SL_Yo_Callback(hObject, eventdata, handles)
% hObject    handle to SL_Yo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.Yo = get (hObject, 'Value');
    guidata (hObject, handles);
refresh_Sliders(hObject, 0);

% --- Executes on slider movement.
function SL_DL_Callback(hObject, eventdata, handles)
% hObject    handle to SL_DL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.Dx = get (hObject, 'Value');
    handles.Dy = handles.Dx;
    guidata (hObject, handles);
refresh_Sliders(hObject, 0);

% --- Executes on slider movement.
function SL_WL_Callback(hObject, eventdata, handles)
% hObject    handle to SL_WL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.Lx = get (hObject, 'Value');
    handles.Ly = handles.Lx;
    
    handles.Xo = handles.Lx/2;
    handles.Yo = handles.Ly/2;

    handles.Dx = handles.Size_X - handles.Lx;
    handles.Dy = handles.Size_Y - handles.Ly;
    
    guidata (hObject, handles);
refresh_Sliders(hObject, 1);

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

% --- Executes on button press in Time_Test.
function Time_Test_Callback(hObject, eventdata, handles)
% hObject    handle to Time_Test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if (handles.Select_Input == 4) DataZ  = handles.FileBed.Z; end
    if (handles.Select_Input == 1) DataZ  = handles.FileBed.Z13; end
    if (handles.Select_Input == 2) DataZ  = handles.FileBed.Z23; end
    if (handles.Select_Input == 3) DataZ  = handles.FileBed.Z33; end

    totTime = run3DAnalysis (handles.Lx, handles.Ly, handles.Xo, handles.Yo, ...
              handles.Dx, handles.Dy, handles.FileBed.X, handles.FileBed.Y,...
              DataZ, 1, handles.SaveFigure, handles.figextension, []);
          
    TimeText = '';
    if (totTime>3600) 
        TimeText = [TimeText ' ' num2str( floor(totTime/3600) ) ' hours'];
        totTime = rem(totTime,3600);
    end
    
    if (totTime>60) 
        TimeText = [TimeText ' ' num2str( floor(totTime/60) ) ' minutes'];
        totTime = rem(totTime,60);
    end
    
    TimeText = [TimeText ' ' num2str( floor(totTime) ) ' seconds'];
          
    msgbox( [ 'The estimated runtime is ' TimeText '.']);
          
% Update handles structure
guidata(hObject, handles);   

% --- Executes on button press in RunTheProgram.
function RunTheProgram_Callback(hObject, eventdata, handles)
% hObject    handle to RunTheProgram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    
    if (handles.Select_Input == 4) 
        DataZ  = handles.FileBed.Z; 
        fileName = 'Whole data';
    end
    if (handles.Select_Input == 1) 
        DataZ  = handles.FileBed.Z13;
        fileName = 'H13';
    end
    if (handles.Select_Input == 2) 
        DataZ  = handles.FileBed.Z23; 
        fileName = 'H23';
    end
    if (handles.Select_Input == 3) 
        DataZ  = handles.FileBed.Z33; 
        fileName = 'H33';
    end
    
    if (get(handles.All_Domain, 'Value') )
        fileName = [ fileName ' - All Domain' ];
    else
        fileName = [ fileName ' - Specific Domain' ];
    end
    
    run3DAnalysis (handles.Lx, handles.Ly, handles.Xo, handles.Yo, ...
    handles.Dx, handles.Dy, handles.FileBed.X, handles.FileBed.Y,...
    DataZ, 0, handles.SaveFigure, handles.figextension, [handles.output_path fileName]);
display('Successful run!!!');
