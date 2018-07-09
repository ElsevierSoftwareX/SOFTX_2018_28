function varargout = BedformsATM_App2_HovmollerAnalysis(varargin)
% BEDFORMSATM_APP2_HOVMOLLERANALYSIS MATLAB code for BedformsATM_App2_HovmollerAnalysis.fig
%      BEDFORMSATM_APP2_HOVMOLLERANALYSIS, by itself, creates a new BEDFORMSATM_APP2_HOVMOLLERANALYSIS or raises the existing
%      singleton*.
%
%      H = BEDFORMSATM_APP2_HOVMOLLERANALYSIS returns the handle to a new BEDFORMSATM_APP2_HOVMOLLERANALYSIS or the handle to
%      the existing singleton*.
%
%      BEDFORMSATM_APP2_HOVMOLLERANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BEDFORMSATM_APP2_HOVMOLLERANALYSIS.M with the given input arguments.
%
%      BEDFORMSATM_APP2_HOVMOLLERANALYSIS('Property','Value',...) creates a new BEDFORMSATM_APP2_HOVMOLLERANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BedformsATM_App2_HovmollerAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BedformsATM_App2_HovmollerAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BedformsATM_App2_HovmollerAnalysis

% Last Modified by GUIDE v2.5 19-Nov-2016 13:09:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BedformsATM_App2_HovmollerAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @BedformsATM_App2_HovmollerAnalysis_OutputFcn, ...
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


% --- Executes just before BedformsATM_App2_HovmollerAnalysis is made visible.
function BedformsATM_App2_HovmollerAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BedformsATM_App2_HovmollerAnalysis (see VARARGIN)

handles.gui = gcf;

% Choose default command line output for BedformsATM_App2_HovmollerAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BedformsATM_App2_HovmollerAnalysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = BedformsATM_App2_HovmollerAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% Select project to analize
%-------------------------------------------------------------------------

% --- Executes on button press in Select_Project.
function Select_Project_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get projects path for uipickfiles
DocumentsPath = pwd;

% Get .MAT files
handles.projectfolder = fullfile(DocumentsPath, 'Projects');
handles.Files = uipickfiles('FilterSpec',handles.projectfolder);

[~,handles.nfiles] = size(handles.Files);

if ( handles.nfiles==0 )
    msgbox('No files had been loaded.');
    return;
end

load (handles.Files{1},'signifLevel','deltaFreq','motherWavelet');

% Read data to show in the main window
set (handles.FreqIncrText   ,'String',num2str(deltaFreq     ));
set (handles.SignifLevelText,'String',num2str(signifLevel   ));
set (handles.NumProfilesText,'String',num2str(handles.nfiles));
set (handles.WavText        ,'String',motherWavelet);

sig   = [];

wb = waitbar(0,'1','Name','Progress...', 'WindowStyle', 'modal');

v = figure('Visible','off');
for fil = 1:handles.nfiles

    waitbar(fil/handles.nfiles,wb,sprintf(['Loading signal ' num2str(fil) '...']));
    
    % Load data to get significant peaks
    load(handles.Files{fil},'fourierPeriod','wltSignif',...
                            'signalDataX','coneOfInfluence');
    % Contour points    
    [C,~] = contour([1:length(signalDataX)],fourierPeriod,wltSignif,[-99,1]);    
    ind = find(C(1,:)~=1);
    for i = 1:length(ind)
        if floor(C(2,ind(i) ))<1 || floor(C(2,ind(i)))> length(signalDataX) 
            continue;
        end
        % Save points inside the cone of influence
        if ( C(1,ind(i)) < coneOfInfluence(floor(C(2,ind(i) ))))
            sig(end+1) = C(1,ind(i));
        end
    end
end
delete(wb);
close(v);

% Get optiman quantity of bins for the histogram
nbin = calcnbins(sig);

[~,M] = hist( sig,nbin );

% Band scale  <> size of a bin
% Upper scale <> median of peaks
handles.bandlength = round(M(2)-M(1));
handles.upperScale = round(median(sig));
handles.lowerScale = handles.upperScale-handles.bandlength;
handles.maxScale   = round(max(sig));

% Define slider
set(handles.sliderScale,'Enable','on');
set(handles.sliderScale,'Max',round(handles.maxScale-handles.bandlength));
set(handles.sliderScale,'SliderStep',[1/( handles.maxScale-handles.bandlength-1 )...
                                     10/( handles.maxScale-handles.bandlength-1 )]);
set(handles.sliderScale,'Value',handles.upperScale-handles.bandlength);

% Define suggested values
set(handles.SugIntLen,'String',num2str(handles.bandlength));
set(handles.SugMaxUS ,'String',num2str(handles.upperScale));

set(handles.IntLen   ,'Enable','on');
set(handles.IntLen   ,'String',num2str(handles.bandlength  ));

set(handles.TextScales, 'String', [num2str(handles.lowerScale)...
                                  ' - ' num2str(handles.upperScale)]);
                              
set(handles.RunHovmoller, 'Enable', 'on');
                            
% Update handles structure
guidata(hObject, handles);

% Define band scale to analize
%-------------------------------------------------------------------------
function IntLen_Callback(hObject, eventdata, handles)
% hObject    handle to IntLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

newlen = round( str2double(get(hObject,'String')));
if ( newlen > handles.maxScale-2 || newlen < 1 )
    set(hObject,'String',num2str(handles.bandlength));
    msgbox(sprintf(['The length entered is out of the allowed range.\n' ...
           'The value should be between 1 and ' num2str(handles.maxScale-2) '.']));
    return;
end

handles.bandlength = newlen;
handles.lowerScale = min(handles.lowerScale, handles.maxScale - newlen);
handles.upperScale = handles.lowerScale + newlen;

set(hObject,'String',num2str(newlen));

set(handles.sliderScale,'Value',handles.lowerScale);
set(handles.sliderScale,'Max',handles.maxScale-handles.bandlength);
set(handles.sliderScale,'SliderStep',[1/( handles.maxScale-handles.bandlength-1 )...
                                     10/( handles.maxScale-handles.bandlength-1 )]);


set(handles.TextScales, 'String', [num2str(handles.upperScale-handles.bandlength)...
                                  ' - ' num2str(handles.upperScale)]);

% Update handles structure
guidata(hObject, handles);

% --- Executes on slider movement.
function sliderScale_Callback(hObject, eventdata, handles)
% hObject    handle to sliderScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.lowerScale = round(get(hObject,'Value'));
handles.upperScale = handles.lowerScale + handles.bandlength;

set(handles.TextScales, 'String', [num2str(handles.upperScale-handles.bandlength)...
                                  ' - ' num2str(handles.upperScale)]);

% Update handles structure
guidata(hObject, handles);

% Run program
%-------------------------------------------------------------------------

% --- Executes on button press in RunHovmoller.
function RunHovmoller_Callback(hObject, eventdata, handles)
% hObject    handle to RunHovmoller (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% quantity of profiles
jnodes = handles.nfiles;

signifLevel = str2double(get(handles.SignifLevelText,'String'));
deltaFreq   = str2double(get(handles.FreqIncrText,'String'));

load(handles.Files{1},'signalDataX','TitleFig','sampleName','sampleUnit');
inodes = length(signalDataX);

path = handles.Files{1};
pos  = find(path==filesep);
% create new folder to hovmoller analysis
handles.fpath = fullfile(path(1:pos(length(pos)-1)-1), 'Hovmoller_Analysis');
[~,~,~] =mkdir(handles.fpath);
handles.fpath = fullfile(handles.fpath, ['Scale_' get(handles.TextScales,'String')]);

% arrays to save data
FileBed_dataMX = zeros(inodes,jnodes); 
FileBed_dataMY = zeros(inodes,jnodes); 
FileBed_dataMZ = zeros(inodes,jnodes);

Lower_Scale = handles.lowerScale;
Upper_Scale = handles.upperScale;

FileBed_dataMV  = zeros(inodes,jnodes);
FileBed_dataMVS = zeros(inodes,jnodes);

wb = waitbar(0,'1','Name','Progress...',...
        'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');
setappdata(wb,'canceling',0);

for jj=1:jnodes
    if getappdata(wb,'canceling')
        break;
    end
    waitbar(jj/jnodes,wb,sprintf(['Working on signal ' num2str(jj) '...']));
    
    jfile = handles.Files{jj};
    load(jfile,'j', 'signalDataX', 'signalData','n_3D_yJ');
    
    sampleFreq = signalDataX(2)-signalDataX(1);
    FileBed_dataMX(:,j) = n_3D_yJ(:,1);
	FileBed_dataMY(:,j) = n_3D_yJ(:,2);
	FileBed_dataMZ(:,j) = n_3D_yJ(:,3);  
	
    load(jfile, 'fourierPeriod','wltPower','motherWavelet', 'wltParameter');

    [wltSignif, ~,~,~,scaleAveSignif,scaleAvg] =...
    statisticsWlt(motherWavelet, wltParameter, signalData,sampleFreq,...
    deltaFreq, signifLevel, wltPower,Lower_Scale, Upper_Scale);

    scaleAvgt       = transpose(scaleAvg);
	wltSignift      = transpose(wltSignif);

    FileBed_dataMV(:,j)=scaleAvgt;

    for ii=1: inodes
        if (FileBed_dataMV (ii,j)<scaleAveSignif)
            FileBed_dataMVS(ii,j)=0.0;
        else
            FileBed_dataMVS(ii,j)=FileBed_dataMV(ii,j);                
        end
    end
end
delete(wb);

% Plot Hovmoller

colorplot = FileBed_dataMZ;
minz = min(min(FileBed_dataMZ));
maxz = max(max(FileBed_dataMZ));
minh = 0.02; % Minimum value to plot in hovmoller
maxh = mean(FileBed_dataMVS(:))+2*std(FileBed_dataMVS(:));
maxh = round(maxh*100)/100;
[N,M] = size(colorplot);
% Auxiliar values
upperg = 16;
upperc = 80;

for i=1:N
    for j=1:M
        if(FileBed_dataMVS(i,j) <= minh) % if value is less than minh, plot in gray
            colorplot(i,j) = ((FileBed_dataMZ(i,j)-minz)*(upperg-1)/(maxz-minz))+1;
        else % plot in colors
            colorplot(i,j) = ((min(FileBed_dataMVS(i,j),maxh)-minh)*(upperc-(upperg+1))/(maxh-minh))+upperg+1;
        end
    end
end

handles.hovfig = figure;

axes('position',[0.10 0.10 0.8 0.7]);
% Define colors to plot
cg = colormap(gray);
ch = colormap(hsv);
ct = [cg(49:64,:);ch];
colormap(ct);
surf(FileBed_dataMX,FileBed_dataMY,FileBed_dataMZ,colorplot,...
    'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud','SpecularStrength',0);
scales = get(handles.TextScales,'String');
waveletf = get(handles.WavText,'String');
Title  = sprintf(['Power Hovmöller\n' TitleFig ' - ' waveletf '\nScale Band: ' scales ' (' sampleUnit ')']);
% Title
set(get(gca,'Title') ,'String',Title);
set(get(gca,'Title') ,'FontSize',20);
set(get(gca,'Title') ,'FontWeight','bold');
% Axes labels
set(get(gca,'XLabel'),'String',[sampleName ' (' sampleUnit ')']);
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'YLabel'),'String','Y(m)');
set(get(gca,'YLabel'),'FontSize',18);

daspect([5 5 1]);
axis tight;
view(2);
camlight left;
caxis([0 80]);
cm = colorbar('location','southoutside');
set(cm, 'xlim', [17 80]);

% colorbar label and ticks
Xt    = zeros(1,8);
Xt(1) = upperg+1;
Xtl   = zeros(1,8);
Xtl(1)= minh;
for i=2:8
    Xt (i)  = Xt (i-1)+(upperc - upperg - 1)/7;
    Xtl(i)  = Xtl(i-1)+(maxh   - minh      )/7;
end
Xt  = round(Xt *100)/100;
Xtl = round(Xtl*100)/100;
set(cm, 'XTick',     Xt);
set(cm, 'XTickLabel',Xtl);
xlabel(cm,['Ave. variance (' sampleUnit '^2',')'],'fontsize',14);

mTextBox = uicontrol('style','text','Units','normalized','FontSize',8,...
                    'String','Plotted by Bedforms-ATM','BackgroundColor',[1 1 1],...
                    'Position',[0.7 0.05 0.25 0.03]);

set(handles.SaveFig, 'Enable', 'on');
set(handles.ExportFig, 'Enable', 'on');
                
handles.X = FileBed_dataMX;
handles.Y = FileBed_dataMY;
handles.Z = FileBed_dataMZ;
handles.H = FileBed_dataMVS;

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SaveFig.
function SaveFig_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~ishandle(handles.hovfig)
    msgbox('No figure to save.');
    return;
end

str = {'jpg','tiff','PDF'};
[choice,ok] = listdlg('PromptString','Select a format to save:',...
                      'SelectionMode','single',...
                      'ListString',str);
if ~ok
    return;
end
figureformat = str{choice};

set(handles.hovfig,'paperunits','centimeters')
set(handles.hovfig,'papersize',[28,19])
set(handles.hovfig,'paperposition',[0.05,0.05,29*.95,20*.95]) 
switch figureformat
    case 'PDF'
        print (handles.hovfig,'-dpdf','-r300',handles.fpath)
    case 'jpg'
        print (handles.hovfig,'-djpeg','-r300',handles.fpath)
    case 'tiff'
        print (handles.hovfig,'-dtiff','-r300',handles.fpath)
end

% --- Executes on button press in ExportFig.
function ExportFig_Callback(hObject, eventdata, handles)
% hObject    handle to ExportFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~ishandle(handles.hovfig)
    msgbox('No figure to export.');
    return;
end

str = {'TecPlot','img'};
[choice,ok] = listdlg('PromptString','Select a format to export:',...
                      'SelectionMode','single',...
                      'ListString',str);
if ~ok
    return;
end
figureformat = str{choice};

switch figureformat
    case 'img'
        saveas(handles.hovfig,handles.fpath,'fig');
    case 'TecPlot'
        FileTec = [handles.fpath '.dat'];
        FileTecSave=fopen(FileTec,'wt');
        string = ['TITLE = "WAVELET ANALYSIS"','\n'];
        fprintf(FileTecSave, string);
        fprintf(FileTecSave,'VARIABLES = "X" "Y" "Z" "SAS"  \n');
        
        [rowS,colS]=size(handles.X);
        string = ['ZONE T="Band ' num2str(handles.lowerScale) ' to ' num2str(handles.upperScale) '1/m' '"' ' I=' num2str(rowS) ', J=' num2str(colS) ', F=POINT \n'];
        fprintf(FileTecSave, string);
        for j=1: colS
            for i=1: rowS
                fprintf(FileTecSave, '%6.5f %6.5f %6.10f %6.20f \n', handles.X(i,j), handles.Y(i,j), handles.Z(i,j), handles.H(i,j));
            end
        end
        fclose(FileTecSave);
end

% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(handles.gui);



                              
