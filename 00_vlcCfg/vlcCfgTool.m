function varargout = vlcCfgTool(varargin)
% VLCCFGTOOL MATLAB code for vlcCfgTool.fig
%      VLCCFGTOOL, by itself, creates a new VLCCFGTOOL or raises the existing
%      singleton*.

% Last Modified by GUIDE v2.5 27-Jul-2019 16:09:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vlcCfgTool_OpeningFcn, ...
                   'gui_OutputFcn',  @vlcCfgTool_OutputFcn, ...
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

% --- Executes just before vlcCfgTool is made visible.
function vlcCfgTool_OpeningFcn(hObject, eventdata, handles, varargin)
addpath('fcn');

% Choose default command line output for vlcCfgTool
handles.output = hObject;

% Explanation picture for QRX
matlabImage = imread('data/expl.png');
image(handles.explAxis,matlabImage)
axis off
axis image

% create the listener for the lateral movement of TX (in RX frame) slider
handles.latListener = addlistener(handles.lat_slider,'ContinuousValueChange', @(hFigure,eventdata) latContValCallback(hObject,eventdata));
% create the listener for the forward movement of TX (in RX frame) slider
handles.fwdListener = addlistener(handles.fwd_slider,'ContinuousValueChange', @(hFigure,eventdata) fwdContValCallback(hObject,eventdata));
% create the listener for the heading of the TX (in RX frame) slider
handles.hdgListener = addlistener(handles.hdg_slider,'ContinuousValueChange', @(hFigure,eventdata) hdgContValCallback(hObject,eventdata));

% Set slider steps 
set(handles.lat_slider, 'SliderStep', [0.001, 0.01]);
set(handles.fwd_slider, 'SliderStep', [0.001, 0.01]);
set(handles.hdg_slider, 'SliderStep', [0.0005, 0.005]);

% Remove expl, spot, LUT and sys test axes for cooler looks
set(handles.explAxis,'YTickLabel',[]);
set(handles.explAxis,'XTickLabel',[]);
set(handles.spotQrxAxis,'YTickLabel',[]);
set(handles.spotQrxAxis,'XTickLabel',[]);
set(handles.thetaLutAxis,'YTickLabel',[]);
set(handles.thetaLutAxis,'XTickLabel',[]);
set(handles.sysTestAxis,'YTickLabel',[]);
set(handles.sysTestAxis,'XTickLabel',[]);

% initiate polar plot with nothing (this is like a bug/wont-fix situation 
% from MATLAB side so we needed a workaround)
axes(handles.TxPlrPatAxis)
polarplot(0,0);
load('vlcCfgTxPlrPatTool/data/vlcCfgTxPlrPat_stdPwrLed.mat');
ax = gca;
handles.h_TxPlrPatAxis = ax;
handles.txPlrPatNum = vlcCfgTxPlrPatNum;
handles.txPlrPatFit = vlcCfgTxPlrPatFit; % for later, saving it to file
polarplot(deg2rad(vlcCfgTxPlrPatNum(1,:)),vlcCfgTxPlrPatNum(2,:),'g','LineWidth',2)
hold on
handles.h_plr_L = plot(handles.h_TxPlrPatAxis,0,0,'b','LineWidth',1);
handles.h_plr_R = plot(handles.h_TxPlrPatAxis,0,0,'r','LineWidth',1);
hold off
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
ax.ThetaLim = [-90 90];
ax.RTick = [0:0.2:1];

% flag that tells if the fov/LUT is drawn
handles.fov_drawn_flag = 0;

% init graphics
handles.dDL = 03.00e-3;
handles.dSL = 00.25e-3;
handles.dFL = 02.25e-3;
handles.dPX = 01.50e-3;
handles.dRes = 01.00e-5;
handles.lum2pkpk = 10000;
handles.lat = 2;
handles.fwd = 10;
handles.hdg = 0;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
vlcCfgTool_sysTest(handles);
xlabel(handles.thetaLutAxis,'Theta (deg)')
ylabel(handles.thetaLutAxis,'Spot Position (p.u.)')
[~, handles.lutng, handles.lutps] = vlcCfgTool_fovLut(handles);
handles.fov_drawn_flag = 1;
hold(handles.thetaLutAxis,'on')
handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
hold(handles.thetaLutAxis,'off')
if ((handles.lutng>rad2deg(theta_l))||(handles.lutng>rad2deg(theta_r))||(handles.lutps<rad2deg(theta_l))||(handles.lutps<rad2deg(theta_r)))
    set(handles.qrx1Lum_txt,'String','Out of FoV!');
    set(handles.qrx2Lum_txt,'String','Out of FoV!');
else
    [qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
    set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
    set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vlcCfgTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Slider callback function
function latContValCallback(hFigure,eventdata)
handles = guidata(hFigure);
val = get(handles.lat_slider,'Value');
set(handles.lat_txt,'String',num2str(val));

% Slider callback function
function fwdContValCallback(hFigure,eventdata)
handles = guidata(hFigure);
val = get(handles.fwd_slider,'Value');
set(handles.fwd_txt,'String',num2str(val));

% Slider callback function
function hdgContValCallback(hFigure,eventdata)
handles = guidata(hFigure);
val = get(handles.hdg_slider,'Value');
set(handles.hdg_txt,'String',num2str(val));

% --- Outputs from this function are returned to the command line.
function varargout = vlcCfgTool_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TX Load Polar Pattern %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ldTxPlrPat_button_Callback(hObject, eventdata, handles)
[fn, pt] = uigetfile;
load(strcat(pt,fn));
% TODO: check if the file loaded is NOT a txPlrPat array
handles.txPlrPatNum = vlcCfgTxPlrPatNum;
handles.txPlrPatFit = vlcCfgTxPlrPatFit;
polarplot(handles.h_TxPlrPatAxis,deg2rad(vlcCfgTxPlrPatNum(1,:)),vlcCfgTxPlrPatNum(2,:),'g','LineWidth',2)
ax = handles.h_TxPlrPatAxis;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
ax.ThetaLim = [-90 90];
ax.RTick = [0:0.2:1];
hold(handles.h_TxPlrPatAxis,'on');
handles.h_plr_L = plot(handles.h_TxPlrPatAxis,0,0,'b','LineWidth',1);
handles.h_plr_R = plot(handles.h_TxPlrPatAxis,0,0,'r','LineWidth',1);
hold(handles.h_TxPlrPatAxis,'off');
[qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TX LED Max Lumens %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maxLum_edit_Callback(hObject, eventdata, handles)
maxLum = str2double(get(hObject,'String'));
handles.maxLum=maxLum;
[qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
guidata(hObject,handles);
function maxLum_edit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','1000')
maxLum = str2double(get(hObject,'String'));
handles.maxLum=maxLum;
guidata(hObject,handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% QRX Lum2PkPk Conversion %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lum2pkpk_edit_Callback(hObject, eventdata, handles)
lum2pkpk = str2double(get(hObject,'String'));
handles.lum2pkpk=lum2pkpk;
[qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
guidata(hObject,handles);
function lum2pkpk_edit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','10000');
lum2pkpk = str2double(get(hObject,'String'));
handles.lum2pkpk=lum2pkpk;
guidata(hObject,handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% QRX Lens Diameter %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dDL_edit_Callback(hObject, eventdata, handles)
dDL = str2double(get(hObject,'String'));
if(dDL < 2*handles.dPX)
    errordlg('Lens diameter smaller than 2*dPX does not make sense...')
    set(hObject,'String','03.00e-3')
    dDL = str2double(get(hObject,'String'));
end
handles.dDL=dDL;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
[~, handles.lutng, handles.lutps] = vlcCfgTool_fovLut(handles);
handles.fov_drawn_flag = 1;
hold(handles.thetaLutAxis,'on')
handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
hold(handles.thetaLutAxis,'off')
[qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
guidata(hObject,handles);
function dDL_edit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','03.00e-3')
dDL = str2double(get(hObject,'String'));
handles.dDL=dDL;
guidata(hObject,handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% QRX Lens-Sensor Dist %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dSL_edit_Callback(hObject, eventdata, handles)
dSL = str2double(get(hObject,'String'));
handles.dSL=dSL;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
[~, handles.lutng, handles.lutps] = vlcCfgTool_fovLut(handles);
handles.fov_drawn_flag = 1;
hold(handles.thetaLutAxis,'on')
handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
hold(handles.thetaLutAxis,'off')
guidata(hObject,handles);
function dSL_edit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','00.25e-3');
dSL = str2double(get(hObject,'String'));
handles.dSL=dSL;
guidata(hObject,handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% QRX Focal Length %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dFL_edit_Callback(hObject, eventdata, handles)
dFL = str2double(get(hObject,'String'));
handles.dFL=dFL;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
[~, handles.lutng, handles.lutps] = vlcCfgTool_fovLut(handles);
handles.fov_drawn_flag = 1;
hold(handles.thetaLutAxis,'on')
handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
hold(handles.thetaLutAxis,'off')
guidata(hObject,handles);
function dFL_edit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','02.25e-3')
dFL = str2double(get(hObject,'String'));
handles.dFL=dFL;
guidata(hObject,handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% QRX Pixel Dimension %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dPX_edit_Callback(hObject, eventdata, handles)
dPX = str2double(get(hObject,'String'));
handles.dPX=dPX;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
[~, handles.lutng, handles.lutps] = vlcCfgTool_fovLut(handles);
handles.fov_drawn_flag = 1;
hold(handles.thetaLutAxis,'on')
handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
hold(handles.thetaLutAxis,'off')
guidata(hObject,handles);
function dPX_edit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','01.50e-3')
dPX = str2double(get(hObject,'String'));
handles.dPX=dPX;
guidata(hObject,handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% QRX Metric Resolution %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dRes_edit_Callback(hObject, eventdata, handles)
dRes = str2double(get(hObject,'String'));
handles.dRes=dRes;
guidata(hObject,handles);
function dRes_edit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','01.00e-5')
dRes = str2double(get(hObject,'String'));
handles.dRes=dRes;
guidata(hObject,handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% FoV Sweep Resolution %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function thetaSwpRes_edit_Callback(hObject, eventdata, handles)
thetaSwpRes = str2double(get(hObject,'String'));
handles.thetaSwpRes = thetaSwpRes;
guidata(hObject,handles);
function thetaSwpRes_edit_CreateFcn(hObject, eventdata, handles)
set(hObject,'String','0.1')
thetaSwpRes = str2double(get(hObject,'String'));
handles.thetaSwpRes = thetaSwpRes;
guidata(hObject,handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Lateral Slider Callback %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lat_slider_Callback(hObject, eventdata, handles)
lat = get(hObject,'Value');
handles.lat = lat;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
vlcCfgTool_sysTest(handles);
if(handles.fov_drawn_flag == 1)
    hold(handles.thetaLutAxis,'on')
    delete(handles.hL);
	delete(handles.hR);
    handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
    handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
	hold(handles.thetaLutAxis,'off')
end
hold(handles.h_TxPlrPatAxis,'on');
delete(handles.h_plr_L);
delete(handles.h_plr_R);
handles.h_plr_L = polarplot(handles.h_TxPlrPatAxis,[deg2rad(handles.hdg)+theta_l;deg2rad(handles.hdg)+theta_l],[0;1],'b');
handles.h_plr_R = polarplot(handles.h_TxPlrPatAxis,[deg2rad(handles.hdg)+theta_r;deg2rad(handles.hdg)+theta_r],[0;1],'r');
hold(handles.h_TxPlrPatAxis,'off');
if ((handles.lutng>rad2deg(theta_l))||(handles.lutng>rad2deg(theta_r))||(handles.lutps<rad2deg(theta_l))||(handles.lutps<rad2deg(theta_r)))
    set(handles.qrx1Lum_txt,'String','Out of FoV!');
    set(handles.qrx2Lum_txt,'String','Out of FoV!');
else
    [qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
    set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
    set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
end
guidata(hObject,handles);
function lat_slider_CreateFcn(hObject, eventdata, handles)
lat = 2;
set(hObject,'Value',lat);
handles.lat = lat;
guidata(hObject,handles);
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Forward Slider Callback %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fwd_slider_Callback(hObject, eventdata, handles)
fwd = get(hObject,'Value');
handles.fwd = fwd;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
vlcCfgTool_sysTest(handles);
if(handles.fov_drawn_flag == 1)
	hold(handles.thetaLutAxis,'on')
    delete(handles.hL);
	delete(handles.hR);
    handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
    handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
	hold(handles.thetaLutAxis,'off')
end
hold(handles.h_TxPlrPatAxis,'on');
delete(handles.h_plr_L);
delete(handles.h_plr_R);
handles.h_plr_L = polarplot(handles.h_TxPlrPatAxis,[deg2rad(handles.hdg)+theta_l;deg2rad(handles.hdg)+theta_l],[0;1],'b');
handles.h_plr_R = polarplot(handles.h_TxPlrPatAxis,[deg2rad(handles.hdg)+theta_r;deg2rad(handles.hdg)+theta_r],[0;1],'r');
hold(handles.h_TxPlrPatAxis,'off');
if ((handles.lutng>rad2deg(theta_l))||(handles.lutng>rad2deg(theta_r))||(handles.lutps<rad2deg(theta_l))||(handles.lutps<rad2deg(theta_r)))
    set(handles.qrx1Lum_txt,'String','Out of FoV!');
    set(handles.qrx2Lum_txt,'String','Out of FoV!');
else
    [qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
    set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
    set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
end    
guidata(hObject,handles);
function fwd_slider_CreateFcn(hObject, eventdata, handles)
fwd = 10;
set(hObject,'Value',fwd);
handles.fwd = fwd;
guidata(hObject,handles);
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Heading Slider Callback %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdg_slider_Callback(hObject, eventdata, handles)
hdg = get(hObject,'Value');
handles.hdg = hdg;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
vlcCfgTool_sysTest(handles);
if(handles.fov_drawn_flag == 1)
    hold(handles.thetaLutAxis,'on')
    delete(handles.hL);
	delete(handles.hR);
    handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
    handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
	hold(handles.thetaLutAxis,'off')
end
hold(handles.h_TxPlrPatAxis,'on');
delete(handles.h_plr_L);
delete(handles.h_plr_R);
handles.h_plr_L = polarplot(handles.h_TxPlrPatAxis,[deg2rad(hdg)+theta_l;deg2rad(hdg)+theta_l],[0;1],'b');
handles.h_plr_R = polarplot(handles.h_TxPlrPatAxis,[deg2rad(hdg)+theta_r;deg2rad(hdg)+theta_r],[0;1],'r');
hold(handles.h_TxPlrPatAxis,'off');
if ((handles.lutng>rad2deg(theta_l))||(handles.lutng>rad2deg(theta_r))||(handles.lutps<rad2deg(theta_l))||(handles.lutps<rad2deg(theta_r)))
    set(handles.qrx1Lum_txt,'String','Out of FoV!');
    set(handles.qrx2Lum_txt,'String','Out of FoV!');
else
    [qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
    set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
    set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
end
guidata(hObject,handles);
function hdg_slider_CreateFcn(hObject, eventdata, handles)
hdg = 0;
set(hObject,'Value',hdg);
handles.hdg = hdg;
guidata(hObject,handles);
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% LUT/FOV Calc Callback %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function calcFov_button_Callback(hObject, eventdata, handles)
[~, handles.lutng, handles.lutps] = vlcCfgTool_fovLut(handles);
handles.fov_drawn_flag = 1;
hold(handles.thetaLutAxis,'on')
handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
hold(handles.thetaLutAxis,'off')
guidata(hObject,handles);
function calcFov_button_CreateFcn(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Config Save Callback %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cfgSv_button_Callback(hObject, eventdata, handles)
vlcCfg_rx_dDL = handles.dDL;
vlcCfg_rx_dSL = handles.dSL;
vlcCfg_rx_dFL = handles.dFL;
vlcCfg_rx_dPX = handles.dPX;
[vlcCfg_rx_thetaLut, vlcCfg_rx_dRes, vlcCfg_rx_thetaSwpRes]  = vlcCfgTool_rxThetaLutFit(handles);
vlcCfg_rx_lum2pkpk = handles.lum2pkpk;
vlcCfg_tx_maxLum = handles.maxLum;
vlcCfg_tx_plrPatFit = handles.txPlrPatFit;
vlcCfg_tx_plrPatNum = handles.txPlrPatNum;
answer = inputdlg('Enter filename for the config file','Config Filename',[1 50],{'vlcCfg_<explanation>.mat'});
save(strcat('data/',answer{1}),'vlcCfg_rx_dDL','vlcCfg_rx_dSL','vlcCfg_rx_dFL','vlcCfg_rx_dPX','vlcCfg_rx_thetaLut',...
                                  'vlcCfg_rx_lum2pkpk','vlcCfg_rx_thetaSwpRes','vlcCfg_rx_dRes','vlcCfg_tx_maxLum',...
                                  'vlcCfg_tx_plrPatFit','vlcCfg_tx_plrPatNum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Config Load Callback %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cfgLd_button_Callback(hObject, eventdata, handles)
[fn, pt] = uigetfile;
load(strcat(pt,fn));

handles.maxLum = vlcCfg_tx_maxLum;
handles.txPlrPatFit = vlcCfg_tx_plrPatFit;
handles.txPlrPatNum = vlcCfg_tx_plrPatNum;
polarplot(handles.h_TxPlrPatAxis,deg2rad(vlcCfg_tx_plrPatNum(1,:)),vlcCfg_tx_plrPatNum(2,:),'g','LineWidth',2)
ax = handles.h_TxPlrPatAxis;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
ax.ThetaLim = [-90 90];
ax.RTick = [0:0.2:1];
hold(handles.h_TxPlrPatAxis,'on');
handles.h_plr_L = plot(handles.h_TxPlrPatAxis,0,0,'b','LineWidth',1);
handles.h_plr_R = plot(handles.h_TxPlrPatAxis,0,0,'r','LineWidth',1);
hold(handles.h_TxPlrPatAxis,'off');

handles.dDL = vlcCfg_rx_dDL;
handles.dSL = vlcCfg_rx_dSL;
handles.dFL = vlcCfg_rx_dFL;
handles.dPX = vlcCfg_rx_dPX;
handles.lum2pkpk = vlcCfg_rx_lum2pkpk;
% handles.dRes = vlcCfg_rx_dRes;
% handles.thetaSwpRes = vlcCfg_rx_thetaSwpRes;
set(handles.dDL_edit,'String',num2str(handles.dDL,'%2.2e'))
set(handles.dSL_edit,'String',num2str(handles.dSL,'%2.2e'))
set(handles.dFL_edit,'String',num2str(handles.dFL,'%2.2e'))
set(handles.dPX_edit,'String',num2str(handles.dPX,'%2.2e'))
set(handles.dRes_edit,'String',num2str(handles.dRes,'%2.2e'))
set(handles.lum2pkpk_edit,'String',num2str(handles.lum2pkpk))
handles.lat = 2;
handles.fwd = 10;
handles.hdg = 0;
[theta_l, theta_r] = vlcCfgTool_spotQrx(handles);
handles.theta_l=theta_l;
handles.theta_r=theta_r;
vlcCfgTool_sysTest(handles);
xlabel(handles.thetaLutAxis,'Theta (deg)')
ylabel(handles.thetaLutAxis,'Spot Position (p.u.)')
[~, handles.lutng, handles.lutps] = vlcCfgTool_fovLut(handles);
handles.fov_drawn_flag = 1;
hold(handles.thetaLutAxis,'on')
handles.hL = plot(handles.thetaLutAxis,rad2deg([handles.theta_l handles.theta_l]), [-2 2],'b','LineWidth',1);
handles.hR = plot(handles.thetaLutAxis,rad2deg([handles.theta_r handles.theta_r]), [-2 2],'r','LineWidth',1);
hold(handles.thetaLutAxis,'off')
if ((handles.lutng>rad2deg(theta_l))||(handles.lutng>rad2deg(theta_r))||(handles.lutps<rad2deg(theta_l))||(handles.lutps<rad2deg(theta_r)))
    set(handles.qrx1Lum_txt,'String','Out of FoV!');
    set(handles.qrx2Lum_txt,'String','Out of FoV!');
else
    [qrx1Lum, qrx2Lum] = vlcCfgTool_rxGainLum2PkPk(handles);
    set(handles.qrx1Lum_txt,'String',num2str(abs(qrx1Lum*handles.lum2pkpk)));
    set(handles.qrx2Lum_txt,'String',num2str(abs(qrx2Lum*handles.lum2pkpk)));
end
% Update handles structure
guidata(hObject, handles);
