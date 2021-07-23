function varargout = RetrospectiveValveTracker(varargin)
% RETROSPECTIVEVALVETRACKER MATLAB code for RetrospectiveValveTracker.fig
%      RETROSPECTIVEVALVETRACKER, by itself, creates a new RETROSPECTIVEVALVETRACKER or raises the existing
%      singleton*.
%
%      H = RETROSPECTIVEVALVETRACKER returns the handle to a new RETROSPECTIVEVALVETRACKER or the handle to
%      the existing singleton*.
%
%      RETROSPECTIVEVALVETRACKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RETROSPECTIVEVALVETRACKER.M with the given input arguments.
%
%      RETROSPECTIVEVALVETRACKER('Property','Value',...) creates a new RETROSPECTIVEVALVETRACKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RetrospectiveValveTracker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RetrospectiveValveTracker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RetrospectiveValveTracker

% Last Modified by GUIDE v2.5 03-Dec-2020 11:21:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RetrospectiveValveTracker_OpeningFcn, ...
                   'gui_OutputFcn',  @RetrospectiveValveTracker_OutputFcn, ...
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

% --- Executes just before RetrospectiveValveTracker is made visible.
function RetrospectiveValveTracker_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RetrospectiveValveTracker (see VARARGIN)

% Choose default command line output for RetrospectiveValveTracker
handles.output = hObject;
if isempty(varargin{1})
    varargin = varargin(2:end);
    numArg = nargin-4;
else
    numArg = nargin-3;
end
fprintf('Starting up Retrospective Valve Tracker.\n');
if numArg>=0 && (ischar(varargin{1}) || isstring(varargin{1}))
    if exist(varargin{1},'file') && contains(varargin{1},'.mat')
        oldHandles = load(varargin{1});
        oldHandles=oldHandles.handles;
        for fn = fieldnames(oldHandles)'
           if ~isfield(handles,fn) || ~isgraphics(handles.(fn{1}))
               handles.(fn{1})= oldHandles.(fn{1});
           end
        end
    elseif exist(varargin{1},'dir')
        handles.pcvipr= readPCVIPR(varargin{1});
        if numArg>=5 && isnumeric(varargin{5})
            handles.pcvipr.mag = circshift(handles.pcvipr.mag,[0,0,0,varargin{5}]);
            handles.pcvipr.cd = circshift(handles.pcvipr.cd,[0,0,0,varargin{5}]);
            handles.pcvipr.velX = circshift(handles.pcvipr.velX,[0,0,0,varargin{5}]);
            handles.pcvipr.velY = circshift(handles.pcvipr.velY,[0,0,0,varargin{5}]);
            handles.pcvipr.velZ = circshift(handles.pcvipr.velZ,[0,0,0,varargin{5}]);
        end
        handles.extraShift = zeros(3,1);
        if numArg>=2 && (ischar(varargin{2}) || isstring(varargin{2}))
            [handles.twoCh,handles.threeCh,handles.fourCh] = readLongAxisDir(varargin{2});
        end
        if numArg>=3 && (ischar(varargin{3}) || isstring(varargin{3})) && exist(varargin{3},'dir')
            handles.lvot = readDicomList(dir(fullfile(varargin{3},'*.dcm')));
        end
        if numArg>=4 && (ischar(varargin{4}) || isstring(varargin{4})) && exist(varargin{4},'dir')
            handles.rvot = readDicomList(dir(fullfile(varargin{4},'*.dcm')));
        end
        if numArg>=6 && (ischar(varargin{6}) || isstring(varargin{6}))
            handles.saveFile = varargin{6};
        else
            handles.saveFile = '';
        end
        handles.gridSize = [128, 64, 256, 256];
        handles.imgSize = nan(4,2);
        handles.zoomFactor = ones(4,4);
        handles.center = nan(4,4,2);
        handles.im1 = mean(handles.pcvipr.cd,4);
        handles.center(:,1,:) = repmat(cat(3,size(handles.im1,1)/2, size(handles.im1,2)/2),4,1,1);
        handles.imgSize(1,:) = [size(handles.im1,1), size(handles.im1,2)];
        handles.time=1;
        handles.slice=repmat(round(size(handles.im1,3)/2),4);
        handles.slider2.Value = handles.slice(1);
        handles.slider2.Max = size(handles.im1,3);
        handles.valveBoundaries = nan(4,4,20,2,2);
        handles.planeValveBoundaryPoints = nan(4,20,10,2);
        handles.planeValveBoundaryCurve = nan(4,20,30,2);
        handles.p = nan(4,4,3);
        handles.R = nan(4,4,3,3);
        handles.p(:,1,:) = reshape(repmat(handles.pcvipr.headerPos',4,1),4,1,3);
        R = repmat(handles.pcvipr.spacing',3,1).*reshape(handles.pcvipr.orientation,3,3);
        handles.R(1,1,:,:) = reshape(R,1,1,3,3);
        handles.R(2,1,:,:) = reshape(R,1,1,3,3);
        handles.R(3,1,:,:) = reshape(R,1,1,3,3);
        handles.R(4,1,:,:) = reshape(R,1,1,3,3);
        handles.valveDirection=ones(4,1);
        handles.numViews = [2;1;2;1];
        if ~isfield(handles,'lvot'); handles.numViews(3)=1; end
        handles.windowWidth = ones(4,4).*100;
        handles.windowWidth(:,2) = 25;
        handles.windowLevel = ones(4,4).*50;
        handles.windowLevel(:,2) = 0;
        handles = setActiveValve(handles,1);
        set(handles.axes1, 'ButtonDownFcn', {@axesClickdown,1});
        set(handles.axes2, 'ButtonDownFcn', {@axesClickdown,2});
        set(handles.axes3, 'ButtonDownFcn', {@axesClickdown,3});
        set(handles.axes4, 'ButtonDownFcn', {@axesClickdown,4});
        set(handles.figure1,'windowscrollWheelFcn', @FigureScrollCallback);
        set(handles.figure1,'KeyPressFcn', @FigureKeyCallback);
    end
end
% Update handles structure
guidata(hObject, handles);
handles = updateUI(handles);
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = RetrospectiveValveTracker_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
if ~isempty(handles) % Done button pressed
    varargout{1} = computeFlow(handles,1); % Mitral Valve
    varargout{2} = computeFlow(handles,2); % Tricuspid Valve
    varargout{3} = computeFlow(handles,3); % Aortic Valve
    varargout{4} = computeFlow(handles,4); % Pulmonary Valve
    for fn = fieldnames(handles)'
        if isgraphics(handles.(fn{1}))
            handles = rmfield(handles,fn{1});
        end
    end
    save(handles.saveFile,'handles', '-v7.3');
    closereq();
else % Figure x'd out.
    varargout{1} = [];
    varargout{2} = []; 
    varargout{3} = [];
    varargout{4} = [];
end

function valveFlow = computeFlow(handles,valveNum)
if ~all(isnan(handles.planeValveBoundaryCurve(valveNum,:,:,1))) 
    nT = handles.pcvipr.nT;
    bX = squeeze(handles.planeValveBoundaryCurve(valveNum,:,:,1))+handles.gridSize(2)/2+1; % 20 x 30
    bY = squeeze(handles.planeValveBoundaryCurve(valveNum,:,:,2))+handles.gridSize(2)/2+1;
    vel = zeros(handles.gridSize(2),handles.gridSize(2),nT);
    flow = zeros(handles.pcvipr.nT,1);
    forwardFlow = zeros(handles.pcvipr.nT,1);
    reverseFlow = zeros(handles.pcvipr.nT,1);
    [x,y] = ndgrid(1:handles.gridSize(2),1:handles.gridSize(2));
    for tt = 1:nT
        [vel(:,:,tt),~, ~, ~] = getVel(handles,valveNum,tt);
        mask = inpolygon(x,y,bX(tt,:),bY(tt,:));
        maskF = mask.*(vel(:,:,tt)>=0);
        maskR = mask.*(vel(:,:,tt)<0);
        
        handles = getValvePlane(handles,handles.numViews(valveNum), tt);
        dx=handles.size(valveNum)*handles.pcvipr.spacing(1)/handles.gridSize(2);
        dy=handles.size(valveNum)*handles.pcvipr.spacing(2)/handles.gridSize(2);
        flow(tt) = sum(sum(vel(:,:,tt).*mask))*dx*dy*60/1e6;
        forwardFlow(tt) = sum(sum(vel(:,:,tt).*maskF))*dx*dy*60/1e6;
        reverseFlow(tt) = sum(sum(vel(:,:,tt).*maskR))*dx*dy*60/1e6;
    end
    valveFlow.netFlow = mean(flow);
    valveFlow.regurgFxn = abs(mean(reverseFlow))/mean(flow);
    valveFlow.peakFlow = max(flow);
    indE = round(nT/3):(nT-5);
    indA = (nT-4):nT;
    [valveFlow.eFlow,eInd] = max(flow(indE));
    eInd = indE(eInd);
    valveFlow.aFlow = max(flow(indA));
    valveFlow.eaRatio = max(flow(indE))/max(flow(indA));
    dfdt = (flow(eInd+2)-flow(eInd+1))/handles.pcvipr.dT;
    valveFlow.dt = -valveFlow.eFlow/dfdt;
else
    valveFlow = [];
end
    
function handles = setActiveValve(handles, val)
handles.activeValve = val;
switch val
    case 1 % mitral valve
        handles.im3 = handles.fourCh.img;
        handles.im4 = handles.twoCh.img;
        handles.p(val,3,:) = handles.fourCh.p;
        handles.p(val,4,:) = handles.twoCh.p;
        handles.R(val,3,:,:) = repmat(handles.fourCh.spacing,3,1).*handles.fourCh.R;
        handles.R(val,4,:,:) = repmat(handles.twoCh.spacing,3,1).*handles.twoCh.R;
    case 2 % tricuspid valve
        handles.im3 = handles.fourCh.img;
        handles.im4 = [];
        handles.p(val,3,:) = handles.fourCh.p;
        handles.R(val,3,:,:) = repmat(handles.fourCh.spacing,3,1).*handles.fourCh.R;
    case 3 % aortic valve
        handles.im3 = handles.threeCh.img;
        handles.p(val,3,:) = handles.threeCh.p;
        handles.R(val,3,:,:) = repmat(handles.threeCh.spacing,3,1).*handles.threeCh.R;
        if isfield(handles, 'lvot'); handles.im4 = handles.lvot.img; 
            handles.p(val,4,:) = handles.lvot.p;
            handles.R(val,4,:,:) = repmat(handles.lvot.spacing,3,1).*handles.lvot.R;
        else
            handles.im4 = [];
        end
    case 4 % pulmonary valve
        handles.im3 = handles.rvot.img;
        handles.im4 = [];
        handles.p(val,3,1:3) = handles.rvot.p;
        handles.R(val,3,:,:) = repmat(handles.rvot.spacing,3,1).*handles.rvot.R;
end
if isnan(handles.center(handles.activeValve,3,1))
    handles.center(handles.activeValve,3,:) = [size(handles.im3,1)/2, size(handles.im3,2)/2];
end
if isnan(handles.center(handles.activeValve,4,1))
    handles.center(handles.activeValve,4,:) = [size(handles.im4,1)/2, size(handles.im4,2)/2];
end
handles.imgSize(3:4,:) = [size(handles.im3,1), size(handles.im3,2);...
                          size(handles.im4,1), size(handles.im4,2)];
handles.slider2.Value = handles.slice(val);                     
handles = updateUI(handles);

function FigureKeyCallback(hObject, eventdata)
handles = guidata(hObject);
add=0;
if strcmp(eventdata.Key,'rightarrow') && round(handles.slider1.Value)<handles.slider1.Max
    add=1;
elseif strcmp(eventdata.Key,'leftarrow') && round(handles.slider1.Value)>handles.slider1.Min
    add=-1;
end
if add~=0
    handles.slider1.Value = round(handles.slider1.Value+add);
    handles.time = round(handles.slider1.Value);
    handles = updateUI(handles);
    guidata(hObject,handles);
end

function [twoCh, threeCh, fourCh] = readLongAxisDir(dirName)
dirName = string(dirName);
if numel(dirName)==1
    dicoms = dir(fullfile(dirName,'*.dcm'));
    set1 = readDicomList(dicoms(1:(numel(dicoms)/3)));
    set2 = readDicomList(dicoms((numel(dicoms)/3+1):(2*numel(dicoms)/3)));
    set3 = readDicomList(dicoms((2*numel(dicoms)/3+1):end));
    Rs = [set1.R(:,3),set2.R(:,3),set3.R(:,3)];
    [~,ind2Ch] = min(abs(Rs(3,:)));
    sets = [set1,set2,set3];
    twoCh = sets(ind2Ch);
    sets = sets(setdiff(1:3,ind2Ch));
    Rs = Rs(:,setdiff(1:3,ind2Ch));
    [~,ind4Ch] = min(abs(Rs(1,:)));
    fourCh = sets(ind4Ch);
    threeCh = sets(setdiff(1:2,ind4Ch));
else
    twoCh = readDicomList(dir(fullfile(dirName(1),'*.dcm')));
    threeCh = readDicomList(dir(fullfile(dirName(2),'*.dcm')));
    fourCh = readDicomList(dir(fullfile(dirName(3),'*.dcm')));
end

function img = readDicomList(dicoms)
for ii = 1:numel(dicoms)
    dicomFilePath = fullfile(dicoms(ii).folder, dicoms(ii).name);
    if ii==1
        info = dicominfo(dicomFilePath);
        img.p = double(info.ImagePositionPatient);
        R = double(info.ImageOrientationPatient);
        img.spacing = double([info.PixelSpacing',info.SliceThickness]);
    end
    img.img(:,:,ii) = single(dicomread(dicomFilePath))';
end
R(7:9) = cross(R(1:3),R(4:6));
img.R = reshape(R, [3,3]);

function axesClickdown(hObject, eventdata, axisNum)
handles = guidata(hObject);
if strcmp(handles.figure1.SelectionType,'normal') && axisNum~=2
    handles.figure1.WindowButtonMotionFcn = {@pan,axisNum,eventdata.IntersectionPoint./handles.gridSize(axisNum)};
    handles.figure1.WindowButtonUpFcn = {@buttonUpCallback,axisNum,eventdata.IntersectionPoint./handles.gridSize(axisNum)};
elseif strcmp(handles.figure1.SelectionType,'normal') && axisNum==2
    handles.figure1.WindowButtonMotionFcn = @endPan;
    handles.figure1.WindowButtonUpFcn = {@buttonUpAxis2,eventdata.IntersectionPoint./handles.gridSize(axisNum)};
elseif strcmp(handles.figure1.SelectionType,'alt')
    handles.figure1.WindowButtonMotionFcn = {@windowLevel,axisNum,eventdata.IntersectionPoint./handles.gridSize(axisNum)};
    handles.figure1.WindowButtonUpFcn = @endPan;
end
guidata(hObject, handles);

function lineClickdown(hObject, ~, axisNum)
handles = guidata(hObject);
handles.figure1.WindowButtonMotionFcn = {@moveLine,axisNum};
handles.figure1.WindowButtonUpFcn = '';
guidata(hObject, handles);

function pan(hObject,eventdata,axisNum,startPoint) %#ok<INUSL>
handles = guidata(hObject);
handles.figure1.WindowButtonUpFcn = @endPan;
x = handles.figure1.CurrentPoint(1);
y = handles.figure1.CurrentPoint(2);
allAxes = [handles.axes1, handles.axes2, handles.axes3, handles.axes4];
axes = allAxes(axisNum);
x = (x/handles.figure1.InnerPosition(3) - axes.Position(1))/axes.Position(3);
dx = x-startPoint(1);
y = (y/handles.figure1.InnerPosition(4) - axes.Position(2))/axes.Position(4);
dy = y-(1-startPoint(2));
handles.center(handles.activeValve,axisNum,:) = handles.center(handles.activeValve,axisNum,:) + reshape(handles.gridSize(axisNum).*[dy;-dx],1,1,2);
handles = updateUI(handles);
handles.figure1.WindowButtonMotionFcn = {@pan,axisNum,[x,1-y]};
guidata(hObject, handles);

function windowLevel(hObject,eventdata,axisNum,startPoint) %#ok<INUSL>
handles = guidata(hObject);
handles.figure1.WindowButtonUpFcn = @endPan;
x = handles.figure1.CurrentPoint(1);
y = handles.figure1.CurrentPoint(2);
allAxes = [handles.axes1, handles.axes2, handles.axes3, handles.axes4];
axes = allAxes(axisNum);
x = (x/handles.figure1.InnerPosition(3) - axes.Position(1))/axes.Position(3);
dx = x-startPoint(1);
y = (y/handles.figure1.InnerPosition(4) - axes.Position(2))/axes.Position(4);
dy = y-(1-startPoint(2));
if axisNum==2
    wl=0;
    ww = max(min(handles.windowWidth(handles.activeValve,axisNum) + dx*5, 100),0);
else
    wl = max(min(handles.windowLevel(handles.activeValve,axisNum) + dy*5,100),0);
    ww = max(min(handles.windowWidth(handles.activeValve,axisNum) + dx*5, min(100-wl,wl)*2),0);
end
handles.windowLevel(handles.activeValve,axisNum) = wl;
handles.windowWidth(handles.activeValve,axisNum) = ww;
handles = updateUI(handles);
handles.figure1.WindowButtonMotionFcn = {@windowLevel,axisNum,[x,1-y]};
guidata(hObject, handles);

function moveLine(hObject,~,axisNum) 
handles = guidata(hObject);
handles.figure1.WindowButtonUpFcn = @endPan;
x = handles.figure1.CurrentPoint(1);
y = handles.figure1.CurrentPoint(2);
allAxes = [handles.axes1, handles.axes2, handles.axes3, handles.axes4];
axes = allAxes(axisNum);
x = (x/handles.figure1.InnerPosition(3) - axes.Position(1))/axes.Position(3);
y = (y/handles.figure1.InnerPosition(4) - axes.Position(2))/axes.Position(4);
if axisNum==2
    x2 = (x-0.5)*handles.gridSize(2);
    y2 = (0.5-y)*handles.gridSize(2);
    x0 = squeeze(handles.planeValveBoundaryPoints(handles.activeValve, handles.time, :, 1));
    y0 = squeeze(handles.planeValveBoundaryPoints(handles.activeValve, handles.time, :, 2));
    [~,closerPoint] = min((repmat(y2,sum(~isnan(x0)),1)-x0(~isnan(x0))).^2 + (repmat(x2,sum(~isnan(y0)),1)-y0(~isnan(y0))).^2);
    handles.planeValveBoundaryPoints(handles.activeValve, handles.time, closerPoint,:) = [y2,x2];
else
    x2 = (x-0.5)*handles.imgSize(axisNum,1)/handles.zoomFactor(handles.activeValve,axisNum)+handles.center(handles.activeValve,axisNum,2);
    y2 = (0.5-y)*handles.imgSize(axisNum,2)/handles.zoomFactor(handles.activeValve,axisNum)+handles.center(handles.activeValve,axisNum,1);
    x0 = squeeze(handles.valveBoundaries(handles.activeValve,axisNum, handles.time, :, 1));
    y0 = squeeze(handles.valveBoundaries(handles.activeValve,axisNum, handles.time, :, 2));
    [~,closerPoint] = min(([y2;y2]-x0).^2 + ([x2;x2]-y0).^2);
    handles.valveBoundaries(handles.activeValve, axisNum, handles.time, closerPoint,:) = [y2,x2];
end
handles = updateUI(handles);
handles.figure1.WindowButtonMotionFcn = {@moveLine,axisNum};
guidata(hObject, handles);

function endPan(hObject, eventdata) %#ok<INUSD>
handles = guidata(hObject);
handles.figure1.WindowButtonMotionFcn = '';
handles.figure1.WindowButtonUpFcn = '';
guidata(hObject, handles);

function buttonUpCallback(hObject, eventdata, axisNum, startPoint)  %#ok<INUSL>
handles = guidata(hObject);
handles.figure1.WindowButtonMotionFcn = '';
handles.figure1.WindowButtonUpFcn = '';

x = handles.center(handles.activeValve,axisNum,1) + handles.imgSize(axisNum,1)/handles.zoomFactor(handles.activeValve,axisNum)*(startPoint(2)-0.5)-1;
y = handles.center(handles.activeValve,axisNum,2) + handles.imgSize(axisNum,2)/handles.zoomFactor(handles.activeValve,axisNum)*(startPoint(1)-0.5)-1;
if isnan(handles.valveBoundaries(handles.activeValve,axisNum,handles.time,2,1))
    if isnan(handles.valveBoundaries(handles.activeValve,axisNum,handles.time,1,1))
        handles.valveBoundaries(handles.activeValve,axisNum,handles.time,1,1) = x;
        handles.valveBoundaries(handles.activeValve,axisNum,handles.time,1,2) = y;
    else
        handles.valveBoundaries(handles.activeValve,axisNum,handles.time,2,1) = x;
        handles.valveBoundaries(handles.activeValve,axisNum,handles.time,2,2) = y;
    end 
    handles = updateUI(handles);
end
% fprintf('Click detected on axis %i at (%.2f,%.2f)\n', axisNum, x, y);
guidata(hObject, handles);
function buttonUpAxis2(hObject, eventdata, startPoint)  %#ok<INUSL>
handles = guidata(hObject);
handles.figure1.WindowButtonMotionFcn = '';
handles.figure1.WindowButtonUpFcn = '';

x = handles.gridSize(2)*(startPoint(2)-0.5)-1;
y = handles.gridSize(2)*(startPoint(1)-0.5)-1;
bps = squeeze(handles.planeValveBoundaryPoints(handles.activeValve,handles.time,:,:));
indNan = find(isnan(bps(:,1)),1);
if ~isempty(indNan)
    bps(indNan,1) = x;
    bps(indNan,2) = y;
    handles.planeValveBoundaryPoints(handles.activeValve,handles.time,:,:) = reshape(bps,size(handles.planeValveBoundaryPoints(handles.activeValve,handles.time,:,:)));
    handles = updateUI(handles);
end
guidata(hObject, handles);

function handles = updateUI(handles)
updateAxis(handles, handles.axes1, handles.im1, handles.slice(handles.activeValve), 1);
updateAxis(handles, handles.axes3, handles.im3, handles.time, 3);
updateAxis(handles, handles.axes4, handles.im4, handles.time, 4);
if shouldPlotValvePlane(handles)
    handles = getValvePlane(handles,shouldPlotValvePlane(handles), handles.time);
    handles = plotValvePlane(handles);
else
   delete(handles.axes2.Children);
end

function [vel, midP, sz, norm] = getVel(handles, valveNum, time)
midP = squeeze(handles.midpoint(valveNum,:));

sz = handles.size(valveNum);
norm = squeeze(handles.normal(valveNum,:))';
gridSize = squeeze(handles.gridSize(2));

d1 = [1;0;0];
d1 = cross(d1,norm);
d2 = cross(d1,norm);
d1 = d1./sqrt(d1'*d1);
d2 = d2./sqrt(d2'*d2);
if sum((cross(d1,d2)-norm).^2)>0.1; tmp=d1; d1=d2; d2=tmp; end

[x,y,z] = ndgrid(1:handles.imgSize(1,1),1:handles.imgSize(1,2),1:handles.slider2.Max);
[xx,yy] = ndgrid(linspace(-sz/2,sz/2,gridSize),linspace(-sz/2,sz/2,gridSize));
xq = d1(1)*xx + d2(1)*yy + midP(1);
yq = d1(2)*xx + d2(2)*yy + midP(2);
zq = d1(3)*xx + d2(3)*yy + midP(3);

vx = interpn(x,y,z,single(handles.pcvipr.velX(:,:,:,time)),xq,yq,zq);
vy = interpn(x,y,z,single(handles.pcvipr.velY(:,:,:,time)),xq,yq,zq);
vz = interpn(x,y,z,single(handles.pcvipr.velZ(:,:,:,time)),xq,yq,zq);
vel = -vx*norm(1) - vy*norm(2) - vz*norm(3);

function handles = plotValvePlane(handles)
axes1 = handles.axes1;
axes2 = handles.axes2;


[img, midP, sz, norm] = getVel(handles, handles.activeValve, handles.time);
range = 2*handles.pcvipr.venc;
lowerLim = -handles.windowWidth(handles.activeValve, 2)/100*range/2;
upperLim = handles.windowWidth(handles.activeValve, 2)/100*range/2;
% fprintf('Window limits = (%f, %f)\n', lowerLim, upperLim);
img = max(min(img,upperLim),lowerLim);
handles.figure1.CurrentAxes = axes2;
if isempty(axes2.Children)
    h=image(img, 'CDataMapping','scaled', 'Parent', axes2);
    set(axes2,  'HitTest', 'on', 'ButtonDownFcn', {@axesClickdown,2}, 'XTick', [], 'YTick', [], 'CLim', [lowerLim, upperLim]);
    set(h, 'HitTest', 'off');
    colorbar();
    colormap(gray);
else
    if numel(axes2.Children)>1; delete(axes2.Children(1:(end-1))); end
    axes2.Children(end).CData = img;
    axes2.CLim = [lowerLim, upperLim];
end
hold on;

x = squeeze(handles.planeValveBoundaryPoints(handles.activeValve,handles.time,:,1));
x2 = (x(~isnan(x))/handles.gridSize(2)+0.5)*handles.gridSize(2)+1;
y = squeeze(handles.planeValveBoundaryPoints(handles.activeValve,handles.time,:,2));
y2 = (y(~isnan(y))/handles.gridSize(2)+0.5)*handles.gridSize(2)+1;
if ~isempty(x2)
    plot(y2,x2,'or','ButtonDownFcn',{@lineClickdown,2}, 'MarkerSize', 5, 'MarkerFaceColor','r');
    if numel(x2)>=5
        oldT = 1:(numel(x2)+3);
        newT = linspace(2,numel(x2)+2,30);
        x3 = csapi(oldT,[x2;x2(1:3)],newT);
        y3 = csapi(oldT,[y2;y2(1:3)],newT);
        handles.planeValveBoundaryCurve(handles.activeValve,handles.time,:,1) = ((reshape(x3,1,1,30,1)-1)/handles.gridSize(2) - 0.5)*handles.gridSize(2);
        handles.planeValveBoundaryCurve(handles.activeValve,handles.time,:,2) = ((reshape(y3,1,1,30,1)-1)/handles.gridSize(2) - 0.5)*handles.gridSize(2);
        plot(y3,x3,'-r','ButtonDownFcn',{@lineClickdown,2});
    end
end
hold off;

        
z = handles.slice(handles.activeValve);
zRange = [midP(3)-sz*(1-abs(norm(3)))/2, midP(3)+sz*(1-abs(norm(3)))/2];
if z>zRange(1) && z<zRange(2)
    n = norm(1:2)*sz;
    [x0,y0,m] = getIntersection(midP(1), midP(2), midP(3), norm(1), norm(2), norm(3), z, sz/2);
    handles.figure1.CurrentAxes = axes1;
    hold on;
    x = [m(1),m(1)+n(1),x0];
    x2 = ((x-handles.center(handles.activeValve,1,1)+1)*handles.zoomFactor(handles.activeValve,1)/handles.imgSize(1,1)+0.5)*handles.gridSize(1);
    y = [m(2),m(2)+n(2),y0];
    y2 = ((y-handles.center(handles.activeValve,1,2)+1)*handles.zoomFactor(handles.activeValve,1)/handles.imgSize(1,2)+0.5)*handles.gridSize(1);
    plot(y2(3:4),x2(3:4),'-r');
    quiver(y2(1),x2(1),y2(2)-y2(1),x2(2)-x2(1),'r','AutoScaleFactor', 0.5, 'MaxHeadSize',2);
    hold off;
end

function [x,y,m] = getIntersection(x0,y0,z0,x1,y1,z1,z,r)
dSlice = abs(r*sqrt(1 - z1^2));
dx = r*y1*sqrt(1-abs(z-z0)/dSlice);
dy = -r*x1*sqrt(1-abs(z-z0)/dSlice);
x = [x0-dx, x0+dx];
y = [y0-dy, y0+dy];
m = mean([x',y']);

function [ijk2,xyz] = convertCoords(ijk,p,R,p2,R2)
xyz = R*(ijk-ones(size(ijk))) + repmat(p,1,size(ijk,2));
ijk2 = pinv([[R2,p2];[0,0,0,1]])*[xyz;ones(1,size(xyz,2))];
ijk2 = ijk2(1:(end-1),:)+ones(size(ijk2(1:(end-1),:)));

function [valve1IJKvipr, d, R] = getValveEnds(handles,valveNum,axisNum,time)
p = squeeze(handles.p(valveNum,axisNum,:));
R = squeeze(handles.R(valveNum,axisNum,:,:));
valve1 = squeeze(handles.valveBoundaries(valveNum,axisNum,time,:,:));
p2 = squeeze(handles.p(valveNum,1,:));
R2 = squeeze(handles.R(valveNum,1,:,:));
[valve1IJKvipr,~]=convertCoords([valve1';1,1],p,R,p2+handles.extraShift,R2);
valve1IJKvipr = valve1IJKvipr';
d = (valve1IJKvipr(2,:)-valve1IJKvipr(1,:))';

function handles = getValvePlane(handles,numViews, time)
[valve1IJKvipr,d1,R1] =  getValveEnds(handles,handles.activeValve,3,time);
if numViews==2
    [valve2IJKvipr,d2,~] =  getValveEnds(handles,handles.activeValve,4,time);
    handles.midpoint(handles.activeValve,:) = mean([valve1IJKvipr;valve2IJKvipr],1);
    handles.size(handles.activeValve) = 1.5*sqrt(max(d1'*d1,d2'*d2));
else
    handles.midpoint(handles.activeValve,:) = mean(valve1IJKvipr,1);
    [ijk2,~] = convertCoords([0,0,0;0,0,1]',squeeze(handles.p(handles.activeValve,3,:)),R1, squeeze(handles.p(handles.activeValve,1,:)),squeeze(handles.R(handles.activeValve,1,:,:)));
    d2 = (ijk2(:,2)-ijk2(:,1));
    handles.size(handles.activeValve) = 1.5*sqrt(d1'*d1);
end
norm = cross(d1./sqrt(d1'*d1),d2./sqrt(d2'*d2)).*handles.valveDirection(handles.activeValve);
handles.normal(handles.activeValve,:) = norm./sqrt(norm'*norm);

function val = shouldPlotValvePlane(handles)
val=0;
if ~isnan(handles.valveBoundaries(handles.activeValve,3,handles.time,2,1)) % axis 3 has both points place
    val=1;
    if (handles.numViews(handles.activeValve)==2) && ~isnan(handles.valveBoundaries(handles.activeValve,4,handles.time,2,1))
        val=2;
    end
end

function updateAxis(handles, axes, img, ind, axisNum)
if ~isempty(img)
    lowerLim = handles.windowLevel(handles.activeValve, axisNum)-handles.windowWidth(handles.activeValve, axisNum)/2;
    upperLim = handles.windowLevel(handles.activeValve, axisNum)+handles.windowWidth(handles.activeValve, axisNum)/2;
    img2 = zoomImg(img(:,:,ind),handles.center(handles.activeValve, axisNum,:),handles.zoomFactor(handles.activeValve,axisNum), handles.gridSize(axisNum),lowerLim,upperLim);
    if isempty(axes.Children)
        h=image(img2, 'CDataMapping','scaled', 'Parent', axes);
        set(axes, 'ButtonDownFcn', {@axesClickdown,axisNum}, 'HitTest', 'on', 'XTick', [], 'YTick', []);
        set(h, 'HitTest', 'off');
        colormap(gray);
    else
        if numel(axes.Children)>1; delete(axes.Children(1:(end-1))); end
        axes.Children(end).CData = img2;
    end
    handles.figure1.CurrentAxes = axes;
    hold on;
    x = squeeze(handles.valveBoundaries(handles.activeValve,axisNum,handles.time,:,1));
    x2 = ((x-handles.center(handles.activeValve,axisNum,1)+1)*handles.zoomFactor(handles.activeValve,axisNum)/handles.imgSize(axisNum,1)+0.5)*handles.gridSize(axisNum);
    y = squeeze(handles.valveBoundaries(handles.activeValve,axisNum,handles.time,:,2));
    y2 = ((y-handles.center(handles.activeValve,axisNum,2)+1)*handles.zoomFactor(handles.activeValve,axisNum)/handles.imgSize(axisNum,2)+0.5)*handles.gridSize(axisNum);
    if ~isnan(x(1))
        plot(y2(1),x2(1),'og','ButtonDownFcn',{@lineClickdown,axisNum}, 'MarkerSize', 5, 'MarkerFaceColor','g');
        plot(y2(2),x2(2),'oy','ButtonDownFcn',{@lineClickdown,axisNum}, 'MarkerSize', 5, 'MarkerFaceColor','y');
        plot(y2,x2,'-r','ButtonDownFcn',{@lineClickdown,axisNum}, 'LineWidth', 2);
    end
    hold off;
else
    if ~isempty(axes.Children)
        if numel(axes.Children)>1; delete(axes.Children(1:(end-1))); end
        axes.Children(end).CData = zeros(handles.gridSize(axisNum));
    end 
end

function zoomedImg = zoomImg(img,center,zoomFactor,gridSize,lowerLim,upperLim)
[x,y] = ndgrid(1:size(img,1),1:size(img,2));
[xq,yq] = ndgrid(((-gridSize/2):(gridSize/2))/(zoomFactor*gridSize/size(img,1))+center(1),...
                 ((-gridSize/2):(gridSize/2))/(zoomFactor*gridSize/size(img,2))+center(2));
zoomedImg = max(min(interp2(y,x,img,yq,xq),prctile(img(:),upperLim)),prctile(img(:),lowerLim));

function FigureScrollCallback(hObject,eventdata)
if eventdata.VerticalScrollCount~=0 
    dy = -eventdata.VerticalScrollCount;
    handles = guidata(hObject);

    x = handles.figure1.CurrentPoint(1);
    y = handles.figure1.CurrentPoint(2);
    if x>0 && x<handles.figure1.InnerPosition(3) && y>0 && y<handles.figure1.InnerPosition(4)
        [scrollAxis, xAx, yAx] = whichAxis(handles,x/handles.figure1.InnerPosition(3),y/handles.figure1.InnerPosition(4));
        handles = zoomAxis(handles,scrollAxis, xAx, yAx, dy);
        handles = updateUI(handles);
    end
    guidata(hObject,handles);
end

function handles = zoomAxis(handles,scrollAxis, xAx, yAx, dy)
handles.zoomFactor(handles.activeValve,scrollAxis) = handles.zoomFactor(handles.activeValve,scrollAxis)*(1+dy/10);
handles.center(handles.activeValve,scrollAxis,:) = handles.center(handles.activeValve,scrollAxis,:) + ...
    reshape(0.01*[(xAx-0.5),(yAx-0.5)].*[handles.imgSize(scrollAxis,1),handles.imgSize(scrollAxis,2)]/handles.zoomFactor(handles.activeValve,scrollAxis),1,1,2);

function [scrollAxis, xAx, yAx] = whichAxis(handles,x,y)
scrollAxis = 0;
xAx = 0.5;
yAx = 0.5;
allAxes = [handles.axes1, handles.axes2, handles.axes3, handles.axes4];
for ii =1:numel(allAxes)
    minX = allAxes(ii).Position(1);
    minY = allAxes(ii).Position(2);
    maxX = allAxes(ii).Position(3) + allAxes(ii).Position(1);
    maxY = allAxes(ii).Position(4) + allAxes(ii).Position(2);
    if x>minX && x<maxX && y>minY && y<maxY
        scrollAxis = ii;
        yAx = (x-minX)/(maxX-minX);
        xAx = 1-(y-minY)/(maxY-minY);
    end
end

function handles = shiftPlane(handles,direction)
norm = squeeze(handles.normal(handles.activeValve,:))';
d1 = [1;0;0];
d1 = cross(d1,norm);
d2 = cross(d1,norm);
d1 = d1./sqrt(d1'*d1);
d2 = d2./sqrt(d2'*d2);
if sum((cross(d1,d2)-norm).^2)>0.1; tmp=d1; d1=d2; d2=tmp; end
handles.extraShift = handles.extraShift + direction(1)*squeeze(handles.R(1,1,:,:))*d1 + direction(2)*squeeze(handles.R(1,1,:,:))*d2;
handles = updateUI(handles);

function handles = interpolateAxis(handles,axes)
for ii = 1:numel(axes)
    ax = axes(ii);
    if ax==2
        [x,y,~] = ndgrid(1:handles.gridSize(2),1:handles.gridSize(2),1:handles.pcvipr.nT);
        goodT = find(~isnan(handles.planeValveBoundaryPoints(handles.activeValve,:,5,1)));
        if ~isempty(goodT)
            badT = setdiff(1:handles.pcvipr.nT,goodT);
            img = zeros(size(x));
            for iii=1:numel(goodT)
                xv = squeeze(handles.planeValveBoundaryCurve(handles.activeValve,goodT(iii),:,1)) + handles.gridSize(2)/2 + 1;
                yv = squeeze(handles.planeValveBoundaryCurve(handles.activeValve,goodT(iii),:,2)) + handles.gridSize(2)/2 + 1;
                img(:,:,goodT(iii)) = inpolygon(x(:,:,goodT(iii)), y(:,:,goodT(iii)), xv, yv);
            end
            for iii=1:numel(goodT)
                if goodT(iii) ~= max(goodT)
                    sl1 = goodT(iii);
                    sl2 = goodT(iii+1);
                    N = sl2-sl1-1;
                    img(:,:,sl1:sl2) = morph_binary(img(:,:,sl1),img(:,:,sl2),N);
                else
                    sl1 = goodT(iii);
                    sl2 = goodT(1);
                    r1 = sl1:handles.pcvipr.nT; 
                    r2 = 1:sl2;
                    img(:,:,[r1,r2]) = morph_binary(img(:,:,sl1),img(:,:,sl2),numel(r1)+numel(r2)-2);
                end
            end
            for iii=1:numel(badT)
                B = bwboundaries(img(:,:,badT(iii)));
                if ~isempty(B) && size(B{1},1)>=10
                    x = B{1}(round((1:8)*((size(B{1},1)+1)/8)-1),1)-handles.gridSize(2)/2-1;
                    y = B{1}(round((1:8)*((size(B{1},1)+1)/8)-1),2)-handles.gridSize(2)/2-1;
                    handles.planeValveBoundaryPoints(handles.activeValve,badT(iii),1:8,1) = reshape(x,1,1,8,1);
                    handles.planeValveBoundaryPoints(handles.activeValve,badT(iii),1:8,2) = reshape(y,1,1,8,1);
                    x = B{1}(round((1:30)*((size(B{1},1)+1)/30)-1),1)-handles.gridSize(2)/2-1;
                    y = B{1}(round((1:30)*((size(B{1},1)+1)/30)-1),2)-handles.gridSize(2)/2-1;
                    handles.planeValveBoundaryCurve(handles.activeValve,badT(iii),:,1) = reshape(x,1,1,30,1);
                    handles.planeValveBoundaryCurve(handles.activeValve,badT(iii),:,2) = reshape(y,1,1,30,1);
                end
            end
        end
    else
        if any(~isnan(handles.valveBoundaries(handles.activeValve,ax,:,1,1)))
            if any(isnan(handles.valveBoundaries(handles.activeValve,ax,:,1,1)))
                handles = interpolatePoint(handles,ax,[1,2]);
            end
        end
    end
end

function handles = interpolatePoint(handles,ax,ps)
for ii=1:numel(ps)
    p = ps(ii);
    N = size(handles.valveBoundaries,3);
    boundsX = repmat(squeeze(handles.valveBoundaries(handles.activeValve,ax,:,p,1)),2,1);
    boundsY = repmat(squeeze(handles.valveBoundaries(handles.activeValve,ax,:,p,2)),2,1);
    ind = find(~isnan(boundsX));
    indNan = find(isnan(boundsX));
    boundsX(indNan) = interp1(ind,boundsX(ind),indNan);
    boundsY(indNan) = interp1(ind,boundsY(ind),indNan);
    handles.valveBoundaries(handles.activeValve,ax,:,p,1) = reshape(boundsX(round([(N+1):(1.5*N),(N/2+1):N])),1,1,N,1,1);
    handles.valveBoundaries(handles.activeValve,ax,:,p,2) = reshape(boundsY(round([(N+1):(1.5*N),(N/2+1):N])),1,1,N,1,1);
end

function radiobutton1_Callback(hObject, ~, handles) %#ok<DEFNU>
if hObject.Value
    handles.radiobutton2.Value = 0;
    handles.radiobutton3.Value = 0;
    handles.radiobutton4.Value = 0;
    handles = setActiveValve(handles, 1);
else
    handles.radiobutton1.Value = 1;
end
guidata(hObject,handles);
function radiobutton2_Callback(hObject, ~, handles) %#ok<DEFNU>
if hObject.Value
    handles.radiobutton1.Value = 0;
    handles.radiobutton3.Value = 0;
    handles.radiobutton4.Value = 0;
    handles = setActiveValve(handles, 2);
else
    handles.radiobutton2.Value = 1;
end
guidata(hObject,handles);
function radiobutton3_Callback(hObject, ~, handles) %#ok<DEFNU>
if hObject.Value
    handles.radiobutton1.Value = 0;
    handles.radiobutton2.Value = 0;
    handles.radiobutton4.Value = 0;
    handles = setActiveValve(handles, 3);
else
    handles.radiobutton3.Value = 1;
end
guidata(hObject,handles);
function radiobutton4_Callback(hObject, ~, handles) %#ok<DEFNU>
if hObject.Value
    handles.radiobutton1.Value = 0;
    handles.radiobutton2.Value = 0;
    handles.radiobutton3.Value = 0;
    handles = setActiveValve(handles, 4);
else
    handles.radiobutton4.Value = 1;
end
guidata(hObject,handles);

function slider1_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.time = round(hObject.Value);
handles = updateUI(handles);
guidata(hObject,handles);
function slider1_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function slider2_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.slice(handles.activeValve) = round(hObject.Value);
handles = updateUI(handles);
guidata(hObject,handles);
function slider2_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function interpolatePushButton_Callback(hObject, ~, handles) %#ok<DEFNU>
handles = interpolateAxis(handles,[2,3,4]);
handles = updateUI(handles);
guidata(hObject,handles);

function doneButton_Callback(hObject, ~, handles) %#ok<DEFNU>
guidata(hObject,handles);
uiresume(handles.figure1);

function pushbutton3_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.extraShift = handles.extraShift + [0;1;0];
handles = updateUI(handles);
guidata(hObject,handles);
function pushbutton4_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.extraShift = handles.extraShift - [0;1;0];
handles = updateUI(handles);
guidata(hObject,handles);
function pushbutton5_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.extraShift = handles.extraShift + [1;0;0];
handles = updateUI(handles);
guidata(hObject,handles);
function pushbutton6_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.extraShift = handles.extraShift - [1;0;0];
handles = updateUI(handles);
guidata(hObject,handles);
function pushbutton7_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.extraShift = handles.extraShift + [0;0;1];
handles = updateUI(handles);
guidata(hObject,handles);
function pushbutton8_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.extraShift = handles.extraShift - [0;0;1];
handles = updateUI(handles);
guidata(hObject,handles);

function pushbutton9_Callback(hObject, ~, handles) %#ok<DEFNU>
handles = shiftPlane(handles,[0,-1]);
guidata(hObject,handles);
function pushbutton10_Callback(hObject, ~, handles) %#ok<DEFNU>
handles = shiftPlane(handles,[0,1]);
guidata(hObject,handles);
function pushbutton11_Callback(hObject, ~, handles) %#ok<DEFNU>
handles = shiftPlane(handles,[-1,0]);
guidata(hObject,handles);
function pushbutton12_Callback(hObject, ~, handles) %#ok<DEFNU>
handles = shiftPlane(handles,[1,0]);
guidata(hObject,handles);
function pushbutton13_Callback(hObject, ~, handles) %#ok<DEFNU>
handles.valveDirection(handles.activeValve) = handles.valveDirection(handles.activeValve)*-1;
handles=updateUI(handles);
guidata(hObject,handles);
