%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2016, Kan Xia                                             %
% This is the Programme that made to predict the process of the threephase%
% production of the air-lift system for solid-liquid-gas.                 %
% The Idea of this Programme comming from the Yoshinaga Literature.       %
% while some parts of it has been adjusted to fit the missing infromation %
% that didn't privided by the Author.                                     %
%                                                                         %
% The author of this code is K.xia from Techonology University of Delft.  %
% Contact Information: ryanllyy.g@gmail.com                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = YoshinagaGUI(varargin)
% YOSHINAGAGUI MATLAB code for YoshinagaGUI.fig
%YOSHINAGAGUI, by itself, creates a new YOSHINAGAGUI or raises the existing
%singleton*.
%
%H = YOSHINAGAGUI returns the handle to a new YOSHINAGAGUI or the handle to
%      the existing singleton*.
%
%YOSHINAGAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%function named CALLBACK in YOSHINAGAGUI.M with the given input arguments.
%
%YOSHINAGAGUI('Property','Value',...) creates a new YOSHINAGAGUI or raises 
%the existing singleton*.  Starting from the left, property value pairs are
%applied to the GUI before YoshinagaGUI_OpeningFcn gets called.  An
%unrecognized property name or invalid value makes property application
%stop.  All inputs are passed to YoshinagaGUI_OpeningFcn via varargin.
%
%*See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help YoshinagaGUI

% Last Modified by GUIDE v2.5 19-Sep-2016 20:39:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @YoshinagaGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @YoshinagaGUI_OutputFcn, ...
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


% --- Executes just before YoshinagaGUI is made visible.
function YoshinagaGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to YoshinagaGUI (see VARARGIN)

% Choose default command line output for YoshinagaGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes YoshinagaGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = YoshinagaGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
set(handles.pressureinlet,'string','180');
axes(handles.axes2);
imstr = ['C:\Users\pc\OneDrive - IHC Merwede Holding B.V\desktop\Postprocess program\Yoshinaga Threephase Compare\pic.png']
matlabImage = imread(imstr);
image(matlabImage)
axis off
axis image


% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Geometory Set
set(handles.ambientdepth,'string','5.57');
set(handles.suctionpipe,'string','0.57');
set(handles.transferpipe,'string','6.16');
set(handles.riserdiameter,'string','40');
%Solid Set
set(handles.particlediameter,'string','6');
set(handles.soliddensity,'string','2650');
set(handles.solidflux,'string','0.035');
%Division
set(handles.nr2,'string','200');
set(handles.pressureinlet,'string','180');
axes(handles.axes2);
imstr = ['C:\Users\pc\OneDrive - IHC Merwede Holding B.V\desktop\Postprocess program\Yoshinaga Threephase Compare\pic.png']
matlabImage = imread(imstr);
image(matlabImage)






% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Inputs
D = str2num(get(handles.riserdiameter,'string'));
% [mm] Riser Diameter
L1= str2num(get(handles.ambientdepth,'string'));              
% [m] Ambient Pressure length
L2= str2num(get(handles.suctionpipe,'string'));              
% [m] SL Suction pipe length
L3 = str2num(get(handles.transferpipe,'string')); 
% [m] GSL three phase transport pipe length
Nr =str2num(get(handles.nr2,'string'));
LR = L2+L3;            % length of the left riser
N3 = L3/Nr;            % Devided transport pipe length into small sections
N2 = L2/Nr;            % Devided suction pipe length into small sections
Z3 = [L3:-N3:0];       % From suction inlet to the seasurface
Z2 = [LR:-N2:L3];      % from the sucntion inlet to the gas injection inlet 
Z = [Z2 Z3];           % Total length 
Patm = 100000          % [pa] Atmosphere Pressure 
R = 287.058            % [J/(KG.K)] Gas parameter
T = 293                % [K]  20C degree temperature
Pinlet = str2num(get(handles.pressureinlet,'string'))*1000;        
% [pa] Inlet pressure

apha = L3/(L1+L2);     % [-] Ratio between the riser/downcommer
ds = str2num(get(handles.particlediameter,'string'));              
% [mm] Particle diameter
rho_S = str2num(get(handles.soliddensity,'string'));          
% [kg/m^3] Solid density
rho_L = 1024;          % [kg/m^3] Liquid density

Vis_w = 10^-6;               % [m2/s] Liquid viscosity
g = 9.81;                    % [kg.m/s2] gravity accelearation
mu_vis =1.00162*10E-3;       %20 degree Water Dynamic viscosity (Pa.s) 
mu_G = 1.82*10E-5;           %20 degree Air Dynamic viscosity (Pa.s) 

J_G_inlet = 3.55;       % [m/s] Gas inlet flux(First assign)

C_S = 0.1;              % Solid Concentration (First assign)
C_G = 0.2;              % Gas Concentration (First assign)
C_L = 0.7;              % Liquid Concentration (First assign)    
Cls_L = C_L;            % Initial liquid concentration in two phase

epsi = 0.56;            % Inlet fitting loss
epsie =1;               % Entrance length loss

S = rho_S/rho_L;        % Parameter
arphaf  = 0.2;          % Particle shape effect check 
                        %(https://en.wikipedia.org/wiki/Sediment_transport)
cycle = 500;            % Cycle number

%% Particle settling velocity & Reynolds number and drag coefficient 
%(use wiki method)

D = D/1000;             % mm to m (c) Pipe Diameter
ds = ds/1000;           % mm to m (c) Particel Diameter
Area = 0.25*pi*D^2;             % Pipe Area [m^2]
delta = (rho_S - rho_L)/rho_L;  % Density difference

ust = delta*g*ds^2/(18*Vis_w+sqrt(0.75*delta*g*ds^3));
% [m/s]  (c)  Ferguson and Church, 2004	Settling velocity
Rep = ust*ds/Vis_w ;                                                
% Reynolds particle number
n = (4.7+0.41*Rep^0.75)/(1+0.175*Rep^0.75);                         
% exponent number
uh = ust*Cls_L^n;                                          
% Hindered settling velocity  C_S/(1-C_G)
Cd = 24/Rep*(1+0.173*Rep^0.657)+0.413/(1+1.63*10^4*Rep^-1.09);      
% Particle Drag coefficient

%  jspick = [0.025 0.063 0.095] This is used when Js needs to be specified
 J_S1 = [0.001 0.025 0.063]
%   for s = 1:1  % Cycle to choose the desired J_S
  
J_S = str2num(get(handles.solidflux,'string'));
  
for o = 1:80   % add up step for different J_G        
    J_G_inlet = 0.1*o;        % Each step is 0.1     
    PI = rho_L*g*Z3 + Patm ;  % Pressure at inlet 
    rho_G = PI*(R*T)^(-1);    % Density of the gas 
    J_G = J_G_inlet *Pinlet*(PI).^(-1); 
    % JG at ambient pressure converted to the pressure at the gas inlet 
    J_L = 0.1                 
    % Set the Liquid flux back to 0.1 for each circle
    
for kp = 1:cycle;    
    % Cycle to until the pressure is fullfilled
for kg = 1:cycle;    
    % Cycle to fullfill the Momentum
    
%% Three phase concentration calculation 
    PInew = 0*Z3;       % set the initial pressureto be zero
    
%% Two phase concentration calculation
for j2 = 1:30;          % circulation for 50
uh = ust*Cls_L^n;       % particle hindered settling velocity !!! 
                        % Check something might be not right here
                        
Cls_S = 1-Cls_L;        % Concentration of solid is the 
                        % 1-concentraiton of liquid                       
vr = uh;                % Asign hindered settling velocity to vr
% formula = 
% Cls_L^2 + (-Cls * vr - J_L - J_S)*vr^(-1)*Cls_L+ Cls * J_L*vr^(-1)
% b and c in this formula repersent what is working on the 
b = (-vr - J_L - J_S)*vr^(-1);
c = J_L*vr^(-1); 
p = [1 b c];
answer = roots(p);
% Solve the problem of this abc function
for i = 1:2                    
    % Pick up the only vaild result from the root function
if answer(i,1) < 1
    Cls_L = answer(i,1);    
    % Judge for each answer where it is within the ragnge of (0-1)
end
end
uh = ust*Cls_L^n;
% Apply it for the new hindering settling velocity
end

%% SATO Concentration calculation
% MMSL = rho_S*J_S + rho_L*J_L;
% MMLG = rho_L*J_L +rho_G*J_G;
% MMSLG = rho_S*J_S + rho_L*J_L +rho_G*J_G;

% Initial Setting
MM = rho_S*J_S + rho_L*J_L + mean(rho_G.*J_G); % ParameterMM

alpha_x = mean(rho_G.*J_G)/MM; % A_x

% rhoLS_3 = rho_L * mean(C_L)*(1-mean(C_G))^(-1) + 
% rho_S * mean(C_S)*(1-mean(C_G))^(-1);
rhoLS_3 = rho_L* C_L.*(1-C_G).^(-1) + rho_S * C_S.*(1-C_G).^(-1);

for j = 1:100;

C_G = 0.85*(1+0.4*rho_G.*rhoLS_3.^(-1).*(alpha_x^(-1)-1)...
    +0.6*rho_G.*rhoLS_3.^(-1).*(alpha_x^(-1)-1).*((rhoLS_3.*rho_G.^(-1)...
    +0.4*(alpha_x^(-1)-1)).*(1+0.4*(alpha_x^(-1)-1)).^(-1)).^(0.5)).^(-1);
% Gas concentration calculation

rho3 = rho_S.*C_S + rho_L.*C_L + rho_G.*C_G;
% Mixture density of three phase

rho_A = (rho3.*rhoLS_3.^(-1)).^(1.5).*rhoLS_3;
% Aparent Density 

U_SW = (1-C_S.*(1-C_G).^(-1)).^2.4.*(1-(ds/D)^2)...
    .*sqrt((rho_L.*rho_A.^(-1).*S-1)/(S-1)).*ust;
% Settling velocity

C = 1 + 0.2*exp(-5*C_S.*(1-C_G).^(-1));
% Parameter

US = C.*MM.*rho_A.^(-1) + U_SW;
% XXX Velocity

if  C_S - J_S.*US.^(-1) <10E-2  % Judge whether the C_S is changing anymore 
    C_S = J_S.*US.^(-1);
    C_L = 1- C_S - C_G;
    rhoLS_3 = rho_L* C_L.*(1-C_G).^(-1) + rho_S * C_S.*(1-C_G).^(-1);
    % define new density in three phase
    
%     rhoLS_3 = rho_L * mean(C_L)*(1-mean(C_G))^(-1) + ...
%       rho_S * mean(C_S)*(1-mean(C_G))^(-1);
    break
%   If the change of the Cs is small enough break out and stop loop
else
    C_S = J_S.*US.^(-1);
    C_L = 1- C_S - C_G;
    rhoLS_3 = rho_L* C_L.*(1-C_G).^(-1) + rho_S * C_S.*(1-C_G).^(-1);
%   If the Cs is still changeing then keep running the loop

%     rhoLS_3 = rho_L * mean(C_L)*(1-mean(C_G))^(-1) +...
%       rho_S * mean(C_S)*(1-mean(C_G))^(-1);
end

end

rhoLS_2 = rho_L*Cls_L + rho_S*Cls_S;
% Density of in two phase

%% Calculate Frictional pressure drop in two-phase water-solid flow

Re_LS = rho_L*(J_L + J_S)*D*(mu_vis)^(-1);   
% Pay attention to whether it is dynamic viscosity or kinematic viecosity
lamda_LS = 0.316*Re_LS^(-0.25);                 
% Friction factor liquid
FPD_2_dz = lamda_LS*D^(-1)*rhoLS_2*0.5*(J_L+J_S).^2;         
% Frictional pressure drop in two-phase water-solid flow 

%% Calculate Frictional pressure drop in three-phase water-solid flow

Re_G = rho_G.*J_G*D*mu_G^(-1);                 
% In the formula it use kinematic viscosity for expression.
lamda_G = 0.316*Re_G.^(-0.25);
% Friction factor gas
FPD_G_dz = lamda_G* D^(-1).*rho_G.*J_G.^2*0.5;
% Frictional pressure drop caused by gas

XX = sqrt(FPD_2_dz.*FPD_G_dz.^(-1));
%% Paramter XXX
FPD_3_dz = FPD_2_dz.*(1+21*XX.^(-1)+XX.^(-2));
% Frictional pressure drop in three-phase water-solid flow

%% Calculate Entrance pressure drop in two-phase water-solid flow

EPD_2 = (epsi+epsie)*mean(rhoLS_2)*0.5*(J_L+J_S)^2;              
% Entrance pressure drop

%% Calculate Entrance pressure drop in three-phase water-gas-solid flow

EPD_3 = rhoLS_3(1,1)*0.5.*((J_L+J_S)*(1-C_G(1,1))^(-1))^2 ...
    - rhoLS_2*0.5.*(J_L+J_S).^2;   
% Entrance pressure drop in three phase

%% Momentum Balance Equation

VL_E = J_L*C_L(1,1)^(-1) ;          % Entrance liquid velocity
VS_E = J_S*C_S(1,1)^(-1) ;          % Entrance solid velocity
VG_O = J_G(1,size(J_G,2))*C_G(1,size(C_G,2))^(-1);     
                                    % Entrance gas velocity
VS_O = J_S*C_S(1,size(C_S,2))^(-1);     % Outlet solid velocity
VL_O = J_L*C_L(1,size(C_L,2))^(-1);     % Outlet liquid velocity
J_G_O = J_G(1,size(J_G,2));             % Outlet gas flux
rho_G_O = rho_G(1,size(J_G,2));         % Outlet gas density

G1 = Area*(J_L*rho_L*VL_E+J_S*rho_S*VS_E);                                      
% Inlet momentum
G2 = -Area*(J_G_O*rho_G_O*VG_O+J_S*rho_S*VS_O+J_L*rho_L*VL_O);                  
% Outlet momentum
G3 = -Area*(FPD_2_dz*L2+EPD_2);                                                 
% Pressure lost in two phase and waterinlet
G4 = -Area*(sum(FPD_3_dz.*N3)+EPD_3);                                           
% Pressure lost in three phase and gasinlet
G5 = -Area*(rho_L*Cls_L + rho_S*Cls_S)*g*L2;                                    
% Mixture density in two phase
G6 = -Area*sum((rho_G.*C_G+C_L.*rho_L+C_S.*rho_S).*N3)*g;                       
% Mixture density in three phase
G7 = Area*rho_L*g*L1;                                                           
% Pressure provided by envirnoment

Gtotalnew = G1+G2+G3+G4+G5+G6+G7;                                               
% Momentum balance
        
if  abs(Gtotalnew) < 10                                                        
% Judge if the Gtotal is reaching the balance 0.5 as shrewhold(!!!!)
% very important factor determin whether the solution will converge or not
    PInew = rho3.*Z3*g + Patm ;                                                 
% If yes,asign the new density of mixture in the riser calculate new P(k)
    rho_G = PI*(R*T)^(-1);                                                      
% New air density at different location in the riser
    J_G = J_G_inlet *(Pinlet)*(PI).^(-1);                                       
% Different gas flux at different location in the rise
    break
else
    J_L = J_L + 0.01;                                                           
    % If no, asign a new J_L and put it back to calculation until yes
end
% kg
end
if  abs(PInew(1,1)-PI(1,1)) < 100                                               
% Judge if the P(inlet) is the same with the condition before(!!!!)
% very important factor determin whether the solution will converge or not
% 10 as shrewhold is because pressure is in pa Unit(e.g. 180000 pa)
    J_L_good(o,1) = J_G(1,1);                                                   
% Give the gas flux at the inlet to dispaly
    J_L_good(o,2) = J_L;                                                        
% Give the liquid flux at the inlet to dispaly
    QGstr = ['Qg'];QLstr = ['Ql'];QSstr = ['Qs'];
    QG{o} = J_G(1,1)*Area*1000;   %[dm^3/s]
    QL{o} = J_L*Area*1000;        %[dm^3/s]
    QS{o} = J_S*Area*1000;        %[dm^3/s]
    
    QLG{o} = QL{o}*QG{o}^(-1)
% Give the Gas and Liquid flowrate 
%  Efficiency
    eta(o,1) = rho_L*g*QLG{o}*(L2+L3-L1).*(Patm*log(Pinlet/Patm)).^(-1)
    break
% The calculation is good exit the loop
else
% If the situation doesn't meet the requirement then assign the 
    PI = PInew
% Assign the new pressure to riser 
    rho_G = PI*(R*T)^(-1);
% Recauclate the Gas density
    J_G = J_G_inlet*(Pinlet)*(PI).^(-1);
% Assigned the new J_G
    J_L = 0.1
% Reassigned to the 
end
kp
% display how much circulation taken
end
o
% Display what JG is in calculation
end

eta(1,1)= 0


[pks,locs] = findpeaks(J_L_good(:,2));
QGL={QGstr,QG{1,locs},'L/s';QLstr,QL{1,locs},'L/s';QSstr,QS{1,locs},'L/s'};

%% plot
axes(handles.axes1);
cla;
textstr = ['[',num2str(J_L_good(locs,1)),',',num2str(J_L_good(locs,2)),']']
hold on
grid on
yyaxis left
scatter(J_L_good(locs,1),J_L_good(locs,2),'filled','v')
plot([J_L_good(locs,1),J_L_good(locs,1)],[0,J_L_good(locs,2)],'-k');
plot([0,J_L_good(locs,1)],[J_L_good(locs,2),J_L_good(locs,2)],'-k');
text(J_L_good(locs,1)+0.02,J_L_good(locs,2)+0.2,textstr)
plot(J_L_good(:,1),J_L_good(:,2),'-b');

axis([0 10 0 1])

xlabel('Gas volumetric flux(JG) [m/s]');
ylabel('Liqiud volumetric flux(JL) [m/s]'); 

yyaxis right
plot(J_L_good(:,1), eta(:,1),'-r');
ylabel('Efficiency [-]'); 

set(handles.uitable2,'data',QGL);

% figure
% plot(J_L_good(:,1),J_L_good(:,2));
% axis([0 12.5 0 2])
% grid on
% xlabel('Gas volumetric flux(JG) [m/s]');
% ylabel('Liqiud volumetric flux(JL) [m/s]'); 







function pressureinlet_Callback(hObject, eventdata, handles)
% hObject    handle to pressureinlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pressureinlet as text
%str2double(get(hObject,'String')) 
...returns contents of pressureinlet as a double


% --- Executes during object creation, after setting all properties.
function pressureinlet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pressureinlet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function particlediameter_Callback(hObject, eventdata, handles)
% hObject    handle to particlediameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of particlediameter as text
% str2double(get(hObject,'String')) 
...returns contents of particlediameter as a double


% --- Executes during object creation, after setting all properties.
function particlediameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to particlediameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function soliddensity_Callback(hObject, eventdata, handles)
% hObject    handle to soliddensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of soliddensity as text
% str2double(get(hObject,'String')) 
...returns contents of soliddensity as a double


% --- Executes during object creation, after setting all properties.
function soliddensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to soliddensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function solidflux_Callback(hObject, eventdata, handles)
% hObject    handle to solidflux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of solidflux as text
%str2double(get(hObject,'String')) 
...returns contents of solidflux as a double


% --- Executes during object creation, after setting all properties.
function solidflux_CreateFcn(hObject, eventdata, handles)
% hObject    handle to solidflux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function riserdiameter_Callback(hObject, eventdata, handles)
% hObject    handle to riserdiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of riserdiameter as text
% str2double(get(hObject,'String'))
...returns contents of riserdiameter as a double


% --- Executes during object creation, after setting all properties.
function riserdiameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to riserdiameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor')...
        , get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ambientdepth_Callback(hObject, eventdata, handles)
% hObject    handle to ambientdepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ambientdepth as text
%        str2double(get(hObject,'String')) 
...returns contents of ambientdepth as a double


% --- Executes during object creation, after setting all properties.
function ambientdepth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ambientdepth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
    get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function suctionpipe_Callback(hObject, eventdata, handles)
% hObject    handle to suctionpipe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of suctionpipe as text
%        str2double(get(hObject,'String')) 
...returns contents of suctionpipe as a double


% --- Executes during object creation, after setting all properties.
function suctionpipe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to suctionpipe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function transferpipe_Callback(hObject, eventdata, handles)
% hObject    handle to transferpipe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transferpipe as text
%        str2double(get(hObject,'String')) 
...returns contents of transferpipe as a double


% --- Executes during object creation, after setting all properties.
function transferpipe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transferpipe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uibuttongroup1.
function uibuttongroup1_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when uibuttongroup1 is resized.
function uibuttongroup1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function nr2_Callback(hObject, eventdata, handles)
% hObject    handle to nr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nr2 as text
% str2double(get(hObject,'String')) returns contents of nr2 as a double


% --- Executes during object creation, after setting all properties.
function nr2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal...
 (get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nr3_Callback(hObject, eventdata, handles)
% hObject    handle to text30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text30 as text
% str2double(get(hObject,'String')) returns contents of text30 as a double


% --- Executes during object creation, after setting all properties.
function text30_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal...
 (get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
