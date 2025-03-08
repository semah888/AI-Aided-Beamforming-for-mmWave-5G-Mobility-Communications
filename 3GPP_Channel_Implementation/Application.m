classdef Application < matlab.apps.AppBase
    properties (Access = public)
        UIFigure            matlab.ui.Figure
        ScenarioDropdown    matlab.ui.control.DropDown
        PolarDropdown       matlab.ui.control.DropDown
        LOS_ConditionField  matlab.ui.control.DropDown
        h_BSField           matlab.ui.control.NumericEditField
        h_UTField           matlab.ui.control.NumericEditField
        xBS_0Field          matlab.ui.control.NumericEditField
        yBS_0Field          matlab.ui.control.NumericEditField
        xUT_0Field          matlab.ui.control.NumericEditField
        yUT_0Field          matlab.ui.control.NumericEditField
        Nrx_hField          matlab.ui.control.NumericEditField
        Nrx_vField          matlab.ui.control.NumericEditField
        Ntx_hField          matlab.ui.control.NumericEditField
        Ntx_vField          matlab.ui.control.NumericEditField
        SpeedUEField        matlab.ui.control.NumericEditField
        fcField             matlab.ui.control.NumericEditField
        BField              matlab.ui.control.NumericEditField
        TField              matlab.ui.control.NumericEditField
        UpdatePeriodField   matlab.ui.control.NumericEditField
        TitleSystemParamsLabel  % Titre pour la zone des paramètres du système
        RunButton           matlab.ui.control.Button
        UIAxes              matlab.ui.control.UIAxes
        TitleMetricsLabel   % Titre pour la zone des métriques
    TitleChannelLabel   % Titre pour la matrice du canal
    TitleFigureLabel    % Titre pour la zone des graphiques

        UITableMetrics  matlab.ui.control.Table
    UITableH        matlab.ui.control.Table
        UITable             matlab.ui.control.Table
        % Label fields (added these properties)
        h_BSLabel           matlab.ui.control.Label
        h_UTLabel           matlab.ui.control.Label
        xBS_0Label          matlab.ui.control.Label
        yBS_0Label          matlab.ui.control.Label
        xUT_0Label          matlab.ui.control.Label
        yUT_0Label          matlab.ui.control.Label
        Nrx_hLabel          matlab.ui.control.Label
        Nrx_vLabel          matlab.ui.control.Label
        Ntx_hLabel          matlab.ui.control.Label
        Ntx_vLabel          matlab.ui.control.Label
        SpeedUELabel        matlab.ui.control.Label
        fcLabel             matlab.ui.control.Label
        TLabel              matlab.ui.control.Label
        UpdatePeriodLabel   matlab.ui.control.Label
        Layout              matlab.ui.container.GridLayout % Add this line
            TimeLabel   % Étiquette pour afficher le temps actuel
    LocationLabel % Étiquette pour afficher la position actuelle
    snrLabel   matlab.ui.control.Label
    snrField   matlab.ui.control.NumericEditField
    end    

    methods (Access = public)
        function app = Application()
            % Create UI Figure & Components
            app.UIFigure = uifigure('Name', '3GPP_TR38.901_Generator', 'Position', [100, 100, 900, 600]);
            app.UIFigure.Visible = 'on';  % Make the application visible

            % Scenario Dropdown
            app.ScenarioDropdown = uidropdown(app.UIFigure, ...
                'Position', [20, 550, 50, 30], ...
                'Items', {'UMa', 'UMi', 'RMa'}, ...
                'Value', 'UMa');

            % Line of Sight (LoS) Dropdown
            app.LOS_ConditionField = uidropdown(app.UIFigure, 'Position', [140, 550, 100, 30], 'Items', {'LoS', 'NLoS'}, 'Value', 'LoS');

            % Label and Input Fields (adjusted for spacing)
            % Base Station height
            app.h_BSLabel = uilabel(app.UIFigure, 'Text', 'h_BS (m)', 'Position', [20, 500, 100, 30]);
            app.h_BSField = uieditfield(app.UIFigure, 'numeric', 'Position', [120, 500, 150, 30]);

            % User Terminal height
            app.h_UTLabel = uilabel(app.UIFigure, 'Text', 'h_UT (m)', 'Position', [20, 460, 100, 30]);
            app.h_UTField = uieditfield(app.UIFigure, 'numeric', 'Position', [120, 460, 150, 30]);

            % Base Station x and y coordinates
            app.xBS_0Label = uilabel(app.UIFigure, 'Text', 'xBS_0', 'Position', [20, 420, 100, 30]);
            app.xBS_0Field = uieditfield(app.UIFigure, 'numeric', 'Position', [120, 420, 150, 30]);
            app.yBS_0Label = uilabel(app.UIFigure, 'Text', 'yBS_0', 'Position', [20, 380, 100, 30]);
            app.yBS_0Field = uieditfield(app.UIFigure, 'numeric', 'Position', [120, 380, 150, 30]);

            % User Terminal x and y coordinates
            app.xUT_0Label = uilabel(app.UIFigure, 'Text', 'xUT_0', 'Position', [20, 340, 100, 30]);
            app.xUT_0Field = uieditfield(app.UIFigure, 'numeric', 'Position', [120, 340, 150, 30]);
            app.yUT_0Label = uilabel(app.UIFigure, 'Text', 'yUT_0', 'Position', [20, 300, 100, 30]);
            app.yUT_0Field = uieditfield(app.UIFigure, 'numeric', 'Position', [120, 300, 150, 30]);

            % Antenna Fields
            app.Nrx_hLabel = uilabel(app.UIFigure, 'Text', 'Nrx_h', 'Position', [300, 500, 100, 30]);
            app.Nrx_hField = uieditfield(app.UIFigure, 'numeric', 'Position', [400, 500, 150, 30]);
            app.Nrx_vLabel = uilabel(app.UIFigure, 'Text', 'Nrx_v', 'Position', [300, 460, 100, 30]);
            app.Nrx_vField = uieditfield(app.UIFigure, 'numeric', 'Position', [400, 460, 150, 30]);
            app.Ntx_hLabel = uilabel(app.UIFigure, 'Text', 'Ntx_h', 'Position', [300, 420, 100, 30]);
            app.Ntx_hField = uieditfield(app.UIFigure, 'numeric', 'Position', [400, 420, 150, 30]);
            app.Ntx_vLabel = uilabel(app.UIFigure, 'Text', 'Ntx_v', 'Position', [300, 380, 100, 30]);
            app.Ntx_vField = uieditfield(app.UIFigure, 'numeric', 'Position', [400, 380, 150, 30]);

            % Time-related Fields
            app.TLabel = uilabel(app.UIFigure, 'Text', 'Number Slots', 'Position', [300, 340, 100, 30]);
            app.TField = uieditfield(app.UIFigure, 'numeric', 'Position', [400, 340, 150, 30]);
            app.UpdatePeriodLabel = uilabel(app.UIFigure, 'Text', 'Update Period (s)', 'Position', [300, 300, 100, 30]);
            app.UpdatePeriodField = uieditfield(app.UIFigure, 'numeric', 'Position', [400, 300, 150, 30]);

            % Speed and Frequency Fields
            app.SpeedUELabel = uilabel(app.UIFigure, 'Text', 'Speed UE (km/h)', 'Position', [580, 500, 100, 30]);
            app.SpeedUEField = uieditfield(app.UIFigure, 'numeric', 'Position', [680, 500, 150, 30]);
            app.fcLabel = uilabel(app.UIFigure, 'Text', 'Carrier Frequency (Hz)', 'Position', [580, 460, 150, 30]);
            app.fcField = uieditfield(app.UIFigure, 'numeric', 'Position', [680, 460, 150, 30]);
            app.snrLabel = uilabel(app.UIFigure, 'Text', 'SNR transmission (Linear)', 'Position', [580, 420, 150, 30]);
            app.snrField = uieditfield(app.UIFigure, 'numeric', 'Position', [680, 420, 150, 30]);
          % Run Button
app.RunButton = uibutton(app.UIFigure, 'push', 'Text', 'Run', 'Position', [700, 350, 100, 50]);
app.RunButton.ButtonPushedFcn = @(btn, event) runSimulation(app);


            % Création de la première table (Métriques)
app.UITableMetrics = uitable(app.UIFigure, 'Position', [20, 50, 400, 100]);

% Création de la deuxième table (Matrice H)
app.UITableH = uitable(app.UIFigure, 'Position', [20, 175, 400, 100]); 
            % Axes for Graph
            app.UIAxes = uiaxes(app.UIFigure, 'Position', [450, 75, 400, 200]);
% Étiquette pour le temps actuel
app.TimeLabel = uilabel(app.UIFigure, 'Text', 'Time: 0 s', 'Position', [500, 50, 100, 25]);
app.TimeLabel.FontSize = 14;

% Étiquette pour la position actuelle de l'utilisateur
app.LocationLabel = uilabel(app.UIFigure, 'Text', 'Location: 0 m', 'Position', [600, 50, 100, 25]);
app.LocationLabel.FontSize = 14;

            % Remove the axis labels and grid
app.UIAxes.XTick = [];  % Remove X-axis ticks
app.UIAxes.YTick = [];  % Remove Y-axis ticks
            % Make the application visible

            % Titre pour la zone des métriques
app.TitleMetricsLabel = uilabel(app.UIFigure, 'Text', 'Evaluation Metrics', ...
    'Position', [20, 150, 400, 25], 'FontSize', 14, 'FontWeight', 'bold');

% Titre pour la zone de la matrice du canal
app.TitleChannelLabel = uilabel(app.UIFigure, 'Text', 'Channel Matrix', ...
    'Position', [20, 275, 400, 25], 'FontSize', 14, 'FontWeight', 'bold');

% Titre pour la zone des figures
app.TitleFigureLabel = uilabel(app.UIFigure, 'Text', 'Impulse Response', ...
    'Position', [450, 275, 400, 30], 'FontSize', 14, 'FontWeight', 'bold');
            app.UIFigure.Visible = 'on';

        end
    end

    methods (Access = private)
        function field = addLabeledField(app, labelText, x, y)
            % Ajoute un champ numérique avec un label
            uilabel(app.UIFigure, 'Text', labelText, 'Position', [x, y + 20, 120, 20]);
            field = uieditfield(app.UIFigure, 'numeric', 'Position', [x, y, 100, 30]);
        end
    end
  methods (Access = private)
      function runSimulation(app) %, h_BS_input, h_UT_input, xBS_0_input, yBS_0_input, xUT_0_input, yUT_0_input, Ntx_h_input, Ntx_v_input, Nrx_h_input, Nrx_v_input, LOS_input, UpdatePeriod_input, T_input, Polar_dropdown, Scenario_dropdown, speedUE_input, fc_input, B_input)

    % Get values from UI components
   Sc = app.ScenarioDropdown.Value
            %Polar = str2double(app.PolarDropdown.Value);
            LOS_Condition = app.LOS_ConditionField.Value;
            h_BS = app.h_BSField.Value;
            snr = app.snrField.Value;
            h_UT = app.h_UTField.Value;
            xBS_0 = app.xBS_0Field.Value;
            yBS_0 = app.yBS_0Field.Value;
            xUT_0 = app.xUT_0Field.Value;
            yUT_0 = app.yUT_0Field.Value;
            Nrx_h = app.Nrx_hField.Value;
            Nrx_v = app.Nrx_vField.Value;
            Ntx_h = app.Ntx_hField.Value;
            Ntx_v = app.Ntx_vField.Value;
            speedUE = app.SpeedUEField.Value * 1000 / 3600; % Convert km/h to m/s
            fc = app.fcField.Value;
          %  B = app.BField.Value;
            T = app.TField.Value;
            UpdatePeriod = app.UpdatePeriodField.Value;
if strcmp(Sc, 'UMa')   % Comparaison de chaîne de caractères
    scenario = 0;
elseif strcmp(Sc, 'UMi')
    scenario = 1;
else
    scenario = 2;   % Valeur par défaut si l'option n'est pas 'UMa' ou 'UMi'
end
Polar = 90;
if Polar == 0
 Polar_Slant_Angle = [0];
elseif Polar == 90
 Polar_Slant_Angle = [90];
elseif Polar == 45
 Polar_Slant_Angle = [45 -45];
end
%speedUE = 360;
speedUE = speedUE*1000/3600;
NumSC = 1;
c = 3e8;
lambda = c/fc;
d_rx = lambda/2;
d_tx = lambda/2;
%disp([Phi_LOS_AOA Phi_LOS_AOD Theta_LOS_ZOA Theta_LOS_ZOD])
i = 1;
pathLoss = zeros(1,T); % dB, 1 x T
Scenario = 2;
 velocity = speedUE*[1;0;0];
 %UpdatePeriod = 0.0005; %0.5 ms
 GridDimension = speedUE*UpdatePeriod; % 0.05 m
load("MIMOchannel_details.mat","PropagChanModelVariable");
 d_2D_ref = zeros(1,101);
 AODLOS_ref = zeros(1,101);
 ZODLOS_ref = zeros(1,101);
 % Initialize arrays for storing results
Location = zeros(T, 1);       % Array for location
Chan_db = zeros(T, 1);        % Array for channel gain in dB
SpectralEfficiency = zeros(T, 1); % Array for spectral efficiency
DelaySpread = zeros(T, 1);    % Array for delay spread
DopplerShift = zeros(T, 1);   % Array for Doppler shift
 for t = 1:101
 d_2D_ref(1,t) = PropagChanModelVariable{t}.ScalingValues.DistanceBSMS;
 ZODLOS_ref(1,t) = rad2deg(PropagChanModelVariable{t}.ScalingValues.ZLOSBSToMS);
 AODLOS_ref(1,t) = rad2deg(PropagChanModelVariable{t}.ScalingValues.ALOSBSToMS);
 end
 d_2D = zeros(1,T);
 xUT = zeros(1,T);
 for t = 0:T-1
 xUT(1,t+1) = xUT_0 + speedUE*t*UpdatePeriod;
 d_2D(1,t+1) = sqrt((abs(xUT(1,t+1))-xBS_0)^2+(yUT_0-yBS_0)^2);
 end
 % d_2D = d_2D_ref;
Phi_LOS_AOA = zeros(1,T);
Phi_LOS_AOD = zeros(1,T);
for t = 0:T-1
 vBS = [xBS_0;yBS_0;h_BS];
 vUT = [xUT(1,t+1);yUT_0;h_UT];
 vD = vUT - vBS;
 vA = vBS - vUT;
 if xUT(1,t+1)>=0
 Phi_LOS_AOD(1,t+1) = atand((yUT_0-yBS_0)/(xUT(1,t+1)-xBS_0));
 else
 Phi_LOS_AOD(1,t+1) = -180-atand((yUT_0-yBS_0)/(abs(xUT(1,t+1))-xBS_0));
 end
 Phi_LOS_AOA(1,t+1) = 180-Phi_LOS_AOD(1,t+1);
end
AOD_LOS = Phi_LOS_AOD;
%Phi_LOS_AODFinal = Phi_LOS_AOD
% initial LOS angles from data Q.Li
%Azimuth oriontation 0 at BS and 180 at MS
%Zenith oriontation 90 at BS and 90 at MS
%Phi_LOS_AOD = AODLOS_ref;
%Phi_LOS_AOA = 180-Phi_LOS_AOD;
Theta_LOS_ZOA = 90*ones(1,T+1);
Theta_LOS_ZOD = 90*ones(1,T+1);
d_3D = sqrt(d_2D.^2+(h_BS-h_UT)^2);
if Scenario == 2
 prLOS = exp(-(d_2D-10)/1000);
end
rho = 1; %Aleatoire value between [0 1]
LOSIndicator = ones(1,T);
Ntx = Ntx_h*Ntx_v;
Nrx = Nrx_h*Nrx_v;
Ant_orientation_tx_Az = 0;
Ant_orientation_tx_El = 90;
Ant_orientation_rx_Az = 180;
Ant_orientation_rx_El = 90;
[p,numPolarAngle] = size(Polar_Slant_Angle);
numItr = 1;
Channel_2 = zeros(numItr,numPolarAngle*Nrx,numPolarAngle*Ntx,T,NumSC);
freqChan = zeros(numPolarAngle*Nrx,numPolarAngle*Ntx,NumSC,T);
Ntx = Ntx*numPolarAngle;
Nrx = Nrx*numPolarAngle;
DataSet = zeros(numItr*T,22);
DataSet2 = zeros(numItr*T,26);
i=1;
for itr=1:numItr
DS_final = zeros(1,T);
ASA_final = zeros(1,T);
ASD_final = zeros(1,T);
ZSA_final = zeros(1,T);
ZSD_final = zeros(1,T);
DS = zeros(1,T);
ASD = zeros(1,T);
ZSD = zeros(1,T);
ASA = zeros(1,T);
ZSA = zeros(1,T);
k_factor = zeros(1,T);
SF = zeros(1,T);
initial = 1;
while initial<=T
 LOS_Condition = LOSIndicator(1,initial);
 t= initial;
 while t<T+1 && LOSIndicator(1,t) == LOS_Condition 
 N_slots_model = t;
 t = t+1;
 end
N_slots_model = N_slots_model-initial+1;
if Scenario == 0
 if LOS_Condition == 1
C_ASD_DS = 0.5;
 C_ASA_DS = 0.8;
 C_ASA_SF = -0.4;
 C_ASD_SF = -0.5;
 C_DS_SF = -0.4;
 C_ASD_ASA = 0.4;
 C_ASD_K = -0.2;
 C_ASA_K = -0.3;
 C_DS_K = -0.7;
 C_SF_K = 0.5;
 C_ZSD_DS = 0;
 C_ZSA_DS = 0.2;
 C_ZSA_SF = 0;
 C_ZSD_SF = 0;
 C_ZSD_ZSA = 0;
 C_ZSD_K = 0;
 C_ZSA_K = 0;
 C_ZSA_ASA = 0;
 C_ZSD_ASA = 0;
 C_ZSA_ASD = 0.3;
 C_ZSD_ASD = 0.5;
 dcorr_DS = 7;
 dcorr_K = 15;
 dcorr_SF = 10;
 dcorr_ASA = 8;
 dcorr_ASD = 8;
 dcorr_ZSA = 12;
 dcorr_ZSD = 12;
 mulgDS = -0.24*log10(1+fc/1e9)-7.14;
 sigmalgDS = 0.38;
 mulgASA = -0.08*log10(1+fc/1e9)+1.73;
 sigmalgASA = 0.014*log10(1+fc/1e9)+0.28;
 mulgASD = -0.05*log10(1+fc/1e9)+1.21;
 sigmalgASD = 0.41;
 mulgZSA = -0.1*log10(1+fc/1e9)+0.073;
 sigmalgZSA = -0.04*log10(1+fc/1e9)+0.34;
 mulgZSD = max(-0.21,-14.8*(min(d_2D)/1000)+0.01*abs(h_UT-h_BS)+0.83);
 sigmalgZSD = 0.35;
 muK = 9;
 sigmaK = 5;
 muSF = 0;
 sigmaSF = 4;
 r_taux = 3;
mu_XPR = 9;
sigma_XPR = 3;
N_Cluster = 12;
N_Ray = 20;
C_DS = 5;
C_ASA = 17;
C_ASD = 3;
C_ZSA = 7;
Per_Cluster_Shadowing = 3;
mu_ZOD_offset = 0; 
d_BP = 4*h_BS*h_UT*fc/c;
for t=initial:initial+N_slots_model-1
 if 10 <=d_2D(1,t)<=d_BP
 pathLoss(1,t) = 32.4+21*log10(d_3D(1,t))+20*log10(fc/1e9);
 else
 pathLoss(1,t) = 32.4+40*log10(d_3D(1,t))+20*log10(fc/1e9)-9.5*log10(d_BP^2+(h_BS-h_UT)^2);
 end
end
elseif LOS_Condition == 0
 C_ASD_DS = 0;
 C_ASA_DS = 0.4;
 C_ASA_SF = -0.4;
 C_ASD_SF = 0;
 C_DS_SF = -0.7;
 C_ASD_ASA = 0;
 C_ASD_K = -0.2;
 C_ASA_K = -0.3;
 C_DS_K = -0.7;
 C_SF_K = 0.5;
 C_ZSD_DS = -0.5;
 C_ZSA_DS = 0;
 C_ZSA_SF = 0;
 C_ZSD_SF = 0;
 C_ZSD_ZSA = 0;
 C_ZSD_K = 0;
 C_ZSA_K = 0;
 C_ZSA_ASA = 0.2;
 C_ZSD_ASA = 0;
 C_ZSA_ASD = 0.5;
 C_ZSD_ASD = 0.5;
 dcorr_DS = 10;
 dcorr_K = 15;
 dcorr_SF = 13;
 dcorr_ASA = 9;
 dcorr_ASD = 10;
 dcorr_ZSA = 10;
 dcorr_ZSD = 10;
 mulgDS = -0.24*log10(1+fc/1e9)-6.83;
 sigmalgDS = 0.16*log10(1+fc/1e9)+0.28;
 mulgASA = -0.08*log10(1+fc/1e9)+1.81;
 sigmalgASA = 0.05*log10(1+fc/1e9)+0.3;
 mulgASD = -0.23*log10(1+fc/1e9)+1.53;
 sigmalgASD = 0.11*log10(1+fc/1e9)+0.33;
 mulgZSA = -0.04*log10(1+fc/1e9)+0.92;
 sigmalgZSA = -0.07*log10(1+fc/1e9)+0.41;
 mulgZSD = max(-0.21,-14.8*(min(d_2D)/1000)+0.01*abs(h_UT-h_BS)+0.83);
 sigmalgZSD = 0.35;
 muK = 9;
 sigmaK = 5;
 muSF = 0;
 sigmaSF = 7.82;
 r_taux = 2.1;
mu_XPR = 8;
sigma_XPR = 3;
N_Cluster = 19;
N_Ray = 20;
C_DS = 11;
C_ASA = 22;
C_ASD = 10;
C_ZSA = 7;
Per_Cluster_Shadowing = 3;
mu_ZOD_offset = 0; 
pathLoss_LOS = zeros(1,N_slots_model+1); % dB, 1 x T
d_BP = 4*h_BS*h_UT*fc/c;
for t=initial:initial+N_slots_model-1
 if 10 <=d_2D(1,t)<=d_BP
 pathLoss_LOS(1,t-initial+1) = 32.4+21*log10(d_3D(1,t))+20*log10(fc/1e9);
 else
 pathLoss_LOS(1,t-initial+1) = 32.4+40*log10(d_3D(1,t))+20*log10(fc/1e9)-9.5*log10(d_BP^2+(h_BS-h_UT)^2);
 end
end
pathLoss(1,initial:initial+N_slots_model-1) = max(pathLoss_LOS,22.4+35.3*log10(d_3D)+21.3*log10(fc/1e9)-0.3*(h_UT-1.5));
 end
elseif Scenario == 1
 if LOS_Condition == 0
 C_ASD_DS = 0.4;
 C_ASA_DS = 0.6;
 C_ASA_SF = 0;
 C_ASD_SF = -0.6;
 C_DS_SF = -0.4;
 C_ASD_ASA = 0.4;
 C_ASD_K = -0.2;
 C_ASA_K = -0.3;
 C_DS_K = -0.7;
 C_SF_K = 0.5;
 C_ZSD_DS = -0.5;
 C_ZSA_DS = 0;
 C_ZSA_SF = -0.4;
 C_ZSD_SF = 0;
 C_ZSD_ZSA = 0;
 C_ZSD_K = 0;
 C_ZSA_K = 0;
 C_ZSA_ASA = 0;
 C_ZSD_ASA = 0;
 C_ZSA_ASD = -0.1;
 C_ZSD_ASD = 0.5;
 dcorr_DS = 40;
 dcorr_K = 15;
 dcorr_SF = 50;
 dcorr_ASA = 50;
 dcorr_ASD = 50;
 dcorr_ZSA = 50;
 dcorr_ZSD = 50;
 mulgDS = -6.28-0.204*log10(fc/1e9);
 sigmalgDS = 0.39;
 mulgASA = 2.08-0.27*log10(fc/1e9);
 sigmalgASA = 0.11;
 mulgASD = 1.5-0.1144*log10(fc/1e9);
 sigmalgASD = 0.28;
 mulgZSA = -0.3236*log10(fc/1e9)+1.512;
 sigmalgZSA = 0.16;
 mulgZSD = max(-0.5,-2.1*(min(d_2D)/1000)-0.01*abs(h_UT-1.5)+0.9);
 sigmalgZSD = 0.49;
 muK = 9;
 sigmaK = 5;
 muSF = 0;
 sigmaSF = 6;
 r_taux = 2.3;
mu_XPR = 7;
sigma_XPR = 3;
N_Cluster = 20;
N_Ray = 20;
C_DS = 0;
C_ASA = 15;
C_ASD = 2;
C_ZSA = 7;
Per_Cluster_Shadowing = 3;
mu_ZOD_offset = 7.66*log10(fc/1e9)-5.96-10^((0.208*log10(fc/1e9)-0.782)*log10(max(25,max(d_2D)))+(-0.13*log10(fc/1e9)+2.03)-0.07*(h_UT-1.5)); 
d_BP = 4*h_BS*h_UT*fc/c;
for t=1:T+1
 if 10 <=d_2D(1,t)<=d_BP
 pathLoss_LOS(1,t) = 28+22*log10(d_3D(1,t))+20*log10(fc/1e9);
 else
 pathLoss_LOS(1,t) = 28+40*log10(d_3D(1,t))+20*log10(fc)-9*log10(d_BP^2+(h_BS-h_UT)^2);
 end
end
pathLoss(1,initial:initial+N_slots_model-1) = max(pathLoss_LOS,13.54+39.08*log10(d_3D)+20*log10(fc/1e9)-0.6*(h_UT-1.5)); 
 elseif LOS_Condition == 1
 C_ASD_DS = 0.4;
 C_ASA_DS = 0.8;
 C_ASA_SF = -0.5;
 C_ASD_SF = -0.5;
 C_DS_SF = -0.4;
 C_ASD_ASA = 0;
 C_ASD_K = 0;
 C_ASA_K = -0.2;
 C_DS_K = -0.4;
 C_SF_K = 0;
 C_ZSD_DS = -0.2;
 C_ZSA_DS = 0;
 C_ZSA_SF = -0.8;
 C_ZSD_SF = 0;
 C_ZSD_ZSA = 0;
 C_ZSD_K = 0;
 C_ZSA_K = 0;
 C_ZSA_ASA = 0.4;
 C_ZSD_ASA = -0.3;
 C_ZSA_ASD = 0;
 C_ZSD_ASD = 0.5;
 dcorr_DS = 30;
 dcorr_K = 12;
 dcorr_SF = 37;
 dcorr_ASA = 15;
 dcorr_ASD = 18;
 dcorr_ZSA = 15;
 dcorr_ZSD = 15;
 mulgDS = -6.955-0.0963*log10(fc/1e9);
 sigmalgDS = 0.66;
 mulgASA = 1.81;
 sigmalgASA = 0.2;
 mulgASD = 1.06-0.1144*log10(fc/1e9);
 sigmalgASD = 0.28;
 mulgZSA = 0.95;
 sigmalgZSA = 0.95;
 mulgZSD = max(-0.5,-2.1*(min(d_2D)/1000)-0.01*abs(h_UT-1.5)+0.75);
 sigmalgZSD = 0.4;
 muK = 9;
 sigmaK = 3.5;
 muSF = 0;
 sigmaSF = 4;
 r_taux = 2.5;
mu_XPR = 8;
sigma_XPR = 4;
N_Cluster = 12;
N_Ray = 20;
C_DS = 0;
C_ASA = 11;
C_ASD = 5;
C_ZSA = 7;
Per_Cluster_Shadowing = 3;
mu_ZOD_offset = 0;
d_BP = 4*h_BS*h_UT*fc/c;
for t=1:T+1
 if 10 <=d_2D(1,t)<=d_BP
 pathLoss(1,initial:initial+N_slots_model-1) = 28+22*log10(d_3D(1,t))+20*log10(fc/1e9);
 else
 pathLoss(1,initial:initial+N_slots_model-1) = 28+40*log10(d_3D(1,t))+20*log10(fc/1e9)-9*log10(d_BP^2+(h_BS-h_UT)^2);
 end
end
 end
elseif Scenario ==2
 if LOS_Condition == 1
 C_ASD_DS = 0;
 C_ASA_DS = 0;
 C_ASA_SF = 0;
 C_ASD_SF = 0;
 C_DS_SF = -0.5;
 C_ASD_ASA = 0;
 C_ASD_K = 0;
 C_ASA_K = 0;
 C_DS_K = 0;
 C_SF_K = 0;
 C_ZSD_DS = -0.05;
 C_ZSA_DS = 0.27;
 C_ZSA_SF = -0.17;
 C_ZSD_SF = 0.01;
 C_ZSD_ZSA = -0.07;
 C_ZSD_K = 0;
 C_ZSA_K = -0.02;
 C_ZSA_ASA = 0.24;
 C_ZSD_ASA = -0.2;
 C_ZSA_ASD = -0.14;
 C_ZSD_ASD = 0.73;
 dcorr_DS = 50;
 dcorr_K = 40;
 dcorr_SF = 37;
 dcorr_ASA = 35;
 dcorr_ASD = 25;
 dcorr_ZSA = 15;
 dcorr_ZSD = 15;
 mulgDS = -7.49;
 sigmalgDS = 0.55;
 mulgASA = 1.52;
 sigmalgASA = 0.24;
 mulgASD = 0.9;
 sigmalgASD = 0.38;
 mulgZSA = 0.47;
 sigmalgZSA = 0.4;
 mulgZSD = max(-1,-0.17*(min(d_2D)/1000)-0.01*abs(h_UT-1.5)+0.22);
 sigmalgZSD = 0.34;
 muK = 7;
 sigmaK = 4;
 muSF = 0;
 sigmaSF = 4;
 r_taux = 3.8;
mu_XPR = 12;
sigma_XPR = 4;
N_Cluster = 11;
N_Ray = 20;
C_DS = 0;
C_ASA = 3;
C_ASD = 2;
C_ZSA = 3;
Per_Cluster_Shadowing = 3;
mu_ZOD_offset = 0; 
d_BP = 4*h_BS*h_UT*fc/c;
for t=initial:initial+N_slots_model-1
 if 10 <=d_2D(1,t)<=d_BP
 pathLoss(1,t) = 20*log10(40*pi*d_3D(1,t)*(fc/1e9)/3)+min(0.03*5^1.72,10)*log10(d_3D(1,t))-min(0.044*5^1.72,14.77)+0.002*log10(5)*d_3D(1,t);
 else
 pathLoss(1,t) = 20*log10(40*pi*d_BP*(fc/1e9)/3)+min(0.03*5^1.72,10)*log10(d_BP)-min(0.044*5^1.72,14.77)+0.002*log10(5)*d_BP(1,t)+40*log10(d_3D(1,t)/d_BP);
 end
end
 elseif LOS_Condition == 0
 C_ASD_DS = -0.4;
 C_ASA_DS = 0;
 C_ASA_SF = 0;
 C_ASD_SF = 0.6;
 C_DS_SF = -0.5;
 C_ASD_ASA = 0;
 C_ASD_K = 0;
 C_ASA_K = 0;
 C_DS_K = 0;
 C_SF_K = 0;
 C_ZSD_DS = -0.1;
 C_ZSA_DS = -0.4;
 C_ZSA_SF = -0.25;
 C_ZSD_SF = -0.04;
 C_ZSD_ZSA = -0.27;
 C_ZSD_K = 0;
 C_ZSA_K = 0;
 C_ZSA_ASA = 0.26;
 C_ZSD_ASA = -0.18;
 C_ZSA_ASD = -0.27;
 C_ZSD_ASD = 0.42;
 dcorr_DS = 36;
 dcorr_K = 40;
 dcorr_SF = 120;
 dcorr_ASA = 40;
 dcorr_ASD = 30;
 dcorr_ZSA = 50;
 dcorr_ZSD = 50;
 mulgDS = -7.43;
 sigmalgDS = 0.48;
 mulgASA = 1.52;
 sigmalgASA = 0.13;
 mulgASD = 0.95;
 sigmalgASD = 0.45;
 mulgZSA = 0.58;
 sigmalgZSA = 0.37;
 mulgZSD = max(-1,-0.19*(min(d_2D)/1000)-0.01*abs(h_UT-1.5)+0.28);
 sigmalgZSD = 0.3;
 muK = 0;
 sigmaK = 0;
 muSF = 0;
 sigmaSF = 8;
 r_taux = 1.7;
mu_XPR = 7;
sigma_XPR = 3;
N_Cluster = 10;
N_Ray = 20;
C_DS = 0;
C_ASA = 3;
C_ASD = 2;
C_ZSA = 3;
Per_Cluster_Shadowing = 3;
mu_ZOD_offset = 0; 
d_BP = 4*h_BS*h_UT*fc/c;
pathLoss_LOS = zeros(1,N_slots_model);
for t=initial:initial+N_slots_model-1
 if 10 <=d_2D(1,t)<=d_BP
 pathLoss_LOS(1,t-initial+1) = 20*log10(40*pi*d_3D(1,t)*(fc/1e9)/3)+min(0.03*5^1.72,10)*log10(d_3D(1,t))-min(0.044*5^1.72,14.77)+0.002*log10(5)*d_3D(1,t);
 else
 pathLoss_LOS(1,t-initial+1) = 20*log10(40*pi*d_BP*(fc/1e9)/3)+min(0.03*5^1.72,10)*log10(d_BP)-min(0.044*5^1.72,14.77)+0.002*log10(5)*d_BP(1,t)+40*log10(d_3D(1,t)/d_BP);
 end
end
pathLoss(1,initial:initial+N_slots_model-1) = max(pathLoss_LOS,161.04+7.1*log10(20)+7.5*log10(5)-(24.37-3.7*(5/h_BS)^2)*log10(h_BS)+(43.42-3.1*log10(h_BS))*(log10(d_3D(1,initial:initial+N_slots_model-1))-3)+20*log10(fc/1e9)-(3.2*(log10(11.75*h_UT))^2-4.97));
 end
end
%disp(['Large scale parameters']);
[DS(1,initial:initial+N_slots_model-1),ASD(1,initial:initial+N_slots_model-1),ASA(1,initial:initial+N_slots_model-1),ZSD(1,initial:initial+N_slots_model-1),ZSA(1,initial:initial+N_slots_model-1),k_factor(1,initial:initial+N_slots_model-1),SF(1,initial:initial+N_slots_model-1)] = LSPGenerator(GridDimension,N_slots_model,mulgDS,sigmalgDS,mulgASD,sigmalgASD,mulgASA,sigmalgASA,mulgZSD,sigmalgZSD,mulgZSA,sigmalgZSA,muSF,sigmaSF,muK,sigmaK,C_ASD_DS,C_ASA_DS,C_ASA_SF,C_ASD_SF,C_DS_SF,C_ASD_ASA,C_ASD_K,C_ASA_K,C_DS_K,C_SF_K,C_ZSD_DS,C_ZSA_DS,C_ZSA_SF,C_ZSD_SF,C_ZSD_ZSA,C_ZSD_K,C_ZSA_K,C_ZSA_ASA,C_ZSD_ASA,C_ZSA_ASD,C_ZSD_ASD,dcorr_DS,dcorr_K,dcorr_SF,dcorr_ASA,dcorr_ASD,dcorr_ZSA,dcorr_ZSD);
%Model = LOS_Condition
%disp([' mulgDS mulgASD mulgASA mulgZSD mulgZSA muk muSF'])
%disp(real([mean(log10(DS(1,initial:initial+N_slots_model-1))) mean(log10(ASD(1,initial:initial+N_slots_model-1))) mean(log10(ASA(1,initial:initial+N_slots_model-1))) mean(log10(ZSD(1,initial:initial+N_slots_model-1))) mean(log10(ZSA(1,initial:initial+N_slots_model-1))) mean(k_factor(1,initial:initial+N_slots_model-1)) mean(SF(1,initial:initial+N_slots_model-1))]))
LSP = zeros(N_slots_model,7);
for t = 0:N_slots_model-1
 LSP(t+1,1) = (DS(1,t+1));
 LSP(t+1,2) = (ASD(1,t+1));
 LSP(t+1,3) = (ASA(1,t+1));
 LSP(t+1,4) = (ZSD(1,t+1));
 LSP(t+1,5) = ZSA(1,t+1);
 LSP(t+1,6) = k_factor(1,t+1);
 LSP(t+1,7) = SF(1,t+1);
end
%figure(initial)
%corrplot(LSP)
disp(["Large Scale Parameters"]);
lsp = LSP;
disp(["Small Scale Parameters"]);
Delays = DelaysGenerator(N_slots_model,r_taux,DS(1,initial:initial+N_slots_model-1),N_Cluster,k_factor(1,initial:initial+N_slots_model-1));
Powers = PowerGenerator(Delays,DS(1,initial:initial+N_slots_model-1),r_taux,Per_Cluster_Shadowing,LOS_Condition,k_factor(1,initial:initial+N_slots_model-1));
Powers_db = mag2db(Powers)
Phi_AOA = AOAGenerator(Powers,ASA(1,initial:initial+N_slots_model-1),LOS_Condition,k_factor(1,initial:initial+N_slots_model-1),Phi_LOS_AOA(1,initial:initial+N_slots_model-1));
Phi_AOA_Rays = AOA_RaysGenerator(Phi_AOA,C_ASA,N_Ray);
Phi_AOD = AODGenerator(Powers,ASD(1,initial:initial+N_slots_model-1),LOS_Condition,k_factor(1,initial:initial+N_slots_model-1),Phi_LOS_AOD(1,initial:initial+N_slots_model-1));
Phi_AOD_Rays = AOD_RaysGenerator(Phi_AOD,C_ASD,N_Ray);
Theta_ZOA = ZOAGenerator(Powers,ZSA(1,initial:initial+N_slots_model-1),LOS_Condition,k_factor(1,initial:initial+N_slots_model-1),Theta_LOS_ZOA(1,initial:initial+N_slots_model-1));
Theta_ZOA_Rays = ZOA_RaysGenerator(Theta_ZOA,C_ZSA,N_Ray);
Theta_ZOD = ZODGenerator(Powers,ZSD(1,initial:initial+N_slots_model-1),LOS_Condition,k_factor(1,initial:initial+N_slots_model-1),Theta_LOS_ZOD(1,initial:initial+N_slots_model-1),mu_ZOD_offset);
Theta_ZOD_Rays = ZOD_RaysGenerator(Theta_ZOD,N_Ray,mulgZSD);
XPR = XPRGenerator(N_slots_model,N_Ray,N_Cluster,mu_XPR,sigma_XPR);
for t =initial:initial+N_slots_model-1
 for k=1:NumSC
 [channel_cluster,Chan] = ChannelCoefficientGenerator(t,k,initial,N_slots_model,NumSC,Ntx_v,Ntx_h,Nrx_v,Nrx_h,d_rx,d_tx,Polar_Slant_Angle,N_Cluster,N_Ray,Delays,Powers,Theta_ZOA_Rays,Phi_AOA_Rays,Theta_ZOD_Rays,Phi_AOD_Rays,lambda,XPR,velocity,C_DS,LOS_Condition,k_factor,Phi_LOS_AOA,Phi_LOS_AOD,Theta_LOS_ZOA,Theta_LOS_ZOD,d_3D,UpdatePeriod);
  n = size(channel_cluster);
 Channel_2(itr,:,:,t,k) = Chan;
 end
 %find the DS from the channel
 mean_delay = 0;
 mean_AOA = 0;
 mean_ZOA = 0;
 mean_delay_2 = 0;
 mean_AOA_2 = 0;
 mean_AOD_2 = 0;
 mean_AOD = 0;
 mean_ZOA_2 = 0;
 mean_ZOD_2 = 0;
 mean_ZOD = 0;
 somme = 0;
 for n=1:N_Cluster
 mean_delay = mean_delay+norm(channel_cluster(:,:,n))^2*Delays(n,t-initial+1);
 mean_AOA = mean_AOA+norm(channel_cluster(:,:,n))^2*Phi_AOA(n,t-initial+1);
 mean_AOD = mean_AOD+norm(channel_cluster(:,:,n))^2*Phi_AOD(n,t-initial+1);
 mean_ZOA = mean_ZOA+norm(channel_cluster(:,:,n))^2*Theta_ZOA(n,t-initial+1);
 mean_ZOD = mean_ZOD+norm(channel_cluster(:,:,n))^2*Theta_ZOD(n,t-initial+1);
 mean_delay_2 = mean_delay_2+norm(channel_cluster(:,:,n))^2*Delays(n,t-initial+1)^2;
 mean_AOA_2 = mean_AOA_2+norm(channel_cluster(:,:,n))^2*Phi_AOA(n,t-initial+1)^2;
 mean_AOD_2 = mean_AOD_2+norm(channel_cluster(:,:,n))^2*Phi_AOD(n,t-initial+1)^2;
 mean_ZOA_2 = mean_ZOA_2+norm(channel_cluster(:,:,n))^2*Theta_ZOA(n,t-initial+1)^2;
 mean_ZOD_2 = mean_ZOD_2+norm(channel_cluster(:,:,n))^2*Theta_ZOD(n,t-initial+1)^2;
 somme = somme+norm(channel_cluster(:,:,n))^2;
 %channel_cluster(:,:,n)
 end
 mean_delay_2 = mean_delay_2/somme;
 mean_delay = mean_delay/somme;
 mean_AOA_2 = mean_AOA_2/somme;
 mean_AOA = mean_AOA/somme;
 mean_AOD_2 = mean_AOD_2/somme;
 mean_AOD = mean_AOD/somme;
 mean_ZOA_2 = mean_ZOA_2/somme;
 mean_ZOA = mean_ZOA/somme;
 mean_ZOD_2 = mean_ZOD_2/somme;
 mean_ZOD = mean_ZOD/somme;
 for n=1:N_Cluster
 ASA_final(1,t) = ASA_final(1,t)+(Phi_AOA(n,t-initial+1)-mean_AOA)^2*norm(channel_cluster(:,:,n))^2/somme;
 ASD_final(1,t) = ASD_final(1,t)+(Phi_AOD(n,t-initial+1)-mean_AOD)^2*norm(channel_cluster(:,:,n))^2/somme;
 end
 DS_final(1,t) = sqrt(mean_delay_2-mean_delay^2);
 ASA_final(1,t) = sqrt(mean_AOA_2-mean_AOA^2);
 ASD_final(1,t) = sqrt(mean_AOA_2-mean_AOA^2);
 ZSA_final(1,t) = sqrt(mean_ZOA_2-mean_ZOA^2);
 ZSD_final(1,t) = sqrt(mean_ZOA_2-mean_ZOA^2);
end
initial = initial+N_slots_model
end
end
for t=1:T
    for itr =1:numItr
 k=1;
 chan = zeros(Ntx,Nrx);
 for nt = 1:Ntx
 for nr=1:Nrx
 chan(nt,nr) = Channel_2(itr,nt,nr,t,k);
 end
 end
 % Update location (example: moving UE)
    Location(t) = xUT_0 - speedUE * t * UpdatePeriod;
    H = chan/norm(chan);
    % Update Channel Gain (use your H matrix or channel model here)
    H_norm = norm(H, 'fro'); % Frobenius norm of the channel matrix
    Chan_db(t) = 10 * log10(H_norm^2); % Convert to dB
    noisePower = pathLoss(1,t)
    % Calculate Spectral Efficiency
    % Assuming SNR is calculated based on the Frobenius norm of H:
     % Example SNR calculation (use appropriate noise power)

    SpectralEfficiency(t) = log2(1 + snr*H_norm^2);  % Spectral efficiency (bps/Hz)
    
    
    % Update the plot (e.g., Channel Gain over time)
    stem(app.UIAxes, Delays(:,t), Powers(:,t));
    xlabel(app.UIAxes, 'Location (m)');
    ylabel(app.UIAxes, 'Altitude');
 
    % Update the table with current values
    metrics = table(Location(t), Chan_db(t), SpectralEfficiency(t), ...
                    'VariableNames', {'Location', 'Channel Gain (dB)', 'Spectral Efficiency (bps/Hz)'});
      app.UITableMetrics.Data = metrics;  
    
    
    % Add the matrix H as a string (formatted) to the table
    % Prepare an empty row for space (empty string)
    % **Créer une ligne vide avec 3 colonnes**
   % empty_row = cell2table(cell(1, width(metrics)), 'VariableNames', metrics.Properties.VariableNames);
   
    % Prepare a separate row for the channel matrix
   % H_string = mat2str(H); % Convert matrix H to string
    %H_row = cell2table({H_string, '', ''}, 'VariableNames', metrics.Properties.VariableNames);

    
    % Combine metrics, space, and channel matrix into one table
   % fullTable = [metrics; empty_row; H_row];  % Concatenate the metrics, space, and matrix H
    
    % Update the UITable with the full data (metrics + channel matrix)
   H_string = mat2str(H); % Convertir H en texte
    H_table = table({H_string}, 'VariableNames', {'Channel Matrix H'});
    app.UITableH.Data = H;%_table;
     % Mettre à jour le temps et la position de l'utilisateur


    currentTime = t * UpdatePeriod;  % Temps actuel
    currentLocation = xUT_0 - speedUE * t * UpdatePeriod;  % Position actuelle

    % Afficher ces valeurs dans l'interface
    app.TimeLabel.Text = sprintf('Time: %.2f s', currentTime);
    app.LocationLabel.Text = sprintf('Location: %.2f m', currentLocation);
    % Pause for a short time to simulate real-time plotting (optional)
    pause(10);  % Adjust the pause time as needed
    end
end
%plot(app.UIAxes, Location, Chan_db)
%title(app.UIAxes, 'Graphique de la Transmission'); % Titre du graphique
%    xlabel(app.UIAxes, 'Emplacement (m)'); % Étiquette de l'axe X
%    ylabel(app.UIAxes, 'Gain en dB'); % Étiquette de l'axe Y
%disp(["Impulse response at slot t"])
%t = input("time slot =  ");
%stem(app.UIAxes,Delays(:,t),Powers(:,t),'LineWidth',2);
%xlabel(app.UIAxes,'Delay (s)');
%ylabel(app.UIAxes,'Amplitude (Watt)');
end
  end
end