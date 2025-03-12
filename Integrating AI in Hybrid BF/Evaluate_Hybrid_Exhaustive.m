h_BS = 2.5 ; 
h_UT = 2.5 ; 
xBS_0 = 0;
yBS_0 = 5;
xUT_0 = -25;
yUT_0 = 0;
Nrx_h = 2;
Nrx_v = 2;
Ntx_h = 2;
Ntx_v = 2;
numBeams=3;
Ns = 1;
NumberGrids = 1001;
T = NumberGrids;
Polar = 90; %v_polarization
if Polar == 0
 Polar_Slant_Angle = [0];
elseif Polar == 90
 Polar_Slant_Angle = [90];
elseif Polar == 45
 Polar_Slant_Angle = [45 -45];
end
%speedUE = input( ' UE speed (km/h) = ' ) ;
speedUE = 360;
speedUE = speedUE*1000/3600;
fc = 30e9;
B = 60e6;
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
 UpdatePeriod = 0.0005; %0.5 ms
 GridDimension = speedUE*UpdatePeriod; % 0.05 m
 d_2D_ref = zeros(1,101);
 AODLOS_ref = zeros(1,101);
 ZODLOS_ref = zeros(1,101);
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
LOSIndicator = zeros(1,T);
Ntx = Ntx_h*Ntx_v;
Nrx = Nrx_h*Nrx_v;
Ant_orientation_tx_Az = 0;
Ant_orientation_tx_El = 90;
Ant_orientation_rx_Az = 180;
Ant_orientation_rx_El = 90;
[p,numPolarAngle] = size(Polar_Slant_Angle);
numItr = 150;
Channel_2 = zeros(numItr,numPolarAngle*Nrx,numPolarAngle*Ntx,T,NumSC);
freqChan = zeros(numPolarAngle*Nrx,numPolarAngle*Ntx,NumSC,T);
Ntx = Ntx*numPolarAngle;
Nrx = Nrx*numPolarAngle;
selected = zeros(T,6);
DataSet = zeros(numItr*T,25);
DataSet2 = zeros(numItr*T,22);
 [beamBS_opt,beamAzBS_opt,beamElBS_opt,beamAngleBS,beamAngleAzBS,beamAngleElBS] = getDFTCodebook3D(Ntx_h,Ntx_v,8,8,0.5,0.5);
 [beamUE_opt,beamAzUE_opt,beamElUE_opt,beamAngleUE,beamAngleAzUE,beamAngleElUE] = getDFTCodebook3D(Nrx_h,Nrx_v,8,8,0.5,0.5);
numBeamBS_opt = size(beamBS_opt,2);
numBeamUE_opt = size(beamUE_opt,2);
numBeamAzBS_opt = size(beamAzBS_opt,2);
numBeamElBS_opt = size(beamElBS_opt,2);
numBeamAzUE_opt = size(beamAzUE_opt,2);
numBeamElUE_opt = size(beamElUE_opt,2);
pR = ones(1,T);
SNR_linear = ones(1,T);
dr_Ex_hybrid = zeros(numItr,T);
PowerUE = zeros(numItr,T,numBeamUE_opt);
PowerBS = zeros(numItr,T,numBeamBS_opt);
freqChan = zeros(numPolarAngle*Nrx,numPolarAngle*Ntx,NumSC,T);
dr_Ex_3RF_los_hybrid_360 = zeros(1,T);

for t = 1:T
    for itr = 1:numItr
k=1;
 if itr<101
 load("Channel_matrix_test_traj1_los_speed_360_100itr_final.mat","Channel_2")
 freqChan(:,:,k,t) = Channel_2(itr,:,:,t,k);
 else
 load("Channel_matrix_test_traj1_los_speed_360_50itr_finale.mat","Channel_2")
 freqChan(:,:,k,t) = Channel_2(itr-100,:,:,t,k);
 end
[bpTable_Ex,rPower_Ex,SE_All,beam_Ex,~,~,~] = performBeamTraining_t_hybrid(t,k,beamBS_opt,1:numBeamBS_opt,beamUE_opt,1:numBeamUE_opt,freqChan,pR,SNR_linear,numPolarAngle,numBeams); 
FRF = zeros(Ntx,numBeams);
WRF = zeros(Nrx,numBeams);
for i =1:numBeams
    tb = beam_Ex(1,i);
    rb = beam_Ex(2,i);
    FRF(:,i) = beamBS_opt(:,tb);
    WRF(:,i) = beamUE_opt(:,rb);
end
H = freqChan(:,:,k,t);

[U, Sigma, V] = svd(H);
% Select RF precoder and combiner using singular vectors
FBB = hybrid_precoder_vector(FRF, WRF, H, Ns);
WBB = calculate_digital_combiner(WRF, FRF, H, Ns);

SNR_linear(1,t) = 1;
dr_Ex_hybrid(itr,t) = computeDataRate_Hybrid_t(t,1,FRF,WRF,freqChan,SNR_linear(1,t),FBB,WBB,Ns);
se = dr_Ex_hybrid(itr,t);
 end
 dr = mean(dr_Ex_hybrid(:,t),1)
dr_Ex_3RF_los_hybrid_360(1,t) = mean(dr_Ex_hybrid(:,t),1);
t = t
save("SE_Exhaustive_hybrid_3RF_los_traj1_speed_360.mat","dr_Ex_3RF_los_hybrid_360")
end
save("SE_Exhaustive_hybrid_3RF_los_traj1_speed_360.mat","dr_Ex_3RF_los_hybrid_360")