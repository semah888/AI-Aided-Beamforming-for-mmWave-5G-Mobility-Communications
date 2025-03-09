% This script is to implement beam training strategies proposed in 
% Section III-B in [1] and generate similar simulation results in Fig. 5. 
% The spatially-consistent channel matrices are generated according to the 
% 3GPP spatial consistency Prodecure A in [2]. The main purpose is to show 
% the beam training process, including the construction of codebooks, the 
% acquisition of candidate beam pairs and local beam search. Hence the 
% implementation of the 3GPP channel model is omitted. Readers can 
% implement the Prodecure A themselves, associated with the 3GPP CDL 
% channel model provided by MATLAB. The channel matrices of 1 monte-carlo 
% run is saved as .mat file and loaded in the simulation. The exact 
% simulation results shown in [1] are averaged over 500 monte-carlo 
% simulations. The beam training process is assumed to be noise-free. 
%
% Author: Narengerile
% Date: 13 July 2020
% 
% [1] Narengerile,F.Alsaleem,J.Thomspon,T.Ratnarajah,"Low-Complexity Beam 
% Training for Tracking Spatially Consistent Millimeter Wave Channels",
% PIMRC, 2020.
% [2] 3GPP TR 38.901,"Study on channel model for frequencies from 0.5 to 
% 100 GHz," 2017.
%% Simulation parameters
clear;
close all;
% For repeatibility
stream = RandStream('mt19937ar','Seed',3);
RandStream.setGlobalStream(stream);
% System parameters
fc = 30e9; % carrier frequency, Hz
BW = 60e6; % bandwidth, Hz
speedUE = 120;
speedUE = speedUE*1000/3600;
updatePeriod = 0.25*1e-3; % channel update period, second
NumberGrids = 101;
T = NumberGrids;
directionUE = 180; % horizontal
numSC = 100; % number of subcarriers 
numRFBS = 2; % number of BS RF chains
numRFUE = 2; % number of UE RF chains
nVar_dBm = -174+10*log10(BW)-10*log10(numSC); % noise power per subcarrier
nVar = 10^(nVar_dBm/10)*1e-3; % noise variance, watt
pT = 0.001; % transmit power, watt
disBSUE2D = 1; % BS-UE distance, meter
heightBS = 25; % meter
heightUE = 1.5; % meter
arraySizeBS = [2,2]; 
arraySizeUE = [2,2];
Nt = prod(arraySizeBS);
Nr = prod(arraySizeUE);
Polar = 45; %Cross_polarization
if Polar == 0
 Polar_Slant_Angle = [0];
elseif Polar == 90
 Polar_Slant_Angle = [90];
elseif Polar == 45
 Polar_Slant_Angle = [45 -45];
end
numPolAng= 2;
numItr = 100;
velocityUE = speedUE*[-1;0;0];
velocityUE = repmat(velocityUE,1,T);
cellLayout = getCellLayout(disBSUE2D,heightBS,heightUE); % drop UE in the cell
locationBS = [0;5;2.5];
locationUE_t0 = [550;0;2.5];
locationUE_update = zeros(3,T);
locationUE_update(:,1) = locationUE_t0;
for t = 2:T
 locationUE_update(:,t) = velocityUE(:,t)*updatePeriod+locationUE_update(:,t-1)
end
% UE trajectory
figure(1);
dpUE = locationUE_update(1:2,end)-locationUE_update(1:2,1);
p1 = plot(locationBS(1),locationBS(2),'rpentagram','MarkerFaceColor','r','MarkerSize',20);hold on;
p2 = plot(locationUE_update(1,1),locationUE_update(2,1));hold on;
p3 = plot(locationUE_update(1,end),locationUE_update(2,end));hold on;
quiver(locationUE_update(1,1),locationUE_update(2,1),dpUE(1),dpUE(2),0,'MaxHeadSize',2);hold on;
dp1 = [0,20];dp2 = [20*cosd(30);-20*sind(30)];dp3 = [-20*cosd(30);-20*sind(30)];
%xlim([-100,100]);ylim([-100,100]);
grid on;
N = 20; % number of spatial clusters
disBSUE3D_update = vecnorm(locationUE_update-repmat(locationBS,1,T),2); % 1 x T, meter
pathLoss = 32.4+20*log10(fc*1e-9)+30*log10(disBSUE3D_update); % dB, 1 x T
pR = ones(1,T); % average receive power, watt
SNR_linear = pR./nVar; % actual SNR per subcarrier
SNR_linear = ones(1,T);
%% Beam training process - local beam search
% Exhaustive search 
load("MIMOchannel_details.mat","PropagChanModelVariable");
 d_2D_ref = zeros(1,101);
 AODLOS_ref = zeros(1,101);
 ZODLOS_ref = zeros(1,101);
 for t = 1:101
 d_2D_ref(1,t) = PropagChanModelVariable{t}.ScalingValues.DistanceBSMS
 ZODLOS_ref(1,t) = (PropagChanModelVariable{t}.ScalingValues.ZLOSBSToMS);
 AODLOS_ref(1,t) = (PropagChanModelVariable{t}.ScalingValues.ALOSBSToMS);
 end
Az_anglesBSOpt = zeros(1,3);
El_anglesBSOpt = zeros(1,3);
Az_anglesUEOpt = zeros(1,3);
El_anglesUEOpt = zeros(1,3);
Az_anglesBSProp = zeros(1,3);
El_anglesBSProp = zeros(1,3);
Az_anglesUEProp = zeros(1,3);
El_anglesUEProp = zeros(1,3);
load("Chan_3GPP_Merce_Model_50Itrs.mat",'Channel');
load("Chan_3GPP_Merce_Model_50Itrs_V2.mat","Channel_2");
Channel_matrix = cat(1,Channel,Channel_2);
dr_Prop = zeros(numItr,T,7);
dr_Ex = zeros(numItr,T);
dr_Lc = zeros(numItr,T);
dataRate_prop = zeros(numItr,2,7);
for itr = 1:numItr
 iteration = itr
freqChan = zeros(8,8,numSC,T);
for t = 1:T
 for sc = 1:numSC
 freqChan(:,:,sc,t) = Channel_matrix(itr,:,:,t,sc); %channel_matrix(1,:,:,t,sc);
 end
end
instant = 0;
for t = [20 66]
 instant = instant+1
% DFT codebooks
[beamBSProp,beamAzBSProp,beamElBSProp,beamAngleBSProp,beamAngleAzBSProp,beamAngleElBSProp] = getProposedCodebookBS3D(arraySizeBS(1),arraySizeBS(2),locationBS(1,1),locationBS(2,1),locationUE_update(1,t),locationUE_update(2,t),0.01,0.01,AODLOS_ref(1,t),ZODLOS_ref(1,t));
[beamUEProp,beamAzUEProp,beamElUEProp,beamAngleUEProp,beamAngleAzUEProp,beamAngleElUEProp] = getProposedCodebookUE3D(arraySizeBS(1),arraySizeBS(2),locationUE_update(1,t),locationUE_update(2,t),locationBS(1,1),locationBS(2,1),0.01,0.01,AODLOS_ref(1,t),ZODLOS_ref(1,t));
numAngAzUE = size(beamAzUEProp,2);
numAngElUE = size(beamElUEProp,2);
numBeamBSProp = size(beamBSProp,2);
numBeamUEProp = size(beamUEProp,2);
numBeamAzBSProp = size(beamAzBSProp,2);
numBeamElBSProp = size(beamElBSProp,2);
numBeamAzUEProp = size(beamAzUEProp,2);
numBeamElUEProp = size(beamElUEProp,2);
%% Channel matrices
% Load channel matrices for 1 Monte-Carlo run
% Calculate NLOS path loss
[bpTable_Ex,rPower_Ex,beam_Ex,~,~,~] = performBeamTraining_t(t,1,beamBSProp,1:numBeamBSProp,beamUEProp,1:numBeamUEProp,freqChan,pR,numPolAng); % 1 x 2 x T
i = 1;
for SNR = -30:5:0
 SNR_linear(t) = 10^(SNR/10);
 dataRate_prop(itr,instant,i) = computeDataRate_t(t,1,beam_Ex,freqChan,beamBSProp,beamUEProp,1,SNR_linear,numPolAng);
 i = i+1;
end
end
end
dataRate_Prop_mean = (1/numItr)*sum(dataRate_prop,1);
%save("SE_Prop_03.mat","dataRate_Prop_mean_03");
save("SE_Prop_20_60.mat","dataRate_Prop_mean");
figure()
%plot(dataRate_Prop_mean_03(1,:,end),'b','LineWidth',2);hold on;
snr_vect = -30:5:0;
dr_prop_20 = zeros(1,7);
dr_prop_66 = zeros(1,7);
for i = 1:7
 dr_prop_20(1,i)=dataRate_Prop_mean(1,1,i);
end
for i = 1:7
 dr_prop_66(1,i)=dataRate_Prop_mean(1,2,i);
end
plot(snr_vect,dr_prop_20,'r:>','LineWidth',2);hold on;
plot(snr_vect,dr_prop_66,'b:o','LineWidth',2);hold on;
load("SE_20_50Itr.mat","dataRate20_mean");
plot(snr_vect,dataRate20_mean,'k->');hold on;
load("SE_66_50Itr.mat","dataRate66_mean");
plot(snr_vect,dataRate66_mean,'k-o');
% create a new pair of axes inside current figure
%best beam pair at time t = 20
%[bpTable_Ex20,rPower_Ex20,beam_Ex20,~,~,~] = performBeamTraining_t(20,1,beamBS_opt,1:numBeamBS_opt,beamUE_opt,1:numBeamUE_opt,freqChan,pR,numPolAng) % 1 x 2 x T
%dataRate20 = zeros(1,7);
%i = 1;
%for SNR = -30:5:0
% SNR_linear(20) = 10^(SNR/10);
% dataRate20(1,i) = computeDataRate_t(20,1,beam_Ex20,freqChan,beamBS_opt,beamUE_opt,1,SNR_linear,numPolAng)
% i = i+1;
%end
%plot(SNR_vect,dataRate20,'b->','LineWidth',2);hold on;
%best beam pair at time t = 60
%[bpTable_Ex60,rPower_Ex60,beam_Ex60,~,~,~] = performBeamTraining_t(60,1,beamBS_opt,1:numBeamBS_opt,beamUE_opt,1:numBeamUE_opt,freqChan,pR,numPolAng) % 1 x 2 x T
%dataRate60 = zeros(1,7);
%i = 1;
%for SNR = -30:5:0
 % SNR_linear(60) = 10^(SNR/10);
 % dataRate60(1,i) = computeDataRate_t(60,1,beam_Ex60,freqChan,beamBS_opt,beamUE_opt,1,SNR_linear,numPolAng)
% i = i+1;
%end
%plot(SNR_vect,dataRate60,'g-o','LineWidth',2);hold off;
%legend('Optimal t = 1','Optimal t = 20','Optimal t = 60');