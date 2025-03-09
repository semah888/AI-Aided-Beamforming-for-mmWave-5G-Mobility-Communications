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



% For repeatibility
stream = RandStream('mt19937ar','Seed',3);
RandStream.setGlobalStream(stream);

% System parameters
fc = 30e9; % carrier frequency, Hz
BW = 100e6; % bandwidth, Hz
speedUE = 120; % km/h
speedUE = speedUE/3.6;
updatePeriod = 0.00025; % s
T = 101; % number of locations
directionUE = -90; % horizontal
numSC = 100; % number of subcarriers 
nVar_dBm = -174+10*log10(BW)-10*log10(numSC); % noise power per subcarrier
nVar = 10^(nVar_dBm/10)*1e-3; % noise variance, watt
pT = 0.001; % transmit power, watt
xBS = 0;
yBS = 5;
disBSUE2D = 50; % BS-UE distance, meter
heightBS = 2.5; % meter
heightUE = 2.5; % meter
Ntx_h = 2;
Ntx_v = 2;
Nrx_h = 2;
Nrx_v = 2;
Polar_Slant_Angle = [45 -45];
[~,numPolAng] = size(Polar_Slant_Angle);
N = 20; % number of spatial clusters
load("MIMOchannel_details.mat","PropagChanModelVariable")
xUE = zeros(1,T);
for t = 1:101
    xUE(1,t) = sqrt(PropagChanModelVariable{t}.ScalingValues.DistanceBSMS^2-yBS^2);
end
yUE = 0;
% DFT codebooks
R_theta = (0.001:0.1:3.501);
[n,m] = size(R_theta);
dataRate11 = zeros(1,m);
dataRate12 = zeros(1,m);
dataRate13 = zeros(1,m);
dataRate14 = zeros(1,m);
NumberIteration1 = zeros(1,m);
NumberIteration2 = zeros(1,m);
load("dataRate1Correct.mat","dataRate1");
load("dr10Correct.mat","dataRate10");
load("dr20Correct.mat","dataRate20");
load("dr60Correct.mat","dataRate60");
dataRate1Opt = dataRate1(1,5)*ones(1,m);
dataRate2Opt = dataRate10(1,5)*ones(1,m);
dataRate3Opt = dataRate20(1,5)*ones(1,m);
dataRate4Opt = dataRate60(1,5)*ones(1,m);
NumberIteration_opt = (72*181)^2*ones(1,m);
i = 0;
k=0;
load('MIMOchannel_details.mat','channel_matrix');

for r = 0.001:0.1:3.501
    r = r
    i = i+1;
[beamBS1,beamAzBS1,beamElBS1,~,~] = getProposedCodebookBS3D(Ntx_v,Ntx_h,xBS,yBS,xUE(1,1),yUE,r,r);
[beamUE1,beamAzUE1,beamElUE1,~,~] = getProposedCodebookUE3D(Ntx_v,Ntx_h,xUE(1,1),yUE,xBS,yBS,r,r);
[beamBS2,beamAzBS2,beamElBS2,~,~] = getProposedCodebookBS3D(Ntx_v,Ntx_h,xBS,yBS,xUE(1,10),yUE,r,r);
[beamUE2,beamAzUE2,beamElUE2,~,~] = getProposedCodebookUE3D(Ntx_v,Ntx_h,xUE(1,10),yUE,xBS,yBS,r,r);
[beamBS3,beamAzBS3,beamElBS3,~,~] = getProposedCodebookBS3D(Ntx_v,Ntx_h,xBS,yBS,xUE(1,20),yUE,r,r);
[beamUE3,beamAzUE3,beamElUE3,~,~] = getProposedCodebookUE3D(Ntx_v,Ntx_h,xUE(1,20),yUE,xBS,yBS,r,r);
[beamBS4,beamAzBS4,beamElBS4,~,~] = getProposedCodebookBS3D(Ntx_v,Ntx_h,xBS,yBS,xUE(1,60),yUE,r,r);
[beamUE4,beamAzUE4,beamElUE4,~,~] = getProposedCodebookUE3D(Ntx_v,Ntx_h,xUE(1,60),yUE,xBS,yBS,r,r);
numBeamBS1 = size(beamBS1,2);
numBeamUE1 = size(beamUE1,2);
NumberIteration1(1,i) = numBeamUE1*numBeamBS1;


numBeamBS2 = size(beamBS2,2);
numBeamUE2 = size(beamUE2,2);
NumberIteration2(1,i) = numBeamUE2*numBeamBS2;


numBeamBS3 = size(beamBS3,2);
numBeamUE3 = size(beamUE3,2);


numBeamBS4 = size(beamBS4,2);
numBeamUE4 = size(beamUE4,2);

%% Channel matrices

% Load channel matrices for 1 Monte-Carlo run

% Frequency-domain channel
Nr = Nrx_h*Nrx_v*numPolAng;
Nt = Ntx_h*Ntx_v*numPolAng;
freqChan = zeros(Nr,Nt,numSC,T);
for t = 1:T
    for sc = 1:numSC
        for nr = 1:Nr
            for nt = 1:Nt
                freqChan(nr,nt,sc,t) = channel_matrix(1,nr,nt,t,sc);
            end
        end
    end
end

% Calculate NLOS path loss

pR = ones(1,T);
%best beam pair at time t = 1
[bpTable_Ex1,rPower_Ex1,beam_Ex1,~,~,~] = performBeamTraining_t(1,1,beamBS1,1:numBeamBS1,beamUE1,1:numBeamUE1,freqChan,pR,numPolAng); % 1 x 2 x T
[bpTable_Ex2,rPower_Ex2,beam_Ex2,~,~,~] = performBeamTraining_t(10,1,beamBS2,1:numBeamBS2,beamUE2,1:numBeamUE2,freqChan,pR,numPolAng); % 1 x 2 x T
[bpTable_Ex3,rPower_Ex3,beam_Ex3,~,~,~] = performBeamTraining_t(20,1,beamBS3,1:numBeamBS3,beamUE3,1:numBeamUE3,freqChan,pR,numPolAng); % 1 x 2 x T
[bpTable_Ex4,rPower_Ex4,beam_Ex4,~,~,~] = performBeamTraining_t(60,1,beamBS4,1:numBeamBS4,beamUE4,1:numBeamUE4,freqChan,pR,numPolAng); % 1 x 2 x T
SNR_linear = zeros(1,T);
SNR = 0; %dB
    SNR_linear(1) = 10^(SNR/10);
    SNR_linear(10) = 10^(SNR/10);
    SNR_linear(20) = 10^(SNR/10);
    SNR_linear(60) = 10^(SNR/10);
    dataRate11(1,i) = computeDataRate_t(1,1,beam_Ex1,freqChan,beamBS1,beamUE1,1,SNR_linear(1),numPolAng)
    dataRate12(1,i) = computeDataRate_t(10,1,beam_Ex2,freqChan,beamBS2,beamUE2,1,SNR_linear(10),numPolAng)
    dataRate13(1,i) = computeDataRate_t(20,1,beam_Ex3,freqChan,beamBS3,beamUE3,1,SNR_linear(20),numPolAng)
    dataRate14(1,i) = computeDataRate_t(60,1,beam_Ex4,freqChan,beamBS4,beamUE4,1,SNR_linear(60),numPolAng)
 end
load("dataRate1.mat","dataRate1")
figure();
plot(R_theta,dataRate11,'b-+','LineWidth',2);hold on
plot(R_theta,dataRate1Opt,'k-+','LineWidth',2);hold on
plot(R_theta,dataRate12,'r-o','LineWidth',2);hold on
plot(R_theta,dataRate2Opt,'k-o','LineWidth',2);hold on
figure();
plot(R_theta,dataRate13,'b->','LineWidth',2);hold on
plot(R_theta,dataRate3Opt,'k->','LineWidth',2);hold on
plot(R_theta,dataRate14,'r-<','LineWidth',2);hold on
plot(R_theta,dataRate4Opt,'k-<','LineWidth',2);hold on
figure();
plot(R_theta,NumberIteration1,'g-','LineWidth',2);hold on
%plot(R_theta,NumberIteration2,'b-.','LineWidth',2);
plot(R_theta,NumberIteration_opt,'k:','LineWidth',2);hold on
%yticks([0*1e4 0.5*1e4 1*1e4 1.5*1e4 2*1e4 2.5*1e4 (72*181)^2 (72*181)^2+1*1e4])
%best beam pair at time t = 20
%[bpTable_Ex20,rPower_Ex20,beam_Ex20,~,~,~] = performBeamTraining_t(20,1,beamBS_opt,1:numBeamBS_opt,beamUE_opt,1:numBeamUE_opt,freqChan,pR,numPolAng) % 1 x 2 x T
%dataRate20 = zeros(1,7);
%i = 1;
%for SNR = -30:5:0
%    SNR_linear(20) = 10^(SNR/10);
%    dataRate20(1,i) = computeDataRate_t(20,1,beam_Ex20,freqChan,beamBS_opt,beamUE_opt,1,SNR_linear,numPolAng)
%    i = i+1;
%end
%plot(SNR_vect,dataRate20,'b->','LineWidth',2);hold on;
%best beam pair at time t = 60
%[bpTable_Ex60,rPower_Ex60,beam_Ex60,~,~,~] = performBeamTraining_t(60,1,beamBS_opt,1:numBeamBS_opt,beamUE_opt,1:numBeamUE_opt,freqChan,pR,numPolAng) % 1 x 2 x T
%dataRate60 = zeros(1,7);
%i = 1;
%for SNR = -30:5:0
 %   SNR_linear(60) = 10^(SNR/10);
 %   dataRate60(1,i) = computeDataRate_t(60,1,beam_Ex60,freqChan,beamBS_opt,beamUE_opt,1,SNR_linear,numPolAng)
%    i = i+1;
%end
%plot(SNR_vect,dataRate60,'g-o','LineWidth',2);hold off;
%legend('Optimal t = 1','Optimal t = 20','Optimal t = 60');

