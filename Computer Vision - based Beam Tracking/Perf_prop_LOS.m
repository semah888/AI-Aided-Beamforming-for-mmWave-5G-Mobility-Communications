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
[beamBS1,beamAzBS1,beamElBS1,~,~] = getProposedCodebookBS3D(Ntx_v,Ntx_h,xBS,yBS,xUE(1,20),yUE,0.5,0.5)
[beamUE1,beamAzUE1,beamElUE1,~,~] = getProposedCodebookUE3D(Ntx_v,Ntx_h,xUE(1,20),yUE,xBS,yBS,0.5,0.5)
numBeamBS1 = size(beamBS1,2)
numBeamUE1 = size(beamUE1,2)
numBeamAzBS1 = size(beamAzBS1,2);
numBeamElBS1 = size(beamElBS1,2);
numBeamAzUE1 = size(beamAzUE1,2);
numBeamElUE1 = size(beamElUE1,2);
[beamBS2,beamAzBS2,beamElBS2,~,~] = getProposedCodebookBS3D(Ntx_v,Ntx_h,xBS,yBS,xUE(1,60),yUE,0.5,0.5);
[beamUE2,beamAzUE2,beamElUE2,~,~] = getProposedCodebookUE3D(Ntx_v,Ntx_h,xUE(1,60),yUE,xBS,yBS,0.5,0.5);
numBeamBS2 = size(beamBS2,2)
numBeamUE2 = size(beamUE2,2)
numBeamAzBS2 = size(beamAzBS2,2);
numBeamElBS2 = size(beamElBS2,2);
numBeamAzUE2 = size(beamAzUE2,2);
numBeamElUE2 = size(beamElUE2,2);

%% Channel matrices

% Load channel matrices for 1 Monte-Carlo run
load('MIMOchannel_details.mat','channel_matrix');

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
[bpTable_Ex1,rPower_Ex1,beam_Ex1,~,~,~] = performBeamTraining_t(20,1,beamBS1,1:numBeamBS1,beamUE1,1:numBeamUE1,freqChan,pR,numPolAng); % 1 x 2 x T
dataRate1_Pro = zeros(1,7);
dataRate2_Pro = zeros(1,7);
SNR_vect = zeros(1,7);
i = 1;
SNR_vect = zeros(1,7);
for SNR = -30:5:0
    SNR_vect(i) = SNR;
    SNR_linear = 10^(SNR/10);
    dataRate1_Pro(1,i) = computeDataRate_t(20,1,beam_Ex1,freqChan,beamBS1,beamUE1,1,SNR_linear,numPolAng)
    i = i+1;
end
[bpTable_Ex2,rPower_Ex2,beam_Ex2,~,~,~] = performBeamTraining_t(60,1,beamBS2,1:numBeamBS2,beamUE2,1:numBeamUE2,freqChan,pR,numPolAng); % 1 x 2 x T
SNR_vect = zeros(1,7);
i = 1;
for SNR = -30:5:0
    SNR_vect(i) = SNR;
    SNR_linear = 10^(SNR/10);
    dataRate2_Pro(1,i) = computeDataRate_t(60,1,beam_Ex2,freqChan,beamBS2,beamUE2,1,SNR_linear,numPolAng)
    i = i+1;
end

load("dr20.mat","dataRate20")
for i=1:7
   SNR =  SNR_vect(i);
    SNR_linear = 10^(SNR/10);
    BG = (exp(dataRate20(1,i))-1)/SNR_linear;
    dataRate20(1,i) = abs(log2(1+SNR_linear*BG/norm(freqChan(:,:,1,20))^2))
end
save("dr20Correct.mat","dataRate20")
load("dr_60_opt.mat","dataRate60")
for i=1:7
   SNR =  SNR_vect(i);
    SNR_linear = 10^(SNR/10);
    BG = (exp(dataRate60(1,i))-1)/SNR_linear;
    dataRate60(1,i) = abs(log2(1+SNR_linear*BG/norm(freqChan(:,:,1,60))^2))
end
save("dr60Correct.mat","dataRate60")
plot(SNR_vect,dataRate20,'k:+','LineWidth',2);hold on;
plot(SNR_vect,dataRate1_Pro,'b-+','LineWidth',2);hold on;
plot(SNR_vect,dataRate60,'k:o','LineWidth',2);hold on;
plot(SNR_vect,dataRate2_Pro,'r-o','LineWidth',2);hold off;
% create a new pair of axes inside current figure

