clear;
close all;
% For repeatibility
stream = RandStream('mt19937ar','Seed',3);
RandStream.setGlobalStream(stream);
% System parameters
fc = 30e9; % carrier frequency, Hz
BW = 100e6; % bandwidth, Hz
speedUE = 1; % m/sec
updatePeriod = 0.5; % channel update period, second
timeEpoch = 30; % second
movingDis = speedUE*timeEpoch;
T = timeEpoch/updatePeriod; % number of locations
directionUE = 90; % horizontal
numSC = 64; % number of subcarriers 
nVar_dBm = -174+10*log10(BW)-10*log10(numSC); % noise power per subcarrier
nVar = 10^(nVar_dBm/10)*1e-3; % noise variance, watt
pT = 0.001; % transmit power, watt
disBSUE2D = 50; % BS-UE distance, meter
heightBS = 25; % meter
heightUE = 1.5; % meter
arraySizeBS = [8,8]; 
arraySizeUE = [8,8];
Nt = prod(arraySizeBS);
Nr = prod(arraySizeUE);
N = 23; % number of spatial clusters
velocityUE = speedUE*[cosd(directionUE);sind(directionUE);0];
velocityUE = repmat(velocityUE,1,T);
cellLayout = getCellLayout(disBSUE2D,heightBS,heightUE); % drop UE in the cell
locationBS = cellLayout.locationBS;
locationUE_t0 = cellLayout.locationUE;
locationUE_update = zeros(3,T);
locationUE_update(:,1) = locationUE_t0;
for t = 2:T
    locationUE_update(:,t) = velocityUE(:,t)*updatePeriod+locationUE_update(:,t-1);
end
% Calculate NLOS path loss
disBSUE3D_update = vecnorm(locationUE_update-repmat(locationBS,1,T),2); % 1 x T, meter
pathLoss = 32.4+20*log10(fc*1e-9)+30*log10(disBSUE3D_update); % dB, 1 x T
pR = pT./10.^(pathLoss./10); % average receive power, watt
SNR_linear = pR./nVar; % actual SNR per subcarrier
%CDL parameters
Delays = [0 0.3819 0.4025 0.5888 0.4610 0.5375 0.6708 0.5750 0.7618 0.575 0.7618 1.5375 1.8978 2.2242 2.1718 2.4942 2.5119 3.0582 4.0810 4.4579 4.7966 5.0066 5.3043 9.6586]; 
Power = 10.^([-13.4 0 -2.2 -4 -6 -8.2 -9.9 -10.5 -7.5 -15.9 -6.6 -16.7 -12.4 -15.2 -10.8 -11.3 -12.7 -16.2 -18.3 -18.9 -16.6 -19.9 -29.7]/10); % Watt
Cluster_AOD = [-178.1 -4.2 -4.2 -4.2 90.2 90.2 90.2 121.5 -81.7 158.4 -83 134.8 -153 -172 -129.9 -136 165.4 148.4 132.7 -118.6 -154.1 126.5 -56.2]; %degree
Cluster_ZOD = [50.2 93.2 93.2 93.2 122 122 122 150.2 55.2 26.4 126.4 171.6 151.4 157.2 47.2 40.4 43.3 161.8 10.8 16.7 171.1 22.7 144.9];
Cluster_AOA = [51.3 -152.7 -152.7 152.7 76.6 76.6 76.6 -1.8 -41.9 94.2 51.9 -115.9 26.6 76.6 -7 -23 -47.2 110.4 144.5 155.3 102 -151.8 55.2];
Cluster_ZOA = [125.4 91.3 91.3 91.3 94 94 94 47.1 56 30.1 58.8 26 49.2 143.1 117.4 122.7 123.2 32.6 27.2 15.2 146 150.7 156.1];
XPR = 10; %Cross polarization ratio (dB)
Xi = 10^(XPR/10); %Cross polarization ratio (W)
C_ASD = 5;
C_ZSD = 3;
C_ASA = 11;
C_ZSA = 3;
NRays = 20;
Alpha = [0.0447 -0.0447 0.1413 -0.1413 0.2492 -0.2492 0.3715 -0.3715 0.5129 -0.5129 0.6797 -0.6797 0.8844 -0.8844 1.1481 -1.1481 1.5195 -1.5195 2.1551 -2.1551];
cdl = nrCDLChannel;
cdl.DelayProfile = 'CDL-B';
cdl.AngleScaling = true;
cdl.CarrierFrequency = fc;
cdl.UTDirectionOfTravel = [0; directionUE];
txSize = [8 8 2 1 1];
cdl.TransmitAntennaArray.Size = [8 8 2 1 1];
cdl.ReceiveAntennaArray.Size = [8 8 2 1 1];
cdl.RandomStream = 'mt19937ar with seed';
lambda_v = 0.5;
lambda_h = 0.5;
dg_v = lambda_v*txSize(1); % lambda_v * M
dg_h = lambda_h*txSize(2); % lambda_h * N
cdlInfo = cdl.info
cdl.TransmitAntennaArray.ElementSpacing = [lambda_v lambda_h dg_v dg_h];
cdl.TransmitArrayOrientation = [cdlInfo.AnglesAoD(1) cdlInfo.AnglesZoD(1)-90 0]';
cdl.ReceiveArrayOrientation = [cdlInfo.AnglesAoA(1) cdlInfo.AnglesZoA(1)-90 0]';
figTx = displayChannel(cdl,'LinkEnd','Tx');
datacursormode on;
figRx = displayChannel(cdl,'LinkEnd','Rx');
datacursormode on;
cdl.ChannelFiltering = false;
SR = 15.36e6;
T = SR * 1e-3;
cdl.SampleRate = SR;
cdlinfo = info(cdl);
Nt = cdlinfo.NumInputSignals;
 
txWaveform = complex(randn(T,Nt),randn(T,Nt));
rxWaveform = cdl(txWaveform);
%txWaveform = complex(randn(T,Nt),randn(T,Nt));