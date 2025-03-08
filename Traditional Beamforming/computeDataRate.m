function dataRate = computeDataRate(beamSelection,freqChan,beamBS,beamUE,RF,SNR_linear,numPolAng)
T = 1001;
numSC = 1;
dataRateTemp = zeros(numSC,T);
for t = 1:T  
    f = beamBS(:,beamSelection(:,1,t).');
    %f = repelem(f,numPolAng,1);
    w = beamUE(:,beamSelection(:,2,t).');
    %w = repelem(w,numPolAng,1);
    for sc = 1:numSC
        H = freqChan(:,:,t,sc);
        dataRateTemp(sc,t) = abs(log2(det(eye(RF))+(SNR_linear(t)/RF)*w'*H*f*f'*H'*w));
    end
end
dataRate = mean(dataRateTemp,1); % averaged over subcarriers