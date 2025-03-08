function [beamTable,rPower,beamSelection,numTransmitBeam,numReceiveBeam,numBeamMeasure] = performBeamTraining(tbCodebook,tbIdxVec,rbCodebook,rbIdxVec,freqChan,pR,numPolAng)
T = 1001;
numSC = 1;
numTransmitBeam = size(tbCodebook,2);
numReceiveBeam = size(rbCodebook,2);
numBeamMeasure = numTransmitBeam*numReceiveBeam;
rPower = zeros(numTransmitBeam,numReceiveBeam,T);
beamTable = zeros(numBeamMeasure,2,T); 
for t = 1:T
    rPowerTemp = zeros(numTransmitBeam,numReceiveBeam,numSC);
    for sc = 1:numSC
        H = zeros(4,4);
        for nr = 1:Nr
            for nt = 1:Nt
                H(nr,nt) = freqChan(nr,nt,t,sc); % Nr x Nt
            end
        end
        for tb = 1:numTransmitBeam % Search all beam pairs
            f = tbCodebook(:,tb); % Nt x 1
          n= size(f);
            for rb = 1:numReceiveBeam
                w = rbCodebook(:,rb); % Nr x 1
                m = size(w);
                rData = sqrt(pR(1,t))*w'*H*f; % transmit symbol has unit power
                rPowerTemp(tb,rb,sc) = sum(abs(rData).^2); % sum over all RF chains
            end
        end
    end
    rPower(:,:,t) = mean(rPowerTemp,3); % average over subcarriers
    
    power = rPower(:,:,t);
    for bp = 1:1 % obtain beampair table ordered by receive power
        [tB,rB] = find(power == max(max(power)));
        power(tB,rB) = -Inf;
    end
end
beamSelection = beamTable(1,:,:); % the first beampair
