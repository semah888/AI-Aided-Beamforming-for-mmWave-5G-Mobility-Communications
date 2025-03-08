function [beamTable,rPower,SE_All,beamSelection,numTransmitBeam,numReceiveBeam,numBeamMeasure] = performBeamTraining_t(t,numsc,tbCodebook,tbIdxVec,rbCodebook,rbIdxVec,freqChan,pR,SNR_linear,numPolAng)
[~,~,~,T] = size(freqChan);
numTransmitBeam = size(tbCodebook,2);
numReceiveBeam = size(rbCodebook,2);
numBeamMeasure = numTransmitBeam*numReceiveBeam;
rPower = zeros(numTransmitBeam,numReceiveBeam);
beamTable = zeros(numBeamMeasure,2); 
SE_All = zeros(numTransmitBeam,numReceiveBeam);
x =0;
for sc = 1:numsc
       H = freqChan(:,:,sc,t); % Nr x Nt
       
       for tb = 1:numTransmitBeam % Search all beam pairs
            f = tbCodebook(:,tb); % Nt x 1
            f = repelem(f,numPolAng,1); %Nt*2 x 1
            for rb = 1:numReceiveBeam
                w = rbCodebook(:,rb); % Nr x 1
                w = repelem(w,numPolAng,1); %Nr*numPolAng x 1
              rData = sqrt(pR(t))*w'*H*f; % transmit symbol has unit power
                rPower(tb,rb) = rPower(tb,rb) + (1/numsc)*sum(abs(rData).^2); % sum over all RF chains
                SE_All(tb,rb) = SE_All(tb,rb) + (1/numsc)*abs(log2(det(eye(1)+SNR_linear(t)*w'*H*f*f'*H'*w)));
                x = x+1;
            end
       end
end
   power = rPower;
   for bp = 1:1 % obtain beampair table ordered by receive power
        [tB,rB] = find(power == max(max(power)));
        power(tB,rB) = -Inf;
        beamTable(bp,1) = tbIdxVec(tB);
        beamTable(bp,2) = rbIdxVec(rB);
    end
beamSelection = [tB rB]; % the first beampair
