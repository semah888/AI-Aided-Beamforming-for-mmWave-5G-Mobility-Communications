function [beamTable,rPower,SE_All,beamSelection,numTransmitBeam,numReceiveBeam,numBeamMeasure] = performBeamTraining_t_hybrid(t,numsc,tbCodebook,tbIdxVec,rbCodebook,rbIdxVec,freqChan,pR,SNR_linear,numPolAng,numBeams)
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
   beamSelection = zeros(2,numBeams);
   for bp = 1:1 % obtain beampair table ordered by receive power
       % Number of max indices to find
n = numBeams;

% Flatten the matrix into a single vector and get the indices
[sortedValues, linearIndices] = sort(power(:), 'descend');

% Get the first n indices
topNIndices = linearIndices(1:n);

% Convert the linear indices back to row and column indices
[row, col] = ind2sub(size(power), topNIndices);
       % [tB,rB] = find(power == max(max(power)));
      %  power(tB,rB) = -1000;
       % beamTable(bp,1) = tbIdxVec(tB);
       % beamTable(bp,2) = rbIdxVec(rB);
       for i = 1:n
           beamSelection(1,i) = row(i);
           beamSelection(2,i) = col(i);
       end
        
    end
 % the first beampair
