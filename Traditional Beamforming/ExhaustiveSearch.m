function [SelectBeamTx,SelectBeamRx] = ExhaustiveSearch(tbCodebook,rbCodebook,H)
numTransmitBeam = size(tbCodebook,2);
numReceiveBeam = size(rbCodebook,2);
rPower = zeros(numTransmitBeam,numReceiveBeam);
        for tb = 1:numTransmitBeam % Search all beam pairs
            f = tbCodebook(:,tb); % Nt x 1
            for rb = 1:numReceiveBeam
                w = rbCodebook(:,rb); % Nr x 1
                rData = w'*H*f; % transmit symbol has unit power
                rPower(tb,rb) = abs(rData); % sum over all RF chains
                power = rPower(tb,rb)
            end
        end
      [tB,rB] = find(rPower == max(max(rPower)));
   SelectBeamTx = tbCodebook(:,tB);
   SelectBeamRx = rbCodebook(:,rB);
end
