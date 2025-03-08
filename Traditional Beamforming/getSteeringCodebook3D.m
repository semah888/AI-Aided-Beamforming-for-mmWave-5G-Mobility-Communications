function [beam3D,beamAz,beamEl,beamAngleAz,beamAngleEl] = getSteeringCodebook3D(numElementV,numElementH,numBeamH,numBeamV,eleSpacingV,eleSpacingH)
numAngle = numBeamV*numBeamH; % total number of resolvable angles
% Vertical - resolved elevation angles
virtAngleEl = (0:numBeamV-1)./numBeamV;
av = unique(real(asind(virtAngleEl./eleSpacingV)),'stable');
bv = flip(-av(1:end-1));
beamAngleEl = unique([av bv],'stable');
%beamEl = 1/sqrt(numElementV)*exp(1j*2*pi*((0:numElementV-1)-numElementV/2).'*virtAngleEl);
%beamEl = (1/sqrt(numElementV))*(1j).^fix((0:numElementV-1)'.*mod(((1:numElementV)+numElementV/2),numElementV)./(numElementV/4));
beamEl = (1/sqrt(numElementV))*exp(1j*pi*(0:numElementV-1).'*cos(-pi/(2*numBeamV-1)+(2*(0:numBeamV-1)-1)*pi/(2*numBeamV))); 
% Horizontal - resolved azimuth angles
virtAngleAz = (0:numBeamH-1)./numBeamH;
ah = unique(real(asind(virtAngleAz./eleSpacingH)),'stable');
bh = flip(-ah(1:end-1));
beamAngleAz = unique([ah bh],'stable');
%beamAz = 1/sqrt(numElementH)*exp(1j*2*pi*((0:numElementH-1)-numElementH/2).'*virtAngleAz); 
%beamAz = (1/sqrt(numElementH))*(1j).^fix((0:numElementH-1)'.*mod(((1:numElementH)+numElementH/2),numElementH)./(numElementH/4));
beamAz= (1/sqrt(numElementH))*exp(1j*pi*(0:numElementH-1).'*sin(pi*(2*(0:numBeamH-1)-1)./numBeamH).*sin(-pi/(2*numBeamH)+(2*(0:numBeamH-1)-1)*pi/(2*numBeamH)));
beam3D = zeros(numElementH*numElementV,numAngle);
beamAngle = zeros(2,numAngle);
idx = 1;
x =numel(virtAngleEl)
for nv = 1:numel(virtAngleEl)
    for nh = 1:numel(beamAngleAz)
        beam2D = flip(beamEl(:,nv)).*beamAz(:,nh).'; % numV x numH
        beam3D(:,idx) = beam2D(:); % numElement x 1 [numV x 1;numV x 1;...]
        beamAngle(:,idx) = [beamAngleEl(nv);beamAngleAz(nh)];
        idx = idx+1;
    end
end
