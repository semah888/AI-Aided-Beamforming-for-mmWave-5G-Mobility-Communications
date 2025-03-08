function [beam3D,beamAz,beamEl,beamAngle,beamAngleAz,beamAngleEl] = getDFTCodebook3D_Opt(numElementV,numElementH,numBeamAZ,numBeamEl,eleSpacingV,eleSpacingH)

numAngle = numBeamAZ*numBeamEl; % total number of resolvable angles

% Vertical - resolved elevation angles
virtAngleEl = (0:numBeamEl-1)./numBeamEl;
av = unique(real(asind(-virtAngleEl./eleSpacingV)),'stable');
bv = flip(-av(1:end-1));
beamAngleEl = 2*pi*virtAngleEl;
beamAngleEldeg = rad2deg(beamAngleEl)
beamEl = exp(1j*2*pi*((0:numElementV-1)).'*virtAngleEl);

% Horizontal - resolved azimuth angles
virtAngleAz = (0:numBeamAZ-1)./numBeamAZ;
ah = unique(real(asind(-virtAngleAz./eleSpacingH)),'stable');
bh = flip(-ah(1:end-1));
beamAngleAz = pi/2-pi*virtAngleAz;
beamAngleAzdeg = rad2deg(beamAngleAz)
beamAz = exp(1j*2*pi*((0:numElementH-1)).'*virtAngleAz); 

beam3D = zeros(numElementH*numElementV,numAngle);
beamAngle = zeros(2,numAngle);
idx = 1;
x = numel(virtAngleEl)
y = numel(virtAngleAz)
for nv = 1:numel(virtAngleEl)
    for nh = 1:numel(virtAngleAz)
        beam2D = flip(beamEl(:,nv)).*beamAz(:,nh).'; % numV x numH
        beam3D(:,idx) = beam2D(:); % numElement x 1 [numV x 1;numV x 1;...]
        beamAngle(:,idx) = [beamAngleEl(nv);beamAngleAz(nh)];
        idx = idx+1;
    end
end