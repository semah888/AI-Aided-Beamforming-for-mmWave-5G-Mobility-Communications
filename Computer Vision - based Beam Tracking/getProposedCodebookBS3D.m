function [beam3D,beamAz,beamEl,beamAngle,beamAngleAz,beamAngleEl] = getProposedCodebookBS3D(numElementV,numElementH,xBS,yBS,xUE_t,yUE_t,r_phi_t,r_theta_t,AODLOS_ref_t,ZODLOS_ref_t)
%azimuth angles 

d_t = sqrt((xBS-xUE_t)^2+(yBS-yUE_t)^2);
phi_t = atan((yBS-yUE_t)/(xBS-xUE_t));
phi_t = AODLOS_ref_t;
Eps_phi_t = asin(r_phi_t/abs(yUE_t-yBS));
pa_phi = 1*pi/(180); % 1°
beamAngleAz = (phi_t-Eps_phi_t:pa_phi:phi_t+Eps_phi_t);
[~,numBeamH] = size(beamAngleAz)

% elevation angles
theta_t = 0;
Eps_theta_t = asin(r_theta_t/abs(yUE_t-yBS));
pa_theta = 1*pi/(180); % 1°
beamAngleEl = (theta_t-Eps_theta_t:pa_theta:theta_t+Eps_theta_t);
[~,numBeamV] = size(beamAngleEl)
numAngle = numBeamV*numBeamH; % total number of resolvable angles
virtAngEl = (0:numBeamV-1)/numBeamV;
virtAngAz = (0:numBeamH-1)/numBeamH;
beamEl = exp(1j*(0:numElementV-1).'*beamAngleEl); 
beamAz= exp(1j*(0:numElementH-1).'*beamAngleAz);
beam3D = zeros(numElementH*numElementV,numAngle);
beamAngleEl = rad2deg(beamAngleEl);
beamAngleAz = rad2deg(beamAngleAz);
beamAngle = zeros(2,numAngle);
idx = 1;
for nv = 1:numel(virtAngEl)
    for nh = 1:numel(virtAngAz)
        beam2D = flip(beamEl(:,nv)).*beamAz(:,nh).'; % numV x numH
        beam3D(:,idx) = beam2D(:); % numElement x 1 [numV x 1;numV x 1;...]
        beamAngle(:,idx) = [beamAngleEl(nv);beamAngleAz(nh)];
        idx = idx+1;
    end
end