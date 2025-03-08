function [F_theta,F_phi] = RadiationPatternGenerator1(theta,phi,Polar_Slant_Angle)
%% Use model2 page 24
[N,m] = size(theta);
F_theta = zeros(N,1);
F_phi= zeros(N,1);
theta_3dB = 65;
phi_3dB = 65;
A_max = 30;
SLA_v = 30;
for n=1:N
A_vertical = -min(12*((theta(n,1)-90)/theta_3dB)^2,SLA_v);
A_horizontal = -min(12*(phi(n,1)/phi_3dB)^2,A_max);
A = -min(-(A_vertical+A_horizontal),A_max);
F_theta(n,1) = sqrt(A)*cosd(Polar_Slant_Angle);
F_phi(n,1) = sqrt(A)*sind(Polar_Slant_Angle);
end
end