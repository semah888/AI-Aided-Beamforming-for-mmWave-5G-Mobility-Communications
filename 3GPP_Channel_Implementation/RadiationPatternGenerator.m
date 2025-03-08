function [F_theta,F_phi] = RadiationPatternGenerator(theta,phi,Polar_Slant_Angle)
%% Use model2 page 24
[T,N,M] = size(theta);
[p,numAng] = size(Polar_Slant_Angle);
F_theta = zeros(T,N,M,numAng);
F_phi= zeros(T,N,M,numAng);
theta_3dB = 65;
phi_3dB = 65;
A_max = 30;
SLA_v = 30;
bita = 0; %Same antenna hauteur at BS and UT
for t = 1:T
for n=1:N
    for m =1:M
        psi = angle(sind((theta(t,n,m)))*cos(bita)-cosd((phi(t,n,m)))*cosd((theta(t,n,m)))*sin(bita)+1j*sind((phi(t,n,m)))*sin(bita));
        theta_prim = acos(cosd((phi(t,n,m)))*sind((theta(t,n,m)))*sin(bita)+cos(bita)*cosd((theta(t,n,m))));
        theta_prim = rad2deg(theta_prim);
        if 180<=theta_prim
            theta_prim = 360-theta_prim;
        end
        phi_prim = angle(cosd((phi(t,n,m)))*sind((theta(t,n,m)))*cos(bita)-sin(bita)*cosd((theta(t,n,m)))+1j*sind((phi(t,n,m)))*sind((theta(t,n,m))));
        phi_prim = rad2deg(phi_prim);
        phi_prim = wrapTo180(phi_prim);
        A_vertical = -min(12*((theta_prim-90)/theta_3dB)^2,SLA_v);
        A_horizontal = -min(12*(phi_prim/phi_3dB)^2,A_max);
        A = -min(-(A_vertical+A_horizontal),A_max);
        F_theta(t,n,m,:) = sqrt(A)*cosd(Polar_Slant_Angle')*cos(psi)-sqrt(A)*sind(Polar_Slant_Angle')*sin(psi);
        F_phi(t,n,m,:) = sqrt(A)*cosd(Polar_Slant_Angle')*sin(psi)+sqrt(A)*sind(Polar_Slant_Angle')*cos(psi);
    end
end
end
