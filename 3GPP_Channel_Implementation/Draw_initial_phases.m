function [Theta_theta,Theta_phi,phi_theta,phi_phi] = Draw_initial_phases(N_Cluster,N_Ray)
max = pi;
min = -pi;
Theta_theta = rand(N_Cluster,N_Ray)*(max-min)+min;
Theta_phi = rand(N_Cluster,N_Ray)*(max-min)+min; 
phi_theta = rand(N_Cluster,N_Ray)*(max-min)+min; 
phi_phi = rand(N_Cluster,N_Ray)*(max-min)+min; % [-pi,pi]
end