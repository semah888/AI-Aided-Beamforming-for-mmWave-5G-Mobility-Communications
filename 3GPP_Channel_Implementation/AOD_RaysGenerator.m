function [Phi_AOD_Rays] = AOD_RaysGenerator(Phi_AOD,C_ASD,N_Ray)
[N_Cluster,T] = size(Phi_AOD);
T = T-1;
Alpha = [0.0447;-0.0447;0.1413;-0.1413;0.2492;-0.2492;0.3715;-0.3715;0.5129;-0.5129;0.6797;-0.6797;0.8844;-0.8844;1.1481;-1.1481;1.5195;-1.5195;2.1551;-2.1551];
Phi_AOD_Rays = zeros(N_Ray,N_Cluster,T+1);
for t =0:T
for m=1:N_Ray
    for n = 1:N_Cluster
        Phi_AOD_Rays(m,n,t+1) = Phi_AOD(n,t+1) + C_ASD*Alpha(m,1);
    end
end
end
Phi_AOD_Rays = wrapTo360(Phi_AOD_Rays);
end