function [Phi_AOA_Rays] = AOA_RaysGenerator(Phi_AOA,C_ASA,N_Ray)
[N_Cluster,T] = size(Phi_AOA);
T = T-1;
Alpha = [0.0447;-0.0447;0.1413;-0.1413;0.2492;-0.2492;0.3715;-0.3715;0.5129;-0.5129;0.6797;-0.6797;0.8844;-0.8844;1.1481;-1.1481;1.5195;-1.5195;2.1551;-2.1551];
Phi_AOA_Rays = zeros(N_Ray,N_Cluster,T+1);
for t = 0:T
    for m=1:N_Ray
        for n = 1:N_Cluster
            Phi_AOA_Rays(m,n,t+1) = Phi_AOA(n,t+1) + C_ASA*Alpha(m,1);
        end
    end
end
Phi_AOA_Rays = wrapTo360(Phi_AOA_Rays);