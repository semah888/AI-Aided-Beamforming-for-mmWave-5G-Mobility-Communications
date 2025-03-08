function [Theta_ZOD_Rays] = ZOD_RaysGenerator(Theta_ZOD,N_Ray,mu_lgZSD)
[N_Cluster,T] = size(Theta_ZOD);
T = T-1;
Alpha = [0.0447;-0.0447;0.1413;-0.1413;0.2492;-0.2492;0.3715;-0.3715;0.5129;-0.5129;0.6797;-0.6797;0.8844;-0.8844;1.1481;-1.1481;1.5195;-1.5195;2.1551;-2.1551];
Theta_ZOD_Rays = zeros(N_Ray,N_Cluster,T+1);
for t = 0:T
for m=1:N_Ray
    for n = 1:N_Cluster
        Theta_ZOD_Rays(m,n,t+1) = Theta_ZOD(n,t+1) + (3/8)*(10^mu_lgZSD)*Alpha(m,1);

    end
end
end
Theta_ZOD_Rays = wrapTo360(Theta_ZOD_Rays);
for t = 0:T
    for m=1:N_Ray
        for n = 1:N_Cluster
        if 180<=Theta_ZOD_Rays(m,n,t+1)
            Theta_ZOD_Rays(m,n,t+1) = 360-Theta_ZOD_Rays(m,n,t+1);
        end
        end
    end
end