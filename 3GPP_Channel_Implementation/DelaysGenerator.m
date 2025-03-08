function [Delays] = DelaysGenerator(T,r_taux,DS,N_Cluster,k_factor)
Delays = zeros(N_Cluster,T);
for t = 1:T
for n=1:N_Cluster
    Delays(n,t) = -r_taux*DS(1,t)*rand;
end
Delays(:,t) = sort(Delays(:,t)-min(Delays(:,t)));%/(0.7705-0.0433*k_factor(1,t)+0.0002*k_factor(1,t)^2+0.000017*k_factor(1,t)^3);
end

