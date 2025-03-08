function [Powers] = PowerGenerator(Delays,DS,r_taux,Per_Cluster_Shadowing,LOS_Condition,k_factor)
[N_Cluster,T] = size(Delays);
T = T-1;
Powers_prim = zeros(N_Cluster,T+1);
Powers = zeros(N_Cluster,T+1);
z = zeros(N_Cluster,T+1);
for t = 0:T
    z(:,t+1) = normrnd(0,Per_Cluster_Shadowing^2,[N_Cluster,1]);
    for n = 1:N_Cluster
        Powers_prim(n,t+1) = exp(-Delays(n,t+1)*(r_taux-1)/(r_taux*DS(1,t+1)))*10^(-z(n,t+1)/10);
    end
if LOS_Condition == 1
    Powers(1,t+1) = (1/(10^(k_factor(1,t+1)/10)+1))*(Powers_prim(1,t+1)/sum(Powers_prim(:,t+1)))+(10^(k_factor(1,t+1)/10)/(10^(k_factor(1,t+1)/10)+1));
    for n=2:N_Cluster
        Powers(n,t+1) = (1/(10^(k_factor(1,t+1)/10)+1))*(Powers_prim(n,t+1)/sum(Powers_prim(:,t+1)));
    end
elseif LOS_Condition == 0
    for n =1:N_Cluster
        Powers(n,t+1) = Powers_prim(n,t+1)/sum(Powers_prim(:,t+1));
    end
end
end
end


