function [XPR] = XPRGenerator(T,N_Ray,N_Cluster,mu_XPR,sigma_XPR)
XPR = zeros(T,N_Cluster,N_Ray);
for t = 1:T
X = zeros(N_Cluster,N_Ray);
for n = 1:N_Cluster
    for m=1:N_Ray
        X(n,m) = normrnd(mu_XPR,sigma_XPR^2);
        XPR(t,n,m) = 10^(X(n,m)/10);
    end
end
end
end
