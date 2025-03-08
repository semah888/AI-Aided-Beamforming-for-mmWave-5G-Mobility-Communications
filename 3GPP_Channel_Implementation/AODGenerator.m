function [Phi_AOD] = AODGenerator(Powers,ASD,LOS_Condition,k_factor,Phi_LOS_AOD)
[N,T] = size(Powers);
T = T-1;
Phi_AOD_prim = zeros(N,T+1);
Phi_AOD = zeros(N,T+1);
C_phi_NLOS = zeros(25,1);
C_phi_NLOS(4,1) = 0.779;
C_phi_NLOS(5,1) = 0.860;
C_phi_NLOS(8,1) = 1.018;
C_phi_NLOS(10,1) = 1.090;
C_phi_NLOS(11,1) = 1.123;
C_phi_NLOS(12,1) = 1.146;
C_phi_NLOS(14,1) = 1.190;
C_phi_NLOS(15,1) = 1.211;
C_phi_NLOS(16,1) = 1.226;
C_phi_NLOS(19,1) = 1.273;
C_phi_NLOS(20,1) = 1.289;
C_phi_NLOS(25,1) = 1.358;
C_phi = 0;
for t=0:T
if LOS_Condition == 0
    C_phi = C_phi_NLOS(N,1);
elseif LOS_Condition == 1
    C_phi = C_phi_NLOS(N,1)*(1.1035-0.028*k_factor(1,t+1)-0.002*k_factor(1,t+1)^2+0.0001*k_factor(1,t+1)^3);
end
for n = 1:N
    Phi_AOD_prim(n,t+1) = 2*(ASD(1,t+1)/1.4)*sqrt(-log(Powers(n,t+1)/max(Powers(:,t+1))))/C_phi;
end
Y = zeros(N,1);
X = zeros(N,1);
if LOS_Condition == 0
  for n = 1:N
      X(n,1) = -1+2*rand;
      Y(n,1) = normrnd(0,(ASD(1,t+1)/7)^2);
    Phi_AOD(n,t+1) = X(n,1)*Phi_AOD_prim(n,t+1)+Y(n,1)+Phi_LOS_AOD(1,t+1);
  end
elseif LOS_Condition == 1 
        X(1,1) = -1+2*rand;
        Y(1,1) = normrnd(0,(ASD(1,t+1)/7)^2);
      Phi_AOD(1,t+1) = Phi_LOS_AOD(1,t+1);
  for n = 2:N
      X(n,1) = -1+2*rand;
      Y(n,1) = normrnd(0,(ASD(1,t+1)/7)^2);
    Phi_AOD(n,t+1) = (X(n,1)*Phi_AOD_prim(n,t+1)+Y(n,1))-(X(1,1)*Phi_AOD_prim(1,t+1)+Y(1,1)-Phi_LOS_AOD(1,t+1));
  end
end
end
Phi_AOD = wrapTo360(Phi_AOD);
end
