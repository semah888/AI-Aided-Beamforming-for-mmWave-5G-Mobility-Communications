function [Theta_ZOD] = ZODGenerator(Powers,ZSD,LOS_Condition,k_factor,Theta_LOS_ZOD,ZOD_offset)
[N,T] = size(Powers);
T = T-1;
Theta_ZOD_prim = zeros(N,T+1);
Theta_ZOD = zeros(N,T+1);
C_theta_NLOS = zeros(25,1);
C_theta_NLOS(8,1) = 0.889;
C_theta_NLOS(10,1) = 0.957;
C_theta_NLOS(11,1) = 1.031;
C_theta_NLOS(12,1) = 1.104;
C_theta_NLOS(15,1) = 1.1088;
C_theta_NLOS(19,1) = 1.184;
C_theta_NLOS(20,1) = 1.178;
C_theta_NLOS(25,1) = 1.282;
C_theta = 0;
for t = 0:T
if LOS_Condition == 0
    C_theta = C_theta_NLOS(N,1);
elseif LOS_Condition == 1
    C_theta = C_theta_NLOS(N,1)*(1.3086+0.0339*k_factor(1,t+1)-0.0077*k_factor(1,t+1)^2+0.0002*k_factor(1,t+1)^3);
end
for n = 1:N
    Theta_ZOD_prim(n,t+1) = -ZSD(1,t+1)*log(Powers(n,t+1)/max(Powers(:,t+1)))/C_theta;
end
Y = zeros(N,1);
X = zeros(N,1);
if LOS_Condition == 0
  for n = 1:N
      X(n,1) = -1+2*rand;
      Y(n,1) = normrnd(0,(ZSD(1,t+1)/7)^2);
    Theta_ZOD(n,t+1) = X(n,1)*Theta_ZOD_prim(n,t+1)+Y(n,1)+Theta_LOS_ZOD(1,t+1)+ZOD_offset;
  end
elseif LOS_Condition == 1 
        X(1,1) = -1+2*rand;
        Y(1,1) = normrnd(0,(ZSD(1,t+1)/7)^2);
      Theta_ZOD(1,t+1) = Theta_LOS_ZOD(1,t+1);
  for n = 2:N
      X(n,1) = -1+2*rand;
      Y(n,1) = normrnd(0,(ZSD(1,t+1)/7)^2);
    Theta_ZOD(n,t+1) = (X(n,1)*Theta_ZOD_prim(n,t+1)+Y(n,1))-(X(1,1)*Theta_ZOD_prim(1,t+1)+Y(1,1)-Theta_LOS_ZOD(1,t+1));
  end
end
end
Theta_ZOD = wrapTo360(Theta_ZOD);
for t = 0:T
    for n = 1:N
        if 180<=Theta_ZOD(n,t+1)
            Theta_ZOD(n,t+1) = 360-Theta_ZOD(n,t+1);
        end
    end
end
end

