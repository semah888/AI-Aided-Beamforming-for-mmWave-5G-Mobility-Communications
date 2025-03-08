function [Channel_cluster,Channel] = ChannelCoefficientGenerator(t,k,initial,N_slots_model,numSC,Ntx_v,Ntx_h,Nrx_v,Nrx_h,d_rx,d_tx,Polar_Slant_Angle,NumCluster,NumRay,Delays,Powers,Theta_ZOA_Ray,Phi_AOA_Ray,Theta_ZOD_Ray,Phi_AOD_Ray,lambda,XPR,Velocity,C_DS,LOS_Condition,k_factor,Phi_LOS_AOA,Phi_LOS_AOD,Theta_LOS_ZOA,Theta_LOS_ZOD,d_3D,UpdatePeriod)
Ntx = Ntx_h*Ntx_v;
Nrx = Nrx_h*Nrx_v;
[p,numPolarAngle] = size(Polar_Slant_Angle);
Channel = zeros(numPolarAngle*Nrx,numPolarAngle*Ntx);
Channel_NLOS_cluster_Ray = zeros(numPolarAngle*Nrx,numPolarAngle*Ntx,NumCluster,NumRay);
[F_rx_theta,F_rx_phi] = RadiationPatternGenerator(Theta_ZOA_Ray,Phi_AOA_Ray,Polar_Slant_Angle);
[F_tx_theta,F_tx_phi] = RadiationPatternGenerator(Theta_ZOD_Ray,Phi_AOD_Ray,Polar_Slant_Angle);
[theta_theta,theta_phi,phi_theta,phi_phi] = Draw_initial_phases(NumCluster,NumRay);
X_Rx = zeros(1,Nrx);
[Y_Rx,Z_Rx] = meshgrid(0:Nrx_h-1,0:Nrx_v-1);
RxPositions = [X_Rx;Y_Rx(:).';Z_Rx(:).']*d_rx; % 3*Nrx
X_Tx = zeros(1,Ntx);
[Y_Tx,Z_Tx] = meshgrid(0:Ntx_h-1,0:Ntx_v-1);
TxPositions = [X_Tx;Y_Tx(:).';Z_Tx(:).']*d_tx; % 3*Ntx
 for u = 1:numPolarAngle*Nrx
    for s = 1:numPolarAngle*Ntx
        for n = 1:NumCluster
            for m = 1:NumRay
                r_chap_rx = [sind(Theta_ZOA_Ray(m,n,t-initial+1))*cosd(Phi_AOA_Ray(m,n,t-initial+1));sind(Theta_ZOA_Ray(m,n,t-initial+1))*sind(Phi_AOA_Ray(m,n,t-initial+1));cosd(Theta_ZOA_Ray(m,n,t-initial+1))];
                r_chap_tx = [sind(Theta_ZOD_Ray(m,n,t-initial+1))*cosd(Phi_AOD_Ray(m,n,t-initial+1));sind(Theta_ZOD_Ray(m,n,t-initial+1))*sind(Phi_AOD_Ray(m,n,t-initial+1));cosd(Theta_ZOD_Ray(m,n,t-initial+1))];
                if mod(u,numPolarAngle)==0
                    f_rx_theta = F_rx_theta(m,n,t-initial+1,1);
                    f_rx_phi = F_rx_phi(m,n,t-initial+1,1);
                else
                    f_rx_theta = F_rx_theta(m,n,t-initial+1,2);
                    f_rx_phi = F_rx_phi(m,n,t-initial+1,2);
                end
                 if mod(s,numPolarAngle)==0
                    f_tx_theta = F_tx_theta(m,n,t-initial+1,1);
                    f_tx_phi = F_tx_phi(m,n,t-initial+1,1);
                else
                    f_tx_theta = F_tx_theta(m,n,t-initial+1,2);
                    f_tx_phi = F_tx_phi(m,n,t-initial+1,2);
                 end
                 if numPolarAngle == 2
                     d_rx_u = RxPositions(:,fix((u+1)/numPolarAngle));
                     d_tx_s = TxPositions(:,fix((s+1)/numPolarAngle));
                 elseif numPolarAngle == 1
                d_rx_u = RxPositions(:,u);
                d_tx_s = TxPositions(:,s);
                 end
                Channel_NLOS_cluster_Ray(u,s,n,m) = sqrt(Powers(n,t-initial+1)/NumRay)*[f_rx_theta f_rx_phi]*[exp(1i*theta_theta(n,m)) exp(1i*phi_theta(n,m))/sqrt(XPR(t-initial+1,n,m));exp(1i*theta_phi(n,m))/sqrt(XPR(t-initial+1,n,m)) exp(1i*phi_phi(n,m))]*[f_tx_theta;f_tx_phi]*exp(2i*pi/lambda*r_chap_rx'*d_rx_u)*exp(2i*pi/lambda*r_chap_tx'*d_tx_s)*exp(2i*pi*((r_chap_rx'*Velocity)/lambda)*t*UpdatePeriod);
            end
        end
    end
 end
 Channel_cluster = sum(Channel_NLOS_cluster_Ray,4);
Channel_NLOS = zeros(numPolarAngle*Nrx,numPolarAngle*Ntx);
%if C_DS ==0
    for u = 1:numPolarAngle*Nrx
        for s = 1:numPolarAngle*Ntx
            for n=1:NumCluster
                for m=1:NumRay
                    Channel_NLOS(u,s) = Channel_NLOS(u,s)+Channel_NLOS_cluster_Ray(u,s,n,m)*exp(-1j*2*pi*(n-1)*k/numSC);
                end
            end
        end
    end
%else
 %   for u = 1:numPolarAngle*Nrx
 %       for s = 1:numPolarAngle*Ntx
 %           for n=3:NumCluster
 %               for m=1:NumRay
%                    Channel_NLOS(u,s) = Channel_NLOS(u,s)+Channel_NLOS_cluster_Ray(u,s,n,m)*exp(-1j*2*pi*n*k/numSC);
%                end
 %           end
%        end
 %   end
%    Delays_1_2 = zeros(2,3);
%    for n = 1:2
 %     Delays_1_2(n,1) = Delays(n,t+1);
 %     Delays_1_2(n,2) = Delays(n,t+1)+1.28*C_DS;
 %     Delays_1_2(n,3) = Delays(n,t+1)+2.56*C_DS;
 %   end
%    for u = 1:Nrx
 %   for s = 1:Ntx
 %      for n=1:2
 %          for m=1:8
 %              Channel_NLOS(u,s) = Channel_NLOS(u,s)+Channel_NLOS_cluster_Ray(u,s,n,m)*delta(taux-Delays_1_2(n,1));
 %          end
 %          for m=19:20
 %            Channel_NLOS(u,s) = Channel_NLOS(u,s)+Channel_NLOS_cluster_Ray(u,s,n,m)*delta(taux-Delays_1_2(n,1));
 %          end
 %         for m=9:12
 %              Channel_NLOS(u,s) = Channel_NLOS(u,s)+Channel_NLOS_cluster_Ray(u,s,n,m)*delta(taux-Delays_1_2(n,2));
 %          end
  %         for m=17:18
  %           Channel_NLOS(u,s) = Channel_NLOS(u,s)+Channel_NLOS_cluster_Ray(u,s,n,m)*delta(taux-Delays_1_2(n,2));
  %         end
  %         for m=13:16
  %             Channel_NLOS(u,s) = Channel_NLOS(u,s)+Channel_NLOS_cluster_Ray(u,s,n,m)*delta(taux-Delays_1_2(n,3));
  %         end
%       end
%   end
%    end
%end
    if LOS_Condition == 0
   Channel = Channel_NLOS;
elseif LOS_Condition == 1
    Channel_LOS = zeros(Nrx,Ntx);
    Theta_LOS_ZOA_matrix = zeros(1,1,1);
    Theta_LOS_ZOA_matrix(1,1,1) = Theta_LOS_ZOA(1,t);
    Phi_LOS_AOA_matrix = zeros(1,1,1);
    Phi_LOS_AOA_matrix(1,1,1) = Phi_LOS_AOA(1,t);
    Theta_LOS_ZOD_matrix = zeros(1,1,1);
    Theta_LOS_ZOD_matrix(1,1,1) = Theta_LOS_ZOD(1,t);
    Phi_LOS_AOD_matrix = zeros(1,1,1);
    Phi_LOS_AOD_matrix(1,1,1) = Phi_LOS_AOA(1,t);
    [F_rx_LOS_theta,F_rx_LOS_phi] = RadiationPatternGenerator(Theta_LOS_ZOA_matrix,Phi_LOS_AOA_matrix,Polar_Slant_Angle);
    [F_tx_LOS_theta,F_tx_LOS_phi] = RadiationPatternGenerator(Theta_LOS_ZOD_matrix,Phi_LOS_AOD_matrix,Polar_Slant_Angle);
     r_chap_rx = [sind(Theta_LOS_ZOA(1,t))*cosd(Phi_LOS_AOA(1,t));sind(Theta_LOS_ZOA(1,t))*sind(Phi_LOS_AOA(1,t));cosd(Theta_LOS_ZOA(1,t))];
     r_chap_tx = [sind(Theta_LOS_ZOD(1,t))*cosd(Phi_LOS_AOD(1,t));sind(Theta_LOS_ZOD(1,t))*sind(Phi_LOS_AOD(1,t));cosd(Theta_LOS_ZOD(1,t))];
     for u = 1:numPolarAngle*Nrx
       for s = 1:numPolarAngle*Ntx
           if mod(u,numPolarAngle)==0
                    f_rx_theta = F_rx_LOS_theta(1,1,1,1);
                    f_rx_phi = F_rx_LOS_phi(1,1,1,1);
           else
                    f_rx_theta = F_rx_LOS_theta(1,1,1,2);
                    f_rx_phi = F_rx_LOS_phi(1,1,1,2);
            end
                 if mod(s,numPolarAngle)==0
                    f_tx_theta = F_tx_LOS_theta(1,1,1,1);
                    f_tx_phi = F_tx_LOS_phi(1,1,1,1);
                else
                    f_tx_theta = F_tx_LOS_theta(1,1,1,2);
                    f_tx_phi = F_tx_LOS_phi(1,1,1,2);
                 end
                 if numPolarAngle == 2
                     d_rx_u = RxPositions(:,fix((u+1)/numPolarAngle));
                     d_tx_s = TxPositions(:,fix((s+1)/numPolarAngle));
                 elseif numPolarAngle == 1
                d_rx_u = RxPositions(:,u);
                d_tx_s = TxPositions(:,s);
                 end
          Channel_LOS(u,s) = [f_rx_theta f_rx_phi]*[1 0;0 -1]*[f_tx_theta;f_tx_phi]*exp(-2i*pi*d_3D(1,t)/lambda)*exp(2i*pi/lambda*r_chap_rx'*d_rx_u)*exp(2i*pi/lambda*r_chap_tx'*d_tx_s)*exp(2i*pi*(r_chap_rx'*Velocity/lambda)*t);
        end
     end
Channel = sqrt(1/(10^(k_factor(1,t)/10)+1))*Channel_NLOS+sqrt(10^(k_factor(1,t)/10)/(10^(k_factor(1,t)/10)+1))*(1/sqrt(numSC))*Channel_LOS;
    end
end
