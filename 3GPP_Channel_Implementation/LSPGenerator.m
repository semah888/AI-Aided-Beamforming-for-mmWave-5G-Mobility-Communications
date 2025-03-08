function [DS,ASD,ASA,ZSD,ZSA,K,SF] = LSPGenerator(Updatedist,T,mulgDS,sigmalgDS,mulgASD,sigmalgASD,mulgASA,sigmalgASA,mulgZSD,sigmalgZSD,mulgZSA,sigmalgZSA,muSF,sigmaSF,muK,sigmaK,C_ASD_DS,C_ASA_DS,C_ASA_SF,C_ASD_SF,C_DS_SF,C_ASD_ASA,C_ASD_K,C_ASA_K,C_DS_K,C_SF_K,C_ZSD_DS,C_ZSA_DS,C_ZSA_SF,C_ZSD_SF,C_ZSD_ZSA,C_ZSD_K,C_ZSA_K,C_ZSA_ASA,C_ZSD_ASA,C_ZSA_ASD,C_ZSD_ASD,dcorr_DS,dcorr_K,dcorr_SF,dcorr_ASA,dcorr_ASD,dcorr_ZSA,dcorr_ZSD)
Corr_matrix = zeros(7*(T),7*(T));
%Level 1
C_LSPt = zeros(7,7);
C_LSPt(1,1) = sigmalgDS^2;
C_LSPt(2,2) = sigmalgASD^2;
C_LSPt(3,3) = sigmalgASA^2;
C_LSPt(4,4) = sigmalgZSD^2;
C_LSPt(5,5) = sigmalgZSA^2;
C_LSPt(6,6) = sigmaK^2;
C_LSPt(7,7) = sigmaSF^2;
C_LSPt(1,2) = sigmalgDS*sigmalgASD*C_ASD_DS;
C_LSPt(1,3) = sigmalgDS*sigmalgASA*C_ASA_DS;
C_LSPt(1,4) = sigmalgDS*sigmalgZSD*C_ZSD_DS;
C_LSPt(1,5) = sigmalgDS*sigmalgZSA*C_ZSA_DS;
C_LSPt(1,6) = sigmalgDS*sigmaK*C_DS_K;
C_LSPt(1,7) = sigmalgDS*sigmaSF*C_DS_SF;
C_LSPt(2,1) = sigmalgDS*sigmalgASD*C_ASD_DS;
C_LSPt(2,3) = sigmalgASD*sigmalgASA*C_ASD_ASA;
C_LSPt(2,4) = sigmalgASD*sigmalgZSD*C_ZSD_ASD;
C_LSPt(2,5) = sigmalgASD*sigmalgZSA*C_ZSA_ASD;
C_LSPt(2,6) = sigmalgASD*sigmaK*C_ASD_K;
C_LSPt(2,7) = sigmalgASD*sigmaSF*C_ASD_SF;
C_LSPt(3,1) = sigmalgDS*sigmalgASA*C_ASA_DS;
C_LSPt(3,2) = sigmalgASD*sigmalgASA*C_ASD_ASA;
C_LSPt(3,4) = sigmalgASA*sigmalgZSD*C_ZSD_ASA;
C_LSPt(3,5) = sigmalgASA*sigmalgZSA*C_ZSA_ASA;
C_LSPt(3,6) = sigmalgASA*sigmaK*C_ASA_K;
C_LSPt(3,7) = sigmalgASA*sigmaSF*C_ASA_SF;
C_LSPt(4,1) = sigmalgDS*sigmalgZSD*C_ZSD_DS;
C_LSPt(4,2) = sigmalgASD*sigmalgZSD*C_ZSD_ASD;
C_LSPt(4,3) = sigmalgASA*sigmalgZSD*C_ZSD_ASA;
C_LSPt(4,5) = sigmalgZSD*sigmalgZSA*C_ZSD_ZSA;
C_LSPt(4,6) = sigmalgZSD*sigmaK*C_ZSD_K;
C_LSPt(4,7) = sigmalgZSD*sigmaSF*C_ZSD_SF;
C_LSPt(5,1) = sigmalgDS*sigmalgZSA*C_ZSA_DS;
C_LSPt(5,2) = sigmalgASD*sigmalgZSA*C_ZSA_ASD;
C_LSPt(5,3) = sigmalgASA*sigmalgZSA*C_ZSA_ASA;
C_LSPt(5,4) = sigmalgZSD*sigmalgZSA*C_ZSD_ZSA;
C_LSPt(5,6) = sigmalgZSA*sigmaK*C_ZSA_K;
C_LSPt(5,7) = sigmalgZSA*sigmaSF*C_ZSA_SF;
C_LSPt(6,1) = sigmalgDS*sigmaK*C_DS_K;
C_LSPt(6,2) = sigmalgASD*sigmaK*C_ASD_K;
C_LSPt(6,3) = sigmalgASA*sigmaK*C_ASA_K;
C_LSPt(6,4) = sigmalgZSD*sigmaK*C_ZSD_K;
C_LSPt(6,5) = sigmalgZSA*sigmaK*C_ZSA_K;
C_LSPt(6,7) = sigmaK*sigmaSF*C_SF_K;
C_LSPt(7,1) = sigmalgDS*sigmaSF*C_DS_SF;
C_LSPt(7,2) = sigmalgASD*sigmaSF*C_ASD_SF;
C_LSPt(7,3) = sigmalgASA*sigmaSF*C_ASA_SF;
C_LSPt(7,4) = sigmalgZSD*sigmaSF*C_ZSD_SF;
C_LSPt(7,5) = sigmalgZSA*sigmaSF*C_ZSA_SF;
C_LSPt(7,6) = sigmaK*sigmaSF*C_SF_K;
%Diag of covariance matrix
for t = 0:T-1
    Corr_matrix(7*t+1:7*t+7,7*t+1:7*t+7) = C_LSPt;
end
%Level 2
%rest of cov matrix
for i =1:7*(T)
    for j =1:7*(T)
        if mod(i,7) == mod(j,7)
        k1 = mod(i,7);
        k2 = mod(j,7);
        if k1==1 && k2 ==1
            ti = floor(i/7);
            tj = floor(j/7);
            xi = Updatedist*ti;
            xj = Updatedist*tj;
            Corr_matrix(i,j) = exp(sqrt((xi-xj)^2)/dcorr_DS)*sigmalgDS^2;
           
     elseif k1==2 && k2 ==2
            ti = floor(i/7);
            tj = floor(j/7);
            xi = Updatedist*ti; %cas particulier directional horizontal
            xj = Updatedist*tj;
            Corr_matrix(i,j) = exp(sqrt((xi-xj)^2)/dcorr_ASD)*sigmalgASD^2;
     elseif k1==3 && k2 ==3
            ti = floor(i/7);
            tj = floor(j/7);
            xi = Updatedist*ti; %cas particulier directional horizontal
            xj = Updatedist*tj;
            Corr_matrix(i,j) = exp(sqrt((xi-xj)^2)/dcorr_ASA)*sigmalgASA^2;
     elseif k1==4 && k2 ==4
            ti = (i-k1)/7;
            tj = (j-k1)/7;
            xi = Updatedist*ti; %cas particulier directional horizontal
            xj = Updatedist*tj;
            Corr_matrix(i,j) = exp(sqrt((xi-xj)^2)/dcorr_ZSD)*sigmalgZSD^2;
      elseif k1==5 && k2 ==5
            ti = floor(i/7);
            tj = floor(j/7);
            xi = Updatedist*ti; %cas particulier directional horizontal
            xj = Updatedist*tj;
            Corr_matrix(j,i) = exp(sqrt((xi-xj)^2)/dcorr_ZSA)*sigmalgZSA^2;

       elseif k1==6 && k2 ==6
            ti = floor(i/7);
            tj = floor(j/7);
            xi = Updatedist*ti; %cas particulier directional horizontal
            xj = Updatedist*tj;
            Corr_matrix(i,j) = exp(sqrt((xi-xj)^2)/dcorr_K)*sigmaK^2;
    else
            ti = floor(i/7);
            tj = floor(j/7);
            xi = Updatedist*ti; %cas particulier directional horizontal
            xj = Updatedist*tj;
            Corr_matrix(i,j) = exp(sqrt((xi-xj)^2)/dcorr_SF)*sigmaSF^2;
           
        end
    end
    end
end
%Diag bloc of covariance matrix
for t = 0:T-1
    Corr_matrix(7*t+1:7*t+7,7*t+1:7*t+7) = C_LSPt;
end
E_t = [mulgDS mulgASD mulgASA mulgZSD mulgZSA muK muSF];
E = repmat(E_t,1,T);
%sigma = topdm(Corr_matrix);
if eig(Corr_matrix)>=0
    lsp = mvnrnd(E,Corr_matrix,1);
else
sigma = nearestSPD(Corr_matrix);
%R = chol(sigma);
%k = size(R)
lsp = mvnrnd(E,sigma,1); 
end
DS = zeros(1,T);
ASD = zeros(1,T);
ASA = zeros(1,T);
ZSD = zeros(1,T);
ZSA = zeros(1,T);
K = zeros(1,T);
SF = zeros(1,T);
for t = 0:T-1
    DS(1,t+1)  = 10^lsp(1,7*t+1);
    ASD(1,t+1) = min(10^lsp(1,7*t+2),104);
    ASA(1,t+1) = min(10^lsp(1,7*t+3),104);
    ZSD(1,t+1) = min(10^lsp(1,7*t+4),52);
    ZSA(1,t+1) = min(10^lsp(1,7*t+5),52);
    K(1,t+1) =  lsp(1,7*t+6);
    SF(1,t+1) = lsp(1,7*t+7);
end
end
