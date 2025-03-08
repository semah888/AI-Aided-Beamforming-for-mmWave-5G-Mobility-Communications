load('MIMOchannel_details.mat','PropagChanModelVariable');
T = 100
LSPs = PropagChanModelVariable{1,1}.SSP_Clusters
for t = 0:T
    Delays = PropagChanModelVariable{1,t+1}.SSP_Clusters.AbsoluteDelay

end