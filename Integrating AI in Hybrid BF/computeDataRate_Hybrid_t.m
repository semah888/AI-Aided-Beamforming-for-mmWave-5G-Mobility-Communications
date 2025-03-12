function R = computeDataRate_Hybrid_t(t,sc,FRF,WRF,freqChan,SNR,FBB,WBB,Ns)
    % H: Channel matrix for the user
    % FRF: RF precoder matrix
    % F_digital: Digital precoder vector
    % W_digital: Digital combiner vector
    % noise_power: Noise power (scalar)
    H = freqChan(:,:,sc,t);
    % Compute the equivalent channel
    % Compute overall precoder and combiner
F = FRF * FBB;
W = WRF * WBB;

% Compute the effective channel
H_eff = W' * H * F;

% Compute the achievable spectral efficiency
R = log2(det(eye(Ns) + (SNR / Ns) *abs(H_eff * H_eff')));
end
