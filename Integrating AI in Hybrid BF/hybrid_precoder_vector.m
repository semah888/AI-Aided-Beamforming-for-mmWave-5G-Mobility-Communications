function [F_digital] = hybrid_precoder_vector(FRF, WRF, H, NS)
    % FRF: RF precoder matrix formed with the 8 first beams
    % H: Channel matrix for the user
    % NS: Number of data streams (1 in this case)
    % NRF: Number of RF chains (8 in this case)
    
    % Compute the equivalent channel
    Heq = WRF' * H * FRF;
    
    % Perform SVD on the equivalent channel
    [U, S, V] = svd(Heq);
    
    % Digital precoder is the first right singular vector of the equivalent channel
    F_digital = V(:, 1);
    m = size(F_digital);
    % Normalize the digital precoder
    %FBB = sqrt(NS) * (F_digital / norm(FRF * F_digital, 'fro'));
    
end
