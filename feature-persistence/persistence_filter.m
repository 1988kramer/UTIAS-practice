function [probs, prior] = persistence_filter(pM, pF, lambda, prior, probs, z, t)
    % calculate prior probability of landmark
    % persistence at time t
    nextPrior = exp(-lambda * t) - 1; 
    
    % loop over all landmark probabilities
    % update probabilities based on measurements and time elapsed
    % assumes correspondence IS KNOWN
    for i = 1:size(probs,1)
        yT = 0;
        % if there is a measurement of landmark i, set yT to 1
        if ismember(i, z(3,i))
            yT = 1;
        end
        % compute partial evidence
        probs(i,2) = (pF^yT) * (1 - pF)^(1-yT) ...
                     * (probs(i,2) + probs(i,3) ...
                     * (nextPrior - prior));
        % compute likelihood
        probs(i,3) = (pM^(1 - yT)) * (1 - pM)^yT * probs(i,3);
        % compute evidence
        probs(i,4) = probs(i,2) * probs(i,3) * (1 - nextPrior);
        % compute probabilities
        probs(i,1) = (probs(i,3) / probs(i,4)) * (1 - nextPrior);
    end
    
    prior = nextPrior;
end