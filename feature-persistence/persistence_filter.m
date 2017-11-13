function [probs, prior] = persistence_filter(pM, pF, lambda, prior, ...
                                             probs, z, t, pV)
    % calculate prior probability of landmark
    % persistence at time t
    nextPrior = -exp(-lambda * t) + 1; 
    
    % loop over all landmark probabilities
    % update probabilities based on measurements and time elapsed
    % assumes correspondence IS KNOWN
    for i = 1:size(probs,1)
        yT = 0;
        % if there is a measurement of landmark i, set yT to 1
        if ismember(i, z(3,:))
            yT = 1;
            % if current landmark has not been encountered
            if probs(i,6) == 0
                probs(i,5) = 0;
                probs(i,6) = t;
            end
        end
        % if landmark has been encountered before
        if probs(i,6) > 0
            % compute next prior for given landmark
            nextPrior = -exp(-lambda * (t - probs(i,6))) + 1;
            % compute partial evidence
            probs(i,2) = ((pF^yT) * (1 - pF)^(1-yT)) ...
                         * (probs(i,2) ...
                         + (probs(i,3) * (nextPrior - probs(i,5))));
            % compute likelihood
            probs(i,3) = (pM^(1 - yT)) * ((1 - pM)^yT) * probs(i,3);
            % compute evidence
            probs(i,4) = probs(i,2) + probs(i,3) * (1 - nextPrior);
            % update prior
            probs(i,5) = nextPrior;
            % compute posterior probability
            probs(i,1) = (probs(i,3) / probs(i,4)) * (1 - nextPrior);
        end
        
        % if persistence probability is less than threshold
        % remove from map
        if probs(i,1) < pV
            probs(i,:) = [0 0 1 1 0 0];
        end
    end
    
    prior = nextPrior;
end