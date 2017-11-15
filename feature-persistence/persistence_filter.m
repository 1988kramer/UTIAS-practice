function [probs, misses] = persistence_filter(pM, pF, lambda, marks_in_FOV, ...
                                             probs, z, t, pV, misses) 
    
    % loop over all landmark probabilities
    % update probabilities based on measurements and time elapsed
    for i = 1:size(probs,1)
        if marks_in_FOV(i) == 1 || ismember(i, z(3,:))
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
            
            if yT == 0
                misses(i) = misses(i) + 1;
            end

            % if landmark has been encountered before
            if probs(i,6) > 0
                % compute next prior for given landmark
                % NOTE need to switch this to ln(prior)
                nextPrior = -exp(-lambda * (t - probs(i,6))) + 1;
                % compute ln(partial) evidence
                probs(i,2) = log((pF^yT) * (1 - pF)^(1-yT)) ...
                             + probs(i,2) ...
                             + log(1 + exp(probs(i,3) + log(nextPrior - probs(i,5)) - probs(i,2)));
                % compute ln(likelihood)
                probs(i,3) = log((pM^(1 - yT)) * ((1 - pM)^yT)) + probs(i,3);
                % compute ln(evidence)
                probs(i,4) = probs(i,2) + log(1 + exp(probs(i,3) - (lambda * t) - probs(i,2)));
                % update prior
                probs(i,5) = nextPrior;
                % compute ln(posterior)
                probs(i,1) = probs(i,3) - probs(i,4) - (lambda * t);
                % convert ln(posterior) to posterior
                probs(i,1) = exp(probs(i,1));
                
                % if persistence probability is less than threshold
                % remove from map
                % -1 in probs(i,6) prevents landmark from 
                % being updated again
                if probs(i,1) < pV
                    probs(i,:) = [0 0 1 1 0 -1];
                end
            end
        end
    end
end