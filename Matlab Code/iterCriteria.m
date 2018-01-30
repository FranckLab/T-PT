function [iterCrit,sSize] = iterCriteria(sSize,nTracked,iterCheck)

% Function to check if further iterations need to be run and provice subset
% size for next iteration

minR = 0.6;
maxSS = 5;

if length(sSize) > 1 % Check at iter>1
    nTrackedDiff = diff(nTracked);
    lastDiff = nTrackedDiff(end);
    ratioT = nTrackedDiff(end)/nTrackedDiff(end-1);
    currentSS = sum(sSize == sSize(end));
    uniqueSS = length(unique(sSize));
    
    % Maximum iteration criteria
    if iterCheck == 1   %max iter criteria
        iterCrit = 0;
    end
    
    
    if uniqueSS == 5
        if lastDiff<0.0025 | ratioT<0.1 | currentSS>=maxSS
            iterCrit = 0;
            sSize(end+1) = sSize(end);
        else
            iterCrit = 1;
            sSize(end+1) = sSize(end);
        end
       
    elseif lastDiff<0.04
            sSize(end+1) = sSize(end)/2;
            iterCrit = 1;
    elseif ratioT<minR
            sSize(end+1) = sSize(end)/2;
            iterCrit = 1;
        
    else        
        if currentSS<maxSS
            sSize(end+1) = sSize(end);
            iterCrit = 1;
        else
            sSize(end+1) = sSize(end)/2;
            iterCrit = 1;
        end
    end
    
else %After first iteration
    iterCrit = 1;
    sSize(2) = sSize(1);
end

if sSize(end) ==8
    sSize(end) = 16;
end


end

