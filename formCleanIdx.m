function [outSet] = formCleanIdx(curS)

        curID = 1;  
        outSet = zeros(length(curS),1); 
        for tt = 1:length(curS) %loop over electrodes
            if curS(tt)==-1 %don't change noise
                outSet(tt) = -1; 
            else
                idx = find(curS(1:tt-1)==curS(tt));
                if ~isempty(idx) %this cluster has been seen before, so grab whatever label was used before
                    outSet(tt) = outSet(idx(1)); 
                else %this is the first instance of this cluster
                    outSet(tt) = curID; 
                    curID = curID + 1; 
                end
            end
        end

end