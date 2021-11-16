function [outClust] = greedyNMIalign(sets)
    
    %input: sets: electrodes X cluster schemes
    


    
    %% first make the input sets fit the form 
    %1) label of item 1 = 1 (or -1, and first clustered electrode is
    %labeled as 1)
    %2) label of item i+1 is <= max( of all labels used in items 1:i) + 1
    %3) keep -1 label for "noise" electrodes
    for ii = 1:size(sets,2) %loop over sets
        curID = 1; 
        curS = sets(:,ii);  
        outSet = zeros(length(curS),1); 
        for tt = 1:size(sets,1) %loop over electrodes
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
        sets(:,ii) = outSet; 
    end
        
    %% get the starting best average nmi and use it as a starting seed
    test = mean(reshape(arrayfun(@(x,y) nmi(sets(:,x), sets(:,y)), ...
                reshape(repmat([1:size(sets,2)],[size(sets,2),1]), [size(sets,2)^2,1]), ...
                reshape(repmat([1:size(sets,2)],[size(sets,2),1])', [size(sets,2)^2,1])), ...
                [size(sets,2), size(sets,2)]), 'omitnan');
    [valToBeat, seedI] = max(test); 
    outClust = sets(:,seedI);  
    IDs = unique(sets); %possible cluster IDs
    
    
    check = true; 
    %loop through all trodes and try changing each to each possible ID,
    %recheck aNMI, if better accept change, if not, don't, stop looping
    %when no changes are made for a whole run through the trodes
    while check
        noChange = true; 
        for ii = 1:length(outClust)
            for kk = 1:length(IDs)
                testSet = outClust; 
                testSet(ii) = IDs(kk); 
                test = mean(arrayfun(@(x) nmi(testSet, sets(:,x)), [1:size(sets,2)]), 'omitnan');
                if unique(testSet) == -1
                    'hold on'
                end
                if test > valToBeat
                    noChange = false; 
                    outClust = testSet;
                    valToBeat = test; 
                end
            end
        end
        if noChange
            check = false;
        end
    end

end