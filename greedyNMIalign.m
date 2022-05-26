function [outClust] = greedyNMIalign(sets, maxIter)
    
    %this function takes in a matrix where the rows indicate variables and
    %the columns indicate different clusterings of those variables. All
    %values should be positive integers except for the value -1, which is
    %reserved to indicate unclustered/noise. The input sets are all made to
    %fit a standardized format. Then the average normalized mutual
    %information between each set and all of the other sets is calculated.
    %The set with the best starting average normalized mutual information
    %is used as the seed for the algorithm. A greedy search algorithm loops
    %through the variables changing each one, one at a time, to all of the
    %possible clusters and recalculating normalized mutual information at
    %each step. Changes that increase the normalized mutual information are
    %accepted. This loop continues until no changes are made. The resulting
    %cluster scheme is returned. 

    %input: 
    %   sets: variables X cluster schemes 

    %output:
    %   outClust: a vector with a single cluster scheme that has maximal
    %             average mutual information with all of the input sets
    
    %Adam Dede, adam.osman.dede@gmail.com, Fall 2021
    
    %Ref: Strehl, A., Ghosh, J. (2002). Cluster Ensembles--A knowledge reuse
    %framework for combining multiple partitions. J Machine Learning Research

    
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
                [size(sets,2), size(sets,2)])-eye(size(sets,2)), 'omitnan');
    [~, seedI] = max(test); 
    outClust = sets(:,seedI);  
    valToBeat = mean(arrayfun(@(x) nmi(outClust, sets(:,x)), [1:size(sets,2)]), 'omitnan');
    IDs = unique(sets); %possible cluster IDs
    

    check = true; 
    loopi = 1; 
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
                if test > valToBeat
                    noChange = false; 
                    outClust = testSet;
                    valToBeat = test; 
                end
            end
        end
        loopi = loopi + 1;
        if noChange || loopi > maxIter
            check = false;
        end
    end

end