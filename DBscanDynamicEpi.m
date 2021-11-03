function [idxVals] = DBscanDynamicEpi(corMat, k, d)
    %This function is designed to take as input a correlation matrix for
    %clustering and a value of k to be used as input to the DBscan algorithm.
    %corMat: an n X n symmetric matrix of correlation values
    %k: an integer value well below n 
    
    %output: 
    %idxVals: 1 X n vector of cluster assignments. The -1 cluster indicates
    %those input rows from the corMat that did not get assigned to a cluster
    
    %Adam Dede, adam.osman.dede@gmail.com, Fall 2021
    
    %Refs: 
    %Liu, Zhu, Chen, 2012. J neurosci Methods.
    
    %% make matrix to be clustered
    %convert into a correlation matrix of relationships between connectivity
    %maps rather than between individual channels (Liu et al., 2012)
    distMat = corrcoef(corMat); 
    %convert into a squared distance matrix. This forces positivity, which is
    %necessary for DBscan input and increases contrast
    for ii = 1:d
        distMat = pdist2(distMat, distMat) .^2; 
    end
    if d==0
        distMat = distMat .^2; 
    end
    
    %% find the k-distance for each point
    kdist = zeros(length(distMat),1);
    for ii = 1:length(kdist)
        curRow = distMat(ii,:); 
        curRow = sort(curRow); 
        kdist(ii) = curRow(k); 
    end
    kdist = sort(kdist); 
    if min(kdist)<0
        'WOW THIS IS WEIRD!'
    end
    %use the derivative (running diff) of k-distances to find discontinuity
    %smoothing is important to identify where the initial sustained rise occurs 
    kdist_diff = smoothdata(diff(kdist), 'gaussian', length(kdist)/10); 
    
    %% create a threshold for what steepness counts as real
    [~, max_diff_loc] = max(kdist_diff); 
    temp = kdist_diff; 
    %trim out the area around the max value (usually going to be the end)
    if max_diff_loc<6
        temp(1:10) = [];
    elseif max_diff_loc>length(kdist_diff)-11
       temp(length(kdist_diff)-10:end) = [];      
    else
        temp(max_diff_loc-5:max_diff_loc+5) = []; 
    end
    
    mean_diff = mean(temp); 
    std_diff = std(temp); 
    
    thresh = mean_diff + 2*std_diff; 
    
    %% find candidate inflection points where steepness is maximal
    %where does the second derivative cross 0 in the negative direction? 
    kdist_diff2 = diff(kdist_diff);
    kdist_diff2(kdist_diff2<0) = -1; 
    kdist_diff2(kdist_diff2>0) = 1; 
    candidates = find(diff(kdist_diff2)==-2) + 1;
    %down select to above threshold peaks in steepness
    candidates(kdist_diff(candidates)<thresh) = []; 
    %don't use a peak too near the end
    candidates(candidates>length(kdist)-5) = []; 
    
    %% figure out which steepness peak to use
    if isempty(candidates) 
        %check if the threshold was ever crossed
        if max(kdist_diff) > thresh
            %it is crossed, so just grab the point where it crosses and use it
            peakToUse = find(kdist_diff>thresh,1)-1; 
        else
            %try again without the smoothing to see if it can be salvaged
            kdist_diff = diff(kdist); 
            kdist_diff2 = diff(kdist_diff);
            kdist_diff2(kdist_diff2<0) = -1; 
            kdist_diff2(kdist_diff2>0) = 1; 
            candidates = find(diff(kdist_diff2)==-2) + 1;
            %down select to above threshold peaks in steepness
            candidates(kdist_diff(candidates)<thresh) = []; 
            %don't use a peak too near the end
            candidates(candidates>length(kdist)-5) = []; 
            if isempty(candidates)
                %if there is still no possible peak above threshold, then it
                %suggests that the input matrix is just a bunch of noise.
                %it might be possible to salvage by trying adding one or more 
                %rounds of extra distancing in the initial step above: 
                %distMat = pdist2(distMat, distMat) .^2; 
    
                peakToUse = []; 
            else
                peakToUse = candidates(1); 
            end
        end
    else
        %take the first peak as this will yield the most granular
        %clustering
        peakToUse = candidates(1); 
    end
    
    %% now do the clustering!
    if ~isempty(peakToUse)
        idxVals = dbscan(distMat, kdist(peakToUse), k, 'Distance','precomputed'); 
    else
        idxVals = zeros(length(distMat), 1);
        "warning: no good epsilon value found, returning zeros"
    end
end