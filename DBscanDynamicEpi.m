function [idxVals] = DBscanDynamicEpi(varargin)

%This function is a wrapper for Matlab's builtin DBscan function
%There must be at least one input in the first input position: 
%1) corMat   : a symmetric matrix of correlations between points to be
%              clustered (also works with any symmetric similarity matrix)
%There are four optional inputs (which are assumed to be inputs 2-5): 
%2) k        : The k parameter for input to DBscan [default = 3]
%3) transform: how the input corMat should be transformed before clustering
%              [default = 'raw']
%4) d        : the number of extra loops of distancing to perform; only
%              works with 'mapDistx' as the transformation; [default = 1] 
%5) plotIt   : 0 = do not make summary plots; 1 = make summary plots 
%              [default = 0]

%four main computations are performed: 
%1) transformation of the corMat. The goal is to increase the contrast of
%the input matrix to make for better clustering output. The possible
%transformations are specified as string (or numeric) inputs to the transform varargin 
%as follows:
%  1 OR 'raw'[default]: pdist2(corMat, corMat, 'correlation')
%  2 OR 'map'         : pdist2(corrcoef(corMat), corrcoef(corMat), 'correlation') 
%  3 OR 'mapDist'     : pdist2(corrcoef(corMat), corrcoef(corMat))
%  4 OR 'mapDist2'    : pdist2(corrcoef(corMat), corrcoef(corMat)) .^2
%  5 OR 'mapDistX'    : loops d additional times on the output from the
%                       mapDist2 option
%2) finding of appropriate epsilon value. Using the k-distance values in
%the transformed matrix, an epsilon value is selected that maximizes the
%steepness of the sorted k-distance curve. See discussion of Figure 4 in
%Ester et al. (1986) for more detail. 
%3) application of DBscan 
%4) plotting if requested

%output is a 1-d vector of cluster assignments where -1 means unclustered and
%positive integers indicate different clusters

%the general idea of transforming the space to be clustered from one of
%point-wise connections/similarities/correlations to one of row by column
%connections using corrcoef(corMat) was inspired by Liu et al. (2012).

%NOTE: the way the threshold for step 2 (finding epsilon) is decided effectively 
%assumes that there are 30-100 items to be clustered. This is because it 
%disregards the points that are +-5 positions around the maximum in the 
%derivative of the k-distance curve. For clustering of larger data sets
%this value should be increased. Key edit at lines 154-160

%Adam Dede, adam.osman.dede@gmail.com, Fall 2021

%Refs: 
%Liu, Zhu, Chen, 2012. J neurosci Methods.
%Ester, Kriegel, Sander, Xu, 1996

    %% get out the inputs
    
    switch nargin
        case 1
            corMat = varargin{1}; 
            k = 3;
            transform = 'raw'; 
            d = 1; 
            plotIt = 0;
        case 2
            corMat = varargin{1}; 
            k = varargin{2}; 
            transform = 'raw';
            d = 1; 
            plotIt = 0; 
        case 3
            corMat = varargin{1}; 
            k = varargin{2}; 
            transform = varargin{3};
            d = 1; 
            plotIt = 0; 
        case 4
            corMat = varargin{1}; 
            k = varargin{2}; 
            transform = varargin{3};
            d = varargin{4}; 
            plotIt = 0; 
        case 5
            corMat = varargin{1}; 
            k = varargin{2}; 
            transform = varargin{3};
            d = varargin{4}; 
            plotIt = varargin{5}; 
        otherwise
            warning('Error: at least one input is needed')
            return
    end

    if isnumeric(transform)
        if transform==1
            transform = 'raw';
        elseif transform==2
            transform = 'map'; 
        elseif transform==3
            transform = 'mapDist';
        elseif transform ==4
            transform = 'mapDist2';
        elseif transform > 4
            transform = 'mapDistX';
        end
    end
    if isnumeric(transform)
        warning('transform argument must be a positive integer or a string')
        return
    end
    if isempty(k) || ~isnumeric(k)
        k = 3;
    end
    
    %% transformation step
    % various ways of increasing contrast on the input correlationo matrix
    % convert into a correlation matrix of relationships between connectivity
    % maps rather than between individual channels is inspired by Liu et al., 2012
    switch transform
        case 'raw'
            distMat = pdist2(corMat, corMat, 'correlation');
        case 'map'        
            distMat = pdist2(corrcoef(corMat), corrcoef(corMat), 'correlation'); 
        case 'mapDist'     
            distMat = pdist2(corrcoef(corMat), corrcoef(corMat));
        case 'mapDist2'     
            distMat = pdist2(corrcoef(corMat), corrcoef(corMat)) .^2;
        case 'mapDistX'
            distMat = pdist2(corrcoef(corMat), corrcoef(corMat)) .^2;
            for ii = 1:d
                distMat = pdist2(distMat, distMat).^2; 
            end
    end

    
    %% find appropriate epsilon value: 
        %f
    %% find the k-distance for each point
    kdist = zeros(length(distMat),1);
    for ii = 1:length(kdist)
        curRow = distMat(ii,:); 
        curRow = sort(curRow); 
        kdist(ii) = curRow(k); 
    end
    kdist = sort(kdist); 

    %use the derivative (running diff) of k-distances to find discontinuity
    %smoothing is important to identify where the initial sustained rise occurs 
    kdist_diff = smoothdata(diff(kdist), 'gaussian', length(kdist)/10); 
    
    %% create a threshold for what steepness counts as real
    % KNOWN BUG: this is really best for clustering about 60-100 items.
    % Larger sets won't work so well because of how this threshold is set
    %           Attempted fix: use of trim variable makes cutting off
    %           around the max dynamic
    [~, max_diff_loc] = max(kdist_diff); 
    temp = kdist_diff; 
    trim = round(length(temp)/10); 
    %trim out the area around the max value (usually going to be the end)
    if max_diff_loc<6
        temp(1:trim) = [];
    elseif max_diff_loc>length(kdist_diff)-11
       temp(length(kdist_diff)-trim:end) = [];      
    else
        temp(max_diff_loc-ceil(trim/2):max_diff_loc+floor(trim/2)) = []; 
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
            
            if peakToUse == 0 && max(kdist_diff(2:end)) > thresh
                peakToUse = find(kdist_diff(2:end)>thresh,1);
            elseif peakToUse == 0
                peakToUse = 1; 
            end 
                
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
%       epsis(epsiidx)
%         idxVals = dbscan(distMat, epsis(epsiidx), k, 'Distance','precomputed'); 
%         kdist(peakToUse)
        idxVals = dbscan(distMat, kdist(peakToUse), k, 'Distance','precomputed'); 


        %% plot summary of how epsilon was chosen? 
        if plotIt ==1
            figure
            subplot(231)
            plot(kdist)
            hold on
            yline(kdist(peakToUse))
            xline(peakToUse)
            
            subplot(232)
            plot(kdist_diff)
            hold on 
            yline(thresh, '--', 'LineWidth', 3, 'alpha', .25)
            xline(peakToUse)
        
            subplot(233)
            imagesc(corMat)
        
            subplot(234)
            imagesc(distMat)
        
        
            subplot(235)
            [~, order] = sort(idxVals);
            imagesc(corMat(order, order))
            locs = find(diff(sort(idxVals)));
            for val=1:length(locs)
                yline(locs(val)+.5, '--', 'LineWidth', 3, 'alpha', .5, 'color', 'k')
                xline(locs(val)+.5, '--', 'LineWidth', 3, 'alpha', .5, 'color', 'k')
            end

            subplot(236)
            imagesc(distMat(order, order))
            for val=1:length(locs)
                yline(locs(val)+.5, '--', 'LineWidth', 3, 'alpha', .5, 'color', 'k')
                xline(locs(val)+.5, '--', 'LineWidth', 3, 'alpha', .5, 'color', 'k')
            end
        
        
        
        
        end


    else
        idxVals = ones(length(distMat), 1) * -1;
%         "warning: no good epsilon value found, returning -1"
%         d
    end
end