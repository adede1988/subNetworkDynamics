function goodi = covMatCleaning(varargin)

%This function calculates indices of good data based on covariance matrix
%cleaning. Timeseries data are taken in as either channels X time OR time X
%channels. Data are chopped into 2000 point windows with a window increment
%of 100 points. The channel X channel covariance matrices for each window
%are calcualted. The average covariance matrix is calculated. The euclidean
%distance between each of the window covariance matrices and the average is
%calculated resulting in a vector of distances. These distances are
%z-scored. Any data points that contributed to a window who's covariance
%z-scored distance is beyond a threshold are flagged as noise. 

%inputs: 
    %data:               Chan X time OR time X chan matrix of timeseries data
    %zThresh (optional): threshold for z-scores such that values beyond the
    %                    threshold will be flagged as noise. Default = 3
    %plotIt (optional):  if plotIt==1 then plot some summary info else
    %                    don't. Default = 0

%output: 
    %goodi:              vector with same length as long dimension of data
    %                    input. 1s indicate points below zThresh. 0s
    %                    indicate points beyond zThresh. 

%Adam Dede, adam.osman.dede@gmail.com, Fall 2021


switch nargin
    case 1
        data = varargin{1}; 
        zThresh = 3;
        plotIt = 0; 
    case 2
        data = varargin{1}; 
        zThresh = varargin{2};
        plotIt = 0; 
    case 3
        data = varargin{1}; 
        zThresh = varargin{2};
        plotIt = varargin{3};
    otherwise
        warning('Error: at least one input is needed')
        return
end

if size(data,1) > size(data,2) % if input is time X chan, switch it
    data = data'; 
end

%create start and end indices for a 2000ms window that slides in 100ms
%steps
windowStarts = [1:100:size(data,2)]; 
windowEnds = windowStarts(21:end); 
windowStarts = windowStarts(1:end-20); 

%preallocate for the covariance matrices
covMats = zeros(size(data,1), size(data,1), length(windowEnds));  

%get cov matrices in every window
for win = 1:length(windowStarts)
   covMats(:,:,win) = (data(:,windowStarts(win):windowEnds(win)) * ...
                       data(:,windowStarts(win):windowEnds(win))') ./ 2000;
end

%get an average covMat
covMatAvg = mean(covMats,3); 

%find the sum of squared difference of each epoch from the average
diffs = zeros(length(windowStarts),1);
for win = 1:length(windowStarts)
   cur = covMats(:,:,win);  
   diffs(win) = sum((cur(:) - covMatAvg(:)).^2, 'omitnan');   
end

% z-score the SS differences 
diffsZ = (diffs - mean(diffs( diffs<prctile(diffs,97.5) & diffs>prctile(diffs,2.5)))) /...
    std(diffs( diffs<prctile(diffs,97.5) & diffs>prctile(diffs,2.5))); 

%1=window that starts at the corresponding windowStarts index is good
%0=window that starts at the corresponding windowStarts index is bad
goodWindows = ones(length(windowStarts), 1);
%z-score threshold 
goodWindows(diffsZ>zThresh) = 0; 
goodWindows(diffsZ<-zThresh) = 0; 
%translate back out into indices across individual data points
goodi = ones(size(data,2),1);
badWins = find(goodWindows==0); 
for bw = 1:length(badWins)
   goodi(windowStarts(badWins(bw)):windowEnds(badWins(bw))) = 0;  
end


if plotIt==1
    figure
    subplot(221)
    plot(diffsZ)
    if max(diffsZ) > zThresh
        maxy = max(diffsZ) + .3; 
    else
        maxy = zThresh + .3; 
    end
    if min(diffsZ) < -zThresh
        miny = min(diffsZ) - .3; 
    else 
        miny = -zThresh - .3; 
    end
    ylim([miny, maxy])
    yline(zThresh, '--', 'lineWidth', 3, 'Alpha', .5)
    yline(-zThresh, '--', 'lineWidth', 3, 'Alpha', .5)
    title('epoch differences from average (z-score)')

    subplot(222)
    imagesc(covMatAvg)
    title('average covariance matrix')

    if length(badWins)>0
        subplot(223)
        imagesc(covMats(:,:,badWins(1)))
        title('example high z-score cov mat')
    end
    if length(badWins)>1
        subplot(224)
        imagesc(covMats(:,:,badWins(end)))
        title('another example high z-score')
    end

end

