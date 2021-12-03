function [L, comps, maps, W, result] = GEDclust(varargin)   

%This function performs generalized eigen decomposition (GED) in order to
%find a single time series that maximizes the difference between Sinput and
%Rinput. Rather than calculating the covariance matrices for input to GED
%directly from the entire data series, both Sinput and Rinput are cut into
%2000 point snips and covariance is calculated for each snip. These
%covariance values are then cleaned by removing outlier covariance matrices
%and averaged together. This procedure avoids the analysis being thrown by 
%transient events. Summary statistics, the eigenvalues, and the analytic
%signal of the eigencomponents are returned.

%This function is written with flexibility to work either as a stand alone
%or in a clustering pipeline. In a clustering pipeline, the inputs CurClus
%and clustMap determine which subset of channels to use in the Sinput and
%Rinput matrices. This allows one to loop across clusters dynamically. 

%if using this function as a stand alone, simply leave off curClus and
%clustMap inputs and the entire Sinput and Rinput matrices will be used in
%the GED


%inputs: 
%   Sinput:              channels X timepoints OR timepoints X channels of 
%                        signal
%   Rinput:              channels X timepoints OR timepoints X channels of 
%                        reference
%   onsets (optional):   vector of timepoints at which to slice the data.
%                        default=[1:2000:length(timepoints)]
%   curClus (optional):  integer label of the current cluster to use.
%                        default=1
%   clustMap (optional): vector of cluster labels the same length as 
%                        channels. default = ones(length(channels),1)
%   IDcode (optional):   optional identifying value. Useful for keeping
%                        track of analyses. default = 1
%   bandFrex (optional): optional 2-item vector of lower and upper
%                        frequency bounds used to create Sinput filtered 
%                        signal. Used exclusively for saving descriptive
%                        information in the output results, not for any
%                        calculation, so it really doesn't matter. Useful
%                        if Sinput is a bandpassed version of Rinput. This
%                        is a good spot to record what band was used.
%                        default = [1,2]

%outputs: 
%   L:      eigenvalues 
%   comps:  channels X timepoints eigencomponents resulting from GED
%   maps:   channels X channels component maps for each eigencomponent
%   W:      eigenvectors
%   result: summary statistics about the GED returned as a vector where: 
            %item 1: eigenvalue of top component
            %item 2: IDcode
            %item 3: channel count
            %item 4: eigenValue of 2nd component 
            %item 5: min frequency
            %item 6: max frequency 
            %item 7: mean frequency
            %item 8: sum of eigenValues
            %item 9: rank of S matrix
            %item 10: rank of R matrix
            %item 11: negative eigenvec flag
            %item 12: variance explained by 1st component


%Adam Dede, adam.osman.dede@gmail.com, Fall 2021

switch nargin
    case 1
        warning('Error: needs at least 2 inputs: signal and reference matrices')
        return
    case 2
        Sinput = varargin{1}; 
        Rinput = varargin{2};
        onsets = [1:2000:max(size(Sinput))]; 
        curClus = 1; 
        clustMap = ones(min(size(Sinput)),1);
        IDcode = 1; 
        bandFrex = [1,2]; 
    case 3
        Sinput = varargin{1}; 
        Rinput = varargin{2};
        onsets = varargin{3}; 
        curClus = 1; 
        clustMap = ones(min(size(Sinput)),1);
        IDcode = 1; 
        bandFrex = [1,2]; 
    case 4
        warning('Error: if specifying curClus, then clustMap must also be specified')
        return
    case 5
        Sinput = varargin{1}; 
        Rinput = varargin{2};
        onsets = varargin{3}; 
        curClus = varargin{4}; 
        clustMap = varargin{5};
        IDcode = 1; 
        bandFrex = [1,2];
    case 6
        Sinput = varargin{1}; 
        Rinput = varargin{2};
        onsets = varargin{3}; 
        curClus = varargin{4}; 
        clustMap = varargin{5};
        IDcode = varargin{6}; 
        bandFrex = [1,2];
    case 7
        Sinput = varargin{1}; 
        Rinput = varargin{2};
        onsets = varargin{3}; 
        curClus = varargin{4}; 
        clustMap = varargin{5};
        IDcode = varargin{6}; 
        bandFrex = varargin{7};
    otherwise
        warning('Error: needs at least 2 inputs: signal and reference matrices')
        return
end
if size(Sinput,1) > size(Sinput,2) % if input is time X chan, switch it
    Sinput = Sinput'; 
end

if size(Rinput,1) > size(Rinput,2) % if input is time X chan, switch it
    Rinput = Rinput'; 
end
%  Sinput, Rinput, onsets, curClus, clustMap, IDcode, bandFrex

   
    curChan_i = find(clustMap==curClus);


         % tensor S and R
    Stmp = zeros(length(onsets)-1,length(curChan_i),length(curChan_i));
    Rtmp = Stmp; 

    for segi=1:length(onsets)-1

        snipdat = ( Sinput(curChan_i,onsets(segi):(onsets(segi+1)-1)) );
        rnipdat = ( Rinput(curChan_i,onsets(segi):(onsets(segi+1)-1)) );

        snipdat = snipdat - mean(snipdat,2);
        rnipdat = rnipdat - mean(rnipdat,2);

        Stmp(segi,:,:) = snipdat*snipdat'/(onsets(segi+1) - onsets(segi));
        Rtmp(segi,:,:) = rnipdat*rnipdat'/(onsets(segi+1) - onsets(segi)); 
    end
    % clean S
    meanS = squeeze(mean(Stmp));
    dists = zeros(1,size(Stmp,1));
    for segi=1:size(Stmp,1)
        s = Stmp(segi,:,:);
        dists(segi) = sqrt( sum((s(:)-meanS(:)).^2) );
    end
    %S for GED
    S = squeeze(mean( Stmp(zscore(dists)<3,:,:) ,1));


    % clean R
    meanR = squeeze(mean(Rtmp));
    dists = zeros(1,size(Rtmp,1));
    for segi=1:size(Rtmp,1)
        r = Rtmp(segi,:,:);
        dists(segi) = sqrt( sum((r(:)-meanR(:)).^2) );
    end
    %R for GED
    R = squeeze(mean( Rtmp(zscore(dists)<3,:,:) ,1));

    % regularize R
    gamma = .01;
    Rr = R *(1-gamma) + eye(length(R))*gamma*mean(eig(R));

    % global variance normalize
    S = S / (std(S(:))/std(R(:)));
    
    % eig
    [W,L] = eig(S,Rr);

    % sort and store
    [L,sidx] = sort(diag(L),'descend');
    W = W(:,sidx); 

    %get the results 
    %item 1: eigenvalue of top component
    %item 2: IDcode
    %item 3: channel count
    %item 4: eigenValue of 2nd component 
    %item 5: min frequency
    %item 6: max frequency 
    %item 7: mean frequency
    %item 8: sum of eigenValues
    %item 9: rank of S matrix
    %item 10: rank of R matrix
    %item 11: negative eigenvec flag
    %item 12: variance explained by top component
    result = zeros(12,1); 
    result(1) = L(1); 
    result(2) = IDcode; 
    result(3) = length(curChan_i); 
    result(4) = L(2); 
    result(5) = bandFrex(1); 
    result(6) = bandFrex(2); 
    result(7) = mean(bandFrex); 
    result(8) = sum(abs(L)); 
    result(9) = rank(S); 
    result(10) = rank(R); 
    if sum(L<0)>0
        result(11) = 1;
    else
        result(11) = 0; 
    end
    result(12) = result(1) / result(8); 

    comps = hilbert(W' * Sinput(curChan_i,:));
    maps = W' * S; 

end