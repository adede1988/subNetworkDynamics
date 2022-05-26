function [outSet] = formCleanIdx(curS)
        %useful function for transforming a set of cluster labels into the "canonical form"

        %input:
        %   curS: vector of cluter IDs where -1 indicates noise/
        %         unclustered and positive integer values indicate 
        %         different cluster labels

        %output:
        %   outSet: vector with the same length and information as the input
        %           but with labels transformed to meet Strehl and Ghosh's "canonical form" 

        %suggested by Strehl and Ghosh 2002: 
        %(i) LAMBDA1 = 1; (ii) for all i = 1, ... , n−1 : LAMBDAi+1 <= maxj=1, ... ,i(LAMBDAj)+1.
        %The first constraint enforces that the first object's label is cluster 1. The second constraint
        %assures that the cluster label LAMBDAi+1 of any successive object xi+1 either has a label that
        %occurred before or a label that is one greater than the highest label so far. By allowing
        %only representations that fulfill both constraints, the integer vector representation can be
        %forced to be unique. Transforming the labels into this `canonical form' solves the combining
        %problem if all clusterings are actually the same."

        %Adam Dede, adam.osman.dede@gmail.com, Fall 2021

        %Ref: Strehl, A., Ghosh, J. (2002). Cluster Ensembles--A knowledge reuse
        %framework for combining multiple partitions. J Machine Learning Research

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