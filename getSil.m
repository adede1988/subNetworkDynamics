function [silhouettes, clusti] = getSil(corMat, idx)  
    
    %This function calculates the silouette values for every input variable
    %in a symmetric input matrix (corMat) based on the cluster scheme
    %specified by an integer vector (idx)
    %it returns both the silouette values themselves (silhouettes) and some
    %useful indexing information to help keep track of silouette values
    %(clusti). Silhouette values are returned sorted into clusters and then
    %within clusters in decending order. This facilitates plotting for
    %quick assessment of cluster stability. 

    %That is: plot(getSil(corMat,idx)) is very useful and quickly
    %diagnostic

    %Input: corMat: variable X variable symmetric similarity matrix. This
    %       will work with R squared, covariance, distance, etc. 
    %       idx: cluster labels as integers, length must be same as number
    %       of electrodes, noise = -1, true clusters must be positive
    %       values
    
    %output: silhouettes: 1d vector of silhouette value for each electrode 
    %                       in a cluster, returns zeros for electrodes in 
    %                       the "noise" cluster, sorted by cluster, and
    %                       then within cluster by value. Noise cluster is
    %                       at the end of the vector
    %        clusti: electrode X 2 matrix, first column is electrode
    %                       indicies based on the order from the input, 
    %                       second column is cluster ID from idx input. 
    %                       sorted order is the same as silhouettes
  
    %Adam Dede, adam.osman.dede@gmail.com, Fall 2021
      
    %Ref: Rousseeuw, P.J. (1987). Silhouettes: a graphical aid to the
    %interpretation and validation of cluster analysis. J Comp and Applied
    %Mathematics.


      A = zeros(length(corMat), 1); 
      clusti = zeros(length(corMat), 2); 
      ai = 1; 
      B = zeros(length(corMat), 1); 
      bi = 1;
      clustIDs = unique(idx); 
      noisei = find(idx==-1); 
      clustIDs(clustIDs==-1) = []; %remove noise from silhouette calculation!
      if length(clustIDs)>1
      for ii = 1:length(clustIDs)
         curClus = corMat(idx==clustIDs(ii),:); 
         curClusDists = pdist2(curClus, curClus); 
         curClusDists = mean(curClusDists); 
         A(ai:ai+sum(idx==clustIDs(ii))-1) = curClusDists; 
         clusti(ai:ai+sum(idx==clustIDs(ii))-1,1) = find(idx==clustIDs(ii)); 
         clusti(ai:ai+sum(idx==clustIDs(ii))-1,2) = clustIDs(ii); 
         ai = ai + sum(idx==clustIDs(ii)); 
         

         %need to get the average distance of each point (pp) to all points in
         %other cluster (oc), do for all clusters, then take minimum of those
         %averages
         otherIDs = clustIDs(clustIDs ~= clustIDs(ii)); 
         for pp = 1:size(curClus,1)
             posB = zeros(length(otherIDs), 1); 
             for oc =1:length(otherIDs) 
                curOth = corMat(idx==otherIDs(oc), :); 
                posB(oc) = mean(pdist2(curOth, curClus(pp,:)));  
             end
             B(bi) = min(posB); 
             bi = bi+1; 
         end
          
      end
      clusti(ai:end,1) = noisei; 
      clusti(ai:end,2) = -1; 
      nomin = B - A; 
      denom = max([B,A]')'; 
      silhouettes = nomin ./ denom; 
      %make the clusti ordering go in order of silhouette value
      for ii = 1:length(clustIDs)
          curSil = silhouettes(clusti(:,2) == clustIDs(ii)); 
          curi = clusti(clusti(:,2) == clustIDs(ii), 1); 
          [silhouettes(clusti(:,2) == clustIDs(ii)), order] = sort(curSil, 'descend'); 
          clusti(clusti(:,2) == clustIDs(ii), 1) = curi(order); 
      end
      
      
      else
%         warning('only 1 or 0 clusters in input, returning NaNs')
        silhouettes = nan(length(corMat), 1); 
        clusti = [1:1:length(corMat)]; 
      end

end