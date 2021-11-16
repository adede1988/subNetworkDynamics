function [nmi] = nmi(set1, set2)

    clustIDs = unique(set1); 
    clustIDs2 = unique(set2); 
    clustIDs(clustIDs==-1) = []; 
    clustIDs2(clustIDs2==-1) = []; 
    %calcualte the normalized mutual information
    nominSum = 0; 
    denomh = 0; 
    denoml = 0; 
    n = length(set1);
    for c = 1:length(clustIDs)
        vec1 = zeros(length(set1),1); 
        vec1(set1==clustIDs(c)) = 1; 
        nh = sum(vec1); 
        for cc = 1:length(clustIDs2)
            vec2 = zeros(length(set2),1);
            vec2(set2==clustIDs2(cc)) = 1; 

            nl = sum(vec2); 
            vec3 = zeros(length(set1),1); 
            vec3(vec1==1 & vec2==1) = 1; 
            nhl = sum(vec3); 
            if nhl==0 %correction for zero overlap!
                nhl=.000000001;
            end

            nominSum = nominSum + nhl*log10((n*nhl)/(nh*nl)); 

        end
        denomh = denomh + nh*log10(nh/n); 
    end
    for cc = 1:length(clustIDs2)
            vec2 = zeros(length(set1),1);
            vec2(set2==clustIDs2(cc)) = 1; 
            nl = sum(vec2); 
            denoml = denoml + nl*log10(nl/n); 
    end

    nmi = nominSum / sqrt(denomh*denoml);
    
    %under and over situations will only occur when there's only a noise
    %cluster in one or other of the cluster schemes being considered!
    if nmi < 0 
%         "under!"
%         length(clustIDs)
%         length(clustIDs2)
        nmi = 0; 
    elseif nmi > 1
%         "over!"
%         length(clustIDs)
%         length(clustIDs2)
        nmi = 1; 
    end

end