function [FM,FL] = getListStages(L,M)
    

%Returns list of stages FL and FM in order to perform multi_stage.m 
%We could also automotate the process by trying to find the best
%combination not only out of a predetermined number of stages, but out of 
%all of them. Might sometimes be excessively long to process (especially 
%high order P-M filters).
%Therefore, for the moment have to choose how many stages we desire to 
%treat manually. Nontheless, it matters the order of the stages so we could 
%still process the different cases for a given number of stages i.e.
%if FL = 10 2 8 or FL = 5 4 8, we won't have the same result in terms of MPOS 
%with the same FM = 7 7 3 for a 3-stage design for example. The following
%commented part is a begining of this automated procedure but might be to
%tricky to realize. 


% %Get factorization factors for L and M
% factL = factor(L);
% factM = factor(M);
% 
% length_factL = length(factL);
% length_factM = length(factM);
% 
% %Ask for how many stages we want to try
% %Max number of stages fixed to 6 because it's not necessary to go further
% %in the majority of the cases. Still can change it though.
% nbrStages = input('Choose number of stages desired (1-6): ');
% 
%     if nbrStages < 1 || nbrStages > 6 
%           error('Choose the number of stages accordingly');
%     end
% 
%     
% %Need to stuff ones if the number of stages we picked is bigger then the 
% %number of factorization factors. But we have to limit the identical number
% %of ones to one (we don't want to have a single stage where L = 1 and M = 1).
% 
% max_length_fact =  max(length(factL),length(factM));
% 
% if (nbrStages - max_length_fact) > 0 
%     
%     
%     
%     %METTRE UN ABS?
%     
%     
%     %We then fix the number of stages to prevent this case from happening
%     nbrStages = max_length_fact;
%     
%     factM = [factM, ones(1, nbrStages - length(factM))];
%     factL = [factL, ones(1, nbrStages - length(factL))];
% else    
% 
%     factM = [factM, ones(1, nbrStages - length(factM))];
%     factL = [factL, ones(1, nbrStages - length(factL))];
% 
% end
% 
% %If one among factL and/or factM is longer than the number of stages, have  
% %to limit it/them
% 
% if length(factL) > nbrStages
%     
%     %Generate list of indices
%     %All permutation are taking into account
%     
%     for i = 1:2*(length(factL)-nbrStages)
%         indices(i,:) = i:i+2*(length(factL)-nbrStages)-1;
%         indices_perm((i-1)*length(perms(indices(i,:))) + 1 ...
%             :(i-1)*length(perms(indices(i,:))) + length(perms(indices(i,:))),:) = perms(indices(i,:));
%     end    
%    
%     
%     %Have to create all possible subsets of factL that contains only
%     %nbrStages elements presently
%      
%     factL = repmat(factL,length(indices_perm),1);
%     
% 
%    for i = 1:length(indices_perm)
%        
%        
%        cumprod_first_indices = cumprod(factL(i,indices_perm(i,1:(length_factL-nbrStages))));
%        cumprod_second_indices = cumprod(factL(i,indices_perm(i,(length_factL-nbrStages+1):end)));
%        
%        %Replace by cumulative product of corresonding elements 
%        
%        factL(i,indices_perm(i,1)) = cumprod_first_indices(end);
%        factL(i,indices_perm(i,2:(length_factL-nbrStages))) = deal(0);
%        
%        factL(i,indices_perm(i,1 + length_factL - nbrStages)) = cumprod_second_indices(end);
%        factL(i,indices_perm(i,(length_factL - nbrStages + 2):end)) = deal(0);
%        
% %        for k = 1:length_factL
% %            if factL(i,k) == 0
% %                factL(i,k) = [];
% %            end    
% %        end
%     
%       factL(i,factL == 0) = [];  
%        
%    end    
%     
%    %Count the number of 0 in a line (same number for every line)
%    idx_zeros = factL(1,:) == 0;
%    nbr_zeros = sum(idx_zeros(:));
%    
%    %Shrink the different permutations
%    size_factL = size(factL);
%    
%       
%     
%     
%    factL = reshape(factL,[],size_factL(1))';
%    
% 
% end




%This is the manual part: a case has been made for each desired frequency. 
%It implies that we have to add a new elseif if we have new frequency to
%add. Moreover, the stage choices is personal. You can choose how many
%stages you want. Refer to first explanations at the begining of the file
%regarding some special cases where it might be a bit long to it by hand
%but still feasible.

    
if L == 160 && M == 147  % 44.1 -> 48 / 88.2 -> 96 / 176.4 -> 192 
    
    FM = [7,7,3,1];   
    FL = [4,4,5,2];   
    
elseif L == 2 && M == 1 % 44.1 -> 88.2 / 48 -> 96 / 88.2 -> 176.4     
    
    FM = 1;
    FL = 2;
    
elseif L == 320 && M == 147 % 44.1 -> 147 / 88.2 -> 192
    
    FM = [7,7,1,3,1];   
    FL = [5,4,2,4,2];
    
elseif L == 4 && M == 1 % 44.1 -> 176.4 / 48 -> 192
    
    FM = [1,1];
    FL = [2,2];

elseif L == 640 && M == 147 % 44.1 -> 192
    
    FM = [7,7,1,3];
    FL = [8,8,2,5];

elseif M == 160 && L == 147  % 44.1 <- 48 / 88.2 <- 96 / 176.4 <- 192 
    
    FL = [7,7,1,3];  
    FM = [5,4,4,2];
    
elseif M == 2 && L == 1 % 44.1 <- 88.2 / 48 <- 96 / 88.2 <- 176.4     
    
    FL = 1;
    FM = 2;
    
elseif M == 320 && L == 147 % 44.1 <- 147 / 88.2 <- 192
    
    FL = [];   
    FM = [];
    
elseif M == 4 && L == 1 % 44.1 <- 176.4 / 48 <- 192
    
    FL = [1,1];
    FM = [2,2];

elseif M == 640 && L == 147 % 44.1 <- 192
    
    FL = [7,7,1,3];
    FM = [8,8,2,5];    
end


end