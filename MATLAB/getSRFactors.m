function [L,M,FL,FM] = getSRFactors(Fx,Fy)

% Function that returns the factors L & M needed for fractional SRC
% Fy = L/M Fx


[L, M] = rat(Fy/Fx); 


end