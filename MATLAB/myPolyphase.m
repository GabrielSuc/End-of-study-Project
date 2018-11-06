function polyMatrix = myPolyphase(h,L,M)


% Returns the type 1 polyphase components (based on cf. Vaidyanathan - "Multirate
% Systems and Filter Banks" and adapted by Andrew I. Russell) of the filter's coefficients h, i.e. 
% ek[n] = h[LMn + k] with k in {0,...,LM-1}

% First off, we need to zero-pad h to make its length equal to an
% integer multiple of L*M

N = length(h);
q = fix(N/(L*M)); %Quotient 

if rem(N,L*M) ~=0  
    h = [h, zeros(1, (q+1)*L*M - N)];
end

% Length of polyphase filters
len_poly = ceil(length(h)/(L*M));

% Defining the size of the matrix
polyMatrix = zeros(L*M,len_poly);

% Creating the polyphase components
for i = 1:(L*M)
    polyMatrix(i,:) = h((0+i):L*M:end);
end    

