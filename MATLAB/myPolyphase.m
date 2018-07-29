function PolyMatrix = myPolyphase(b,a,L,M,str)

% Returns the k polyphase matrix PolyMatrix. 

% Depends on the value of str i.e. if str == '1', case seen in Oppenheim 
% online course (moving the recursive part before the expander). Have to 
% choose adequately a & b according to whether we
% are using a FIR or IIR filter. Here, Ek are the z-transform of the
% polyphase components ek = h[Mn + k], for k = 0 ... M-1.
% If str == '2', we are in the case of Russel's paper. Now ek[n] = h[LMn + k],
% for k = 0 ... LM-1.


if str == '1' %FIR/IIR case
    if a == 1 %FIR 
       N = length(b);
       if N > L %Just to make sure
        if rem(N,L) ~= 0
                if L*ceil(N/(L)) > N
                    PolyMatrix = vec2mat([b;zeros(L*ceil(N/(L)) - N,1)],ceil(N/(L)));
                else 
                    PolyMatrix = vec2mat([b;zeros((ceil(N/(L))+1)*L - N,1)],ceil(N/(L)));
                end    
         else                                                                                                      
                PolyMatrix = vec2mat(b,fix(N/(L)));
         end
       else
       error('Length of filter is smaller than L')   
       end 
    else      % IIR    
        error('Have to be a FIR filter')
    end

elseif str == '2' %Russel's case
    %We now move the expander and the decimator the other side of the
    %polyphase components i.e. there are LM polyphase components
    
    if a == 1 
        N = length(b);
        if N >= L*M
            if rem(N,L*M) ~= 0
                if L*M*ceil(N/(L*M)) > N
                    PolyMatrix = vec2mat([b,zeros(1,L*M*ceil(N/(L*M)) - N)],ceil(N/(L*M)));
                else 
                    PolyMatrix = vec2mat([b,zeros(1,(ceil(N/(L*M))+1)*L*M - N)],ceil(N/(L*M)));
                end    
            else                                                                                                      
                PolyMatrix = vec2mat(b,fix(N/(L*M)));
            end
        else
            error('Length of filter is smaller than L*M')
        end  
    end
else 
    error('Select a correct value for str: 1 for FIR/IIR case, 2 for Russells decomposition ')
end    

end