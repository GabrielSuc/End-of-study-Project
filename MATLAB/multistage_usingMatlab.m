function Order = multistage_usingMatlab(L,M,Fsin,Fsout,Fp,Rp,Rs)

%Get the list of factors L & M for the different stages
[FM,FL] = getListStages(L,M);


%Ask for the type of filter we want to be returned
manual = input('1 to return a Parks-McClellan multistage filter, 2 for Elliptic, 3 for combination:');

    if isempty(manual)
          disp('Choose a filter');
    end
    
%Have to inverse Fx and Fy considering the overall system as a decimator or an interpolator
if M > L
     [Fsout, Fsin] = deal(Fsin,Fsout);
end

% if manual ==1 
% 
% % --------------------------- Parks-McClellan -----------------------------
% 
% %Need of specific variables for the design
% a = [1 0];
% f = zeros(length(FM),length(a));
% 
% 
% %Ripples
% Delta1 = (10^((Rp/length(FM))/20)); %Reduce the passband ripple by factor of length(FM) so that passband
% %ripple in the cascade of the length(FM) filters doesn't exeed Delta1
% Delta2 = (10^(-Rs/20));           
% 
% 
% dev = [abs(Delta1 - 1)/(Delta1 + 1) Delta2]; 
% % -------------------------- Filters' parameters --------------------------
% 
% %Band-edge frequencies
% Fp = zeros(1,length(FM));
% Fs = zeros(1,length(FM)+1); %Or FL, same length 
% %Orders
% Order = zeros(length(FM),1);  
% %Final Filter
% H = cell(1,length(FM));
% 
% %First sampling frequency before any stages
% Fs(1) = Fsin;
% 
% for i = 2:(length(FM)+1) 
%     
%     %Defining the others sampling frequencies, after crossing a stage
%     Fs(i) = (Fs(i-1)*FL(i-1))/FM(i-1);
%     Fp(i-1) = (Fs(i-1)/2)*min(1/FL(i-1),1/FM(i-1));
%     
%     %Defining the limit frequencies for the design
%     f(i-1,:) = [Fp(i-1) Fs(i-1)/2];
%     %Creating the P-M filters
%     [Order(i-1),fo,ao,w] = firpmord(f(i-1,:),a,dev,2*pi*Fs(i-1));
%     
%     if Order(i-1) < 3
%         Order(i-1) = 3;
%     end
%     
%     H{i-1} = dfilt.dffirt(FL(i-1)*(firpm(Order(i-1),fo,ao,w))); %One implementation among others dfilt.dffirt
%     
% %     [H_freq,freq(i-1,:)] = freqz(H{i-1},1,250,Fs(i-1));
% %     A(i-1,:) = abs(H_freq);
% end
%     
% % plot (freq',A','LineWidth',1.5);
% 
% %Cascading them
% H_multi = dfilt.cascade(H{1,1:end});
% return




end