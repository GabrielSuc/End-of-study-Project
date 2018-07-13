function Order = multistage_turek(L,M,Fs,Fp,Fsc,Rp,Rs)
%Multistage interpolator based on MIT Turek's thesis


%Get the list of factors L & M for the different stages
[FM,FL] = getListStages(L,M);

%Defining cutoff frequencies
[Fpassband,Fstopband,Order] = deal(zeros(1,length(FL)));

%First frequencies
Fpassband(1,1) = Fp/FL(1);
Fstopband(1,1) = Fsc/FL(1);
    
for i = 2:length(FL)
    
    A = cumprod(FL(1:i-1));
    B = cumprod(FL(1:i));
    
    Fpassband(1,i) = Fp/B(end);
    Fstopband(1,i) = (Fs*A(end) - Fsc)/B(end); %Fs = input sampling frequency i.e. Fx

end


% --------------------------- Parks-McClellan -----------------------------

%Need of specific variables for the design
a = [1 0];
%Ripples
Delta1 = 10^(Rp/20); %Reduce the passband ripple by factor of length(FM) so that passband
%ripple in the cascade of the length(FM) filters doesn't exeed Delta1
Delta2 = 10^(-Rs/20);           


dev = [(Delta1 - 1)/(length(FL)*(Delta1 + 1)) Delta2]; %abs(Delta1 - 1)

for i = 1:length(FL)
   
    Fmax = Fs;
    
    %Defining the limit frequencies for the design
    f = [Fpassband(i)/Fmax Fstopband(i)/Fmax];
    
    %Creating the P-M filters
    [Order(i),fo,ao,w] = firpmord(f,a,dev);
    H = dfilt.dfsymfir(firpm(Order(i),fo,ao,w)); %Folded implementation, reduces multiples by 2
    fvtool(H)
    
    
end    

end