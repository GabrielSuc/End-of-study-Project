function [ord1, ord2,F1,Fsc1] = Multistage_Proakis(Fx,Fp,Fsc,Rp,Rs)
%Return order of desired filter according to example of multistage decomposition 
%of chapter 11 of Proakis'book


%Defining the multistages (Decimator)
D1 = 21; D2 = 7;

%Filters' parameters
%Filter 1
F1 = Fx/D1;
Fsc1 = Fsc - F1;

deltaf1 = (Fsc1 -  Fp)/Fx;

%Ripples
delta11 = (10^(-Rp))/2;
delta21 = (10^(-Rs));

ord1 = round((-10*log10(delta11*delta21)-13)/(14.6*deltaf1)+1);


ord2 = 0;

end