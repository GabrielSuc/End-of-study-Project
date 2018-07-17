function [output_russell, flag] =  russell(z1,p1,k1,L,M,Nl,Nm,xin)

%Takes an IIR filter under the zero-pole-gain form (z,p,k) as input and 
%returns a decomposed version of it according to Russell's method


if (Nl + Nm) > length(p1)
    error('Error: cannot split denominator. Nl + Nm > Nr. poles')
end    
    
syms z k l A

Z1 = zeros(length(z1),1,'sym'); %zeros
P1 = zeros(length(p1),1,'sym'); %poles

H = zeros(length(z1),1,'sym');

%Factorization in a different way
for i =1:length(z1)
    
    
    Z1(i) = 1 - z1(i).*z^(-1);
    P1(i) = 1 - p1(i).*z^(-1);
    
   
    
    
end   
    
% Using the subtitution Russel does, we can create three different
% filters. 

%--------------------------------------------------------------------------
% FIRST FILTER - SPLIT THE POLES BY Nl - Allpolefilter indice L
%--------------------------------------------------------------------------
    
% Select first Nl poles
PoleL = zeros(Nl,1);

for i = 1:Nl
    PoleL(i) = p1(i)^L;
end


% Implementing the filter as SOS structure
Hl = dfilt.df2sos(zp2sos(zeros(Nl,1),PoleL,1));

%fvtool(Hl)

%--------------------------------------------------------------------------
% Second filter FILTER - FIR filter composed of three factors
%--------------------------------------------------------------------------

% We first need to identify the coeeficient ck described in the paper
% Mathematically, we have:

% Select first Nl poles
[Pl, B1] = deal(zeros(Nl,1,'sym'));

for i = 1:Nl
    Pl(i) = p1(i);
end

% Select last Nm poles
[Pm, B2] = deal(zeros(Nm,1,'sym'));

Np = Nl + Nm;

for i = (Nl + 1):Np
    Pm(i - Nl) = p1(i);
end

%------------------------------ First Factor ------------------------------

A1 = 0;

for k = 1:Nl
    for l = 0:(L-1)
        A1 = A1 + simplify((vpa(Pl(k))*z)^l);
    end
    B1(k) = A1;
    
    A1 = 0;
end


cumProdL = simplify(cumprod(vpa(B1)));
C1 = simplify(cumProdL(end),'Steps',150);


%----------------------------- Second Factor ------------------------------

A2 = 0;

for k = 1:Nm
    for m = 0:(M-1)
        A2 = A2 + simplify((vpa(Pm(k))*z)^m);
    end
    B2(k) = A2;
    
    A2 = 0;
end


cumProdM = simplify(cumprod(vpa(B2)));
C2 = simplify(cumProdM(end),'Steps',150);


%----------------------------- Third Factor ------------------------------



Z1 = zeros(length(z1),1,'sym');

for i =1:length(z1)
    
    Z1(i) = 1 - vpa(z1(i),2).*z^(1);

end 


C3 = simplify(k1.*prod(Z1),'Steps',150);


%-------------------------- Final FIR (numerator) -------------------------


numFIR = expand(simplify(C1.*C2.*C3,'Steps',150));

num  =  subs(numFIR,z,1/z); %Need to work with z^-1 %THIS IS THE REAL NUMERATOR


Nz = length(z1);

%Get the order of the numerator
Coeff_num = sym2poly(numFIR);


if length(Coeff_num) ~= (Nz + Nl*(L-1) + Nm*(M-1) + 1)
    error('Problem order numerator FIR')
end    

% This technique might sometimes not be possible, so we will have to perform
% a direct implementation of the filter

if L*M > length(Coeff_num)
    
    flag = 1; %create a flag
    output_russell = xin;
    return
end

%--------------------------------------------------------------------------
% <THIRD FILTER - SPLIT THE POLES BY Nm - Allpolefilter indice M
%--------------------------------------------------------------------------
    
% Select first Nm poles
PoleM = zeros(Nm,1);

for i = (Nl + 1):Np
    PoleM(i - Nl) = p1(i)^M;
end

% Implementing the filter as SOS structure
Hm = dfilt.df2sos(zp2sos(zeros(Nm,1),PoleM,1));


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                       Filtering the signal
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                           Through Hl(z)
%--------------------------------------------------------------------------


xout_filterL = filter(Hl,xin);


%--------------------------------------------------------------------------
%                           Through Hn(z)
%--------------------------------------------------------------------------

%We first need the coefficient a and b (the delays)

a = 0;
b = 0.1;

while(rem(b,1)~=0)
    a = a + 1;
    b = (L*a - 1)/M;
end   


disp('----------------------- Delay Coefficients ------------------------')
X = ['a = ', num2str(a), ' and b = ', num2str(b)];
disp(X)
disp('-------------------------------------------------------------------')

%Now need the LM polyphase components ek of hn
    
ek = myPolyphase(Coeff_num,1,L,M,'2');

%Have to filter xin through each branch 

sumBranch = 0;

for i = (L*M):-1:1 %Starting from the LM-1 branch
    %Creating the other branches before summation
    delayedBy_a = delayseq(xout_filterL,(i-1)*a);
    downsamp = downsample(delayedBy_a,M);
    filter_polyphase = filter(ek(i,:),1,downsamp);
    upsamp = upsample(filter_polyphase,L); 

    %Sum 

    sumBranch = sumBranch + upsamp;

    if i > 1
        sumBranch = delayseq(sumBranch,b);
    end

end

%--------------------------------------------------------------------------
%                           Through Hl(z)
%--------------------------------------------------------------------------

flag = 0; %no problem, no flag needed
output_russell = filter(Hm,sumBranch);

end