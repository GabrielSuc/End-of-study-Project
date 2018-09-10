function [output_russell, flag] =  russell(z1,p1,k1,L,M,b_fir,Nl,Nm,xin)

%Takes an IIR filter under the zero-pole-gain form (z,p,k) as input and 
%returns a decomposed version of it according to Russell's method

Np = length(p1);
Nz = length(z1);
%--------------------------------------------------------------------------
                            %First Filter
%--------------------------------------------------------------------------

%Select first Nl poles
poleL = zeros(Nl,1);

for i = 1:Nl
    poleL(i) = p1(i)^L;
end
 

%Implementing the filter as SOS structure
format long e

Hl = dfilt.df2sos(zp2sos(zeros(Nl,1),poleL,1));
fvtool(Hl)


%--------------------------------------------------------------------------
                            %Third Filter
%--------------------------------------------------------------------------

%Select first Nl poles
poleM = zeros(Nm,1);

for i = (Nl + 1):Np
    poleM(i - Nl) = p1(i)^M;
end
 

%Implementing the filter as SOS structure
format long e

Hm = dfilt.df2sos(zp2sos(zeros(Nm,1),poleM,k1));
%[bm,am] =zp2tf(zeros(Nm,1),poleM,1);
fvtool(Hm);
%--------------------------------------------------------------------------
                            %Second Filter
%--------------------------------------------------------------------------


%Convolution of the lists of the pole raised to the power of L

listL = zeros(Nl,1);

for i = 1:Nl
    listL(i) = p1(i);
end

if ~isempty(listL)

    coeffL = power(listL(1),[0:(L-1)]);

    for i = 2:Nl
        coeffL = conv(coeffL, power(listL(i),[0:(L-1)]));
    end
    
else
    coeffL = 1;
end

%Convolution of the lists of the pole raised to the power of M

listM = zeros(Nm,1);

for i = (Nl + 1):Np
    listM(i - Nl) = p1(i);
end

coeffM = power(listM(1),[0:(M-1)]);

for i = 2:Nm
    coeffM = conv(coeffM, power(listM(i),[0:(M-1)]));
end

%Convolution of all the coefficients

num_fir = real(conv(conv(coeffL,coeffM),b_fir));


if length(num_fir) ~= (Nz + Nl*(L-1) + Nm*(M-1) + 1)
    error('Problem order numerator FIR')
end 

Hn = dfilt.dffir(num_fir);

fvtool(Hn)




% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
%                       Filtering the signal
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% 
% --------------------------------------------------------------------------
%                           Through Hl(z)
% --------------------------------------------------------------------------


xout_filterL = filter(Hl,xin);


% --------------------------------------------------------------------------
%                           Through Hn(z)
% --------------------------------------------------------------------------

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
% 
% Now need the LM polyphase components ek of hn
%     
ek = myPolyphase(num_fir,1,L,M,'2');
% 
% Have to filter xin through each branch 

sumBranch = 0;

for i = (L*M):-1:1 %Starting from the LM-1 branch
%     Creating the other branches before summation
    delayedBy_a = delayseq(xout_filterL,(i-1)*a);
    downsamp = downsample(delayedBy_a,M);
    filter_polyphase = filter(ek(i,:),1,downsamp);
    upsamp = upsample(filter_polyphase,L); 

%     Sum 

    sumBranch = sumBranch + upsamp;

    if i > 1
        sumBranch = delayseq(sumBranch,-b);
    end

end

% --------------------------------------------------------------------------
%                           Through Hl(z)
% --------------------------------------------------------------------------
% 
flag = 0; %no problem, no flag needed
output_russell = filter(Hm,sumBranch);

end