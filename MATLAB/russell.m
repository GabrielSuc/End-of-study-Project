function [output_russell, flag] =  russell(z1,p1,k1,L,M,b_fir,Nl,Nm,xin,Fsin)

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
%  %Plots
nbr_samples = 10000;
subplot(8,1,1)
plot((0:1/(Fsin):(nbr_samples-1)/(Fsin)),xin(1:nbr_samples))
title(['Input Signal for Stage ', num2str(i)])
%
% --------------------------------------------------------------------------
%                           Through Hl(z)
% --------------------------------------------------------------------------


xout_filterL = filter(Hl,xin);

subplot(8,1,2)
plot((0:1/(Fsin):(nbr_samples-1)/(Fsin)),xout_filterL(1:nbr_samples))
title('Signal After First All-Pole filter')


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
ek = myPolyphase(num_fir,L,M);
% 
% Have to filter xin through each branch 

sumBranch = 0;

for i = (L*M):-1:1 %Starting from the LM-1 branch
%     Creating the other branches before summation
    delayedBy_a = delayseq(xout_filterL,(i-1)*a);
    
    subplot(8,1,3)
    plot((0:1/(Fsin):(nbr_samples-1)/(Fsin)),delayedBy_a(1:nbr_samples))
    title(['Input Signal Delayed by ', num2str(i),' times a (branch nbr ' , num2str(i), ')'])
    %
    %
    downsamp = downsample(delayedBy_a,M);
    
    subplot(8,1,4)
    plot((0:1/(Fsin/M):(nbr_samples-1)/(Fsin/M)),downsamp(1:nbr_samples))
    title('Delayed Signal Downsampled ')
    %
    %
    filter_polyphase = filter(ek(i,:),1,downsamp);
    
    subplot(8,1,5)
    plot((0:1/(Fsin/M):(nbr_samples-1)/(Fsin/M)),filter_polyphase(1:nbr_samples))
    title('Downsampled Signal filtered ')
    %
    %
    upsamp = upsample(filter_polyphase,L);
    
    subplot(8,1,6)
    plot((0:1/(Fsin*L/M):(nbr_samples-1)/(Fsin*L/M)),upsamp(1:nbr_samples))
    title('Filtered Signal Upsampled ')

%     Sum 

    sumBranch = sumBranch + upsamp;

    if i > 1
        sumBranch = delayseq(sumBranch,-b);
    end

end

subplot(8,1,7)
plot((0:1/(Fsin*L/M):(nbr_samples-1)/(Fsin*L/M)),sumBranch(1:nbr_samples))
title('Signal After Polyphase')

% --------------------------------------------------------------------------
%                           Through Hl(z)
% --------------------------------------------------------------------------
% 
flag = 0; %no problem, no flag needed
output_russell = filter(Hm,sumBranch);

subplot(8,1,8)
plot((0:1/(Fsin*L/M):(nbr_samples-1)/(Fsin*L/M)),output_russell(1:nbr_samples))
title('Output Signal (After second All-Pole Filter)')

end