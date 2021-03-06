function [output_russell, flag, delay_current_stage] =  russell(z1,p1,k1,L,M,b_fir,Nl,Nm,xin,Fmax, delay_previous_stage,nbr_samples)

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

% Group delay for this stage
[gd_ap1,w_ap1] = grpdelay(Hl);

% Have to compute the time-scaled delay with the right
% sampling frequency;
gd_ap1 = gd_ap1./(Fmax/L); %divide by Fsin

% We can convert the digital frequencies to the analog ones
% Have to relate them to the frequency at the output of the stage
f_ap1 = w_ap1.*((Fmax/L)/(2*pi)); 

% Making a single array
delay_ap1 = [f_ap1,gd_ap1];

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

% Group delay for this stage
[gd_ap2,w_ap2] = grpdelay(Hm);

% Have to compute the time-scaled delay with the right
% sampling frequency;
gd_ap2 = gd_ap2./(Fmax/M); %divide by F_intermediary1

% We can convert the digital frequencies to the analog ones
% Have to relate them to the frequency at the output of the stage
f_ap2 = w_ap2.*((Fmax/M)/(2*pi)); 

% Making a single array
delay_ap2 = [f_ap2,gd_ap2];
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

% Group delay for this stage
[gd_fir,w_fir] = grpdelay(Hn);

% Have to compute the time-scaled delay with the right
% sampling frequency;
gd_fir = gd_fir./(Fmax/(L*M)); %divide by Fsin/M

% We can convert the digital frequencies to the analog ones
% Have to relate them to the frequency at the output of the stage
f_fir = w_fir.*((Fmax/(L*M))/(2*pi)); 

% Making a single array
delay_fir = [f_fir,gd_fir];

% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
%                       Filtering the signal
% --------------------------------------------------------------------------
% --------------------------------------------------------------------------
% 
%  %Plots
% subplot(9,1,1)
% plot((0:1/(Fmax/L):(nbr_samples-1)/(Fmax/L)),xin(1:nbr_samples))
% title(['Input Signal for Stage ', num2str(i)])
%
% --------------------------------------------------------------------------
%                           Through Hl(z)
% --------------------------------------------------------------------------


xout_filterL = filter(Hl,xin(2:6));

% subplot(9,1,2)
% plot((0:1/(Fmax/L):(nbr_samples-1)/(Fmax/L)),xout_filterL(1:nbr_samples))
% title('Signal After First All-Pole filter')


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
    
%     subplot(9,1,3)
%     plot((0:1/(Fmax/L):(nbr_samples-1)/(Fmax/L)),delayedBy_a(1:nbr_samples))
%     title(['Input Signal Delayed by ', num2str(i),' times a (branch nbr ' , num2str(i), ')'])
    %
    %
    downsamp = downsample(delayedBy_a,M);
    
%     subplot(9,1,4)
%     plot((0:1/((Fmax/L)/M):(nbr_samples-1)/((Fmax/L)/M)),downsamp(1:nbr_samples))
%     title('Delayed Signal Downsampled ')
    %
    %
    filter_polyphase = filter(ek(i,:),1,downsamp);
    
%     subplot(9,1,5)
%     plot((0:1/((Fmax/L)/M):(nbr_samples-1)/((Fmax/L)/M)),filter_polyphase(1:nbr_samples))
%     title('Downsampled Signal filtered ')
    %
    %
    upsamp = upsample(filter_polyphase,L);
    
%     subplot(9,1,6)
%     plot((0:1/(Fmax/M):(nbr_samples-1)/(Fmax/M)),upsamp(1:nbr_samples))
%     title('Filtered Signal Upsampled ')

%     Sum 

    sumBranch = sumBranch + upsamp;

    if i > 1
        sumBranch = delayseq(sumBranch,-b);
    end
    
%     subplot(9,1,7)
%     plot((0:1/(Fmax/M):(nbr_samples-1)/(Fmax/M)),sumBranch(1:nbr_samples))
%     title('Sum of signals (different branches) with advanced samples (b)')

end

% subplot(9,1,8)
% plot((0:1/(Fmax/M):(nbr_samples-1)/(Fmax/M)),sumBranch(1:nbr_samples))
% title('Signal After Polyphase')

% --------------------------------------------------------------------------
%                           Through Hl(z)
% --------------------------------------------------------------------------
% 
flag = 0; %no problem, no flag needed
output_russell = filter(Hm,sumBranch);

% subplot(9,1,9)
% plot((0:1/(Fmax/M):(nbr_samples-1)/(Fmax/M)),output_russell(1:nbr_samples))
% title('Output Signal (After second All-Pole Filter)')

% --------------------------------------------------------------------------
%                           Overall Delay
% --------------------------------------------------------------------------


% Adapting the frequencies of each filter to match each other
delay_ap1(:,1) = delay_ap1(:,1).*(L/M);
delay_fir(:,1) = delay_fir(:,1).*(L);
%Nothing to do for the second all-pole filter

% Taking into account the delays from the previous stages and adapting the
% frequencies to the current stage
delay_previous_stage(:,1) = delay_previous_stage(:,1).*(L/M);

% Combining the previous stage delays with the current ones
delay_3_filters = linkArray(linkArray(delay_ap1,delay_ap2),delay_fir);
delay_current_stage = linkArray(delay_3_filters, delay_previous_stage);


end