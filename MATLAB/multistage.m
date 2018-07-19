function output = multistage(L,M,Fsin,Fsout,Fp,Rp,Rs,input_signal,bestPerm,manual)

%Get the list of factors L & M for the different stages
FL = bestPerm(1:length(bestPerm)/2);
FM = bestPerm(length(bestPerm)/2 + 1:end);

    
%Have to inverse Fx and Fy considering the overall system as a decimator or an interpolator
if M > L
     [Fsout, Fsin] = deal(Fsin,Fsout);
end

%Defining the differents band-edge frequencies
if Fsin < Fsout
    Fstop = Fsin/2;
else 
    Fstop = Fsout/2;
end

Fpass = Fp;

Fs = Fsin;

% Filtering through the different kind of filters:
% Parks-McClellan
if manual == 1
    
    % ---------------------------------------------------------------------
    %                          Filters' parameters
    % ---------------------------------------------------------------------
    
    % Need of specific variables for the design
    A = [1 0];

    % Ripples
    Delta1 = 10^(Rp/20); 
    Delta2 = 10^(-Rs/20);           

    % Deviation
    % Reduce the passband ripple by factor of length(FL) so that passband
    % ripple in the cascade of the length(FL) filters doesn't exeed Delta1
    dev = [(Delta1 - 1)/(length(FL)*(Delta1 + 1)) Delta2]; 
    
    % Ask if we want to use polyphase
    ispolyphase = input('1 if you want to use polyphase decomposition, 0 otherwise [1]: ');
     
    for i = 1:length(FL)

            % Frequency bands
            Fmax = Fs*FL(1,i);

            Fcutoff = (Fmax)/2 * min(1/FL(1,i),1/FM(1,i));
            
            
            f = [Fpass Fcutoff];

            [Order,fo,ao,w] = firpmord(f,A,dev,Fmax);

            
            % Polyphase decomposition    
            if isempty(ispolyphase)
                 ispolyphase = 1;
            end
            
            if ispolyphase
                
                %We need to have the coefficient a and b of the delays
                a = 0;
                b = 0.1;

                while(rem(b,1)~=0)
                    a = a + 1;
                    b = (FL(1,i)*a - 1)/FM(1,i);
                end   
                
                %Compensate gain
                ek = myPolyphase(FL(1,i)*firpm(Order,fo,ao,w),1,FL(1,i),FM(1,i),'2');

                %Have to filter xin through each branch 
                sumBranch = 0;
                
                if i == 1 
                    signal = input_signal; %true only for the first input
                end
                
                for k = (FL(1,i)*FM(1,i)):-1:1 %Starting from the LM-1 branch
                    %Creating the other branches before summation
                    delayedBy_a = delayseq(signal,(k-1)*a); 
                    downsamp = downsample(delayedBy_a,FM(1,i));
                    %Implement polyphase components as folded structures
                    filter_polyphase = filter(ek(k,:),1,downsamp);
                    upsamp = upsample(filter_polyphase,FL(1,i)); 

                    %Sum 

                    sumBranch = sumBranch + upsamp;

                    if k > 1
                        sumBranch = delayseq(sumBranch,-b);
                    end
                end
                

                signal = sumBranch;
                
            else
                
                % Creating the filter
                Filter = dfilt.dfsymfir(FL(1,i)*firpm(Order,fo,ao,w));
               
                if i == 1 
                    signal = input_signal; %true only for the first input
                end
                
                %Upsampling
                signal = upsample(signal,FL(1,i));
                %Filtering
                input_filtered = filter(Filter,signal);
                %Downsampling
                signal = downsample(input_filtered,FM(1,i));
                
            end    
            
         
            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * FL(1,i))/FM(1,i);
            
            
    end         

    %signal = signal;%/max(signal); %_To prevent data from clipping when writing file
    
    
% Elliptic    
elseif manual == 2
    
    % ---------------------------------------------------------------------
    %                          Filters' parameters
    % ---------------------------------------------------------------------
    
    %Need to adapt the ripple in the passband
    Rp = Rp/length(FL);   
        
    % Ask if we want to use polyphase
    ispolyphase = input('1 if you want to use polyphase decomposition, 0 otherwise [1]: ');
     
    for i = 1:length(FL)

%             % Frequency bands
                        % Frequency bands
            Fmax = Fs*FL(1,i);

            Fcutoff = (Fmax)/2 * min(1/FL(1,i),1/FM(1,i));
            

            % Getting the order of the filters
            [Order,Wp] = ellipord(Fpass/(Fmax/2),Fcutoff/(Fmax/2),Rp,Rs);
            [z_ellip,p_ellip,k_ellip] = ellip(Order,Rp,Rs,Wp);
            
            % Polyphase decomposition    
            if isempty(ispolyphase)
                 ispolyphase = 1;
            end
            
            if ispolyphase
                
                %Russell

                Np = length(p_ellip); %Number of poles

                %Have to choose Nl and Nm -> put a munual stuff directly inside function 
                %The gain in efficiency is the same regarless of how many poles of H(z)
                %areassigned to Hl(z) ans how many bto Hm(z). However, complex-conjugate
                %pole pairs should not be separated since this would make the filtr
                %coefficients complex. Therefore, choose Nl even and Nm odd if length
                %filter odd.

                if rem(Np,2) ~=0
                
                disp(['Np = ', num2str(Np)])    
                        
                Nl = input('Choose even number of pole Nl < Np to be assigned to coefficient L: ');

                    if isempty(Nl)
                          error('Choose Nl');

                    elseif Nl > Np
                        error('Nl has to be smaller than Np')

                    elseif rem(Nl,2) ~= 0
                        error('Nl has to be even')
                    end

                else

                    disp(['Np = ', num2str(Np)])    
                        
                    Nl = input('Choose even number of pole Nl < Np to be assigned to coefficient L: ');
                    
                    if isempty(Nl)
                          error('Choose Nl');

                    elseif Nl > Np
                        error('Nl has to be smaller than Np')

                    elseif rem(Nl,2) ~= 0
                        error('Nl has to be even')
                    end

                end

                Nm = Np - Nl;
    
                if i == 1 
                    signal = input_signal; %true only for the first input
                end
                
                [signal, flag] = russell(z_ellip,p_ellip,FL(1,i)*k_ellip,FL(1,i),FM(1,i),Nl,Nm,signal);

                
                %In case the above technique is not feasible 
                if flag == 1
                    
                    %Creating the Elliptic Filters
                    Filter = dfilt.df2sos(zp2sos(z_ellip,p_ellip,FL(1,i)*k_ellip));

                    %Upsampling
                    signal = upsample(signal,FL(1,i));
                    %Filtering
                    input_filtered = filter(Filter,signal);
                    %Downsampling
                    signal = downsample(input_filtered,FM(1,i));
                    
                end
                
            else
                
                %Creating the Elliptic Filters
                Filter = dfilt.df2sos(zp2sos(z_ellip,p_ellip,FL(1,i)*k_ellip));
               
                if i == 1 
                    signal = input_signal; %true only for the first input
                end
                
                fvtool(Filter)
                
                %Upsampling
                signal = upsample(signal,FL(1,i));
                %Filtering
                input_filtered = filter(Filter,signal);
                %Downsampling
                signal = downsample(input_filtered,FM(1,i));
                
%                 output = signal;
%                 return
                
            end    
            
         
            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * FL(1,i))/FM(1,i);
            
            
    end
    
% Schuessler    
else








end




output = signal;




end