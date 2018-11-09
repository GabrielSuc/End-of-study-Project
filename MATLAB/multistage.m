function output = multistage(Fsin,Fsout,Fp,Rp,Rs,input_signal,bestPerm,filter_choice, multistage_method)

%Get the list of factors L & M for the different stages
FL = bestPerm(1:length(bestPerm)/2);
FM = bestPerm(length(bestPerm)/2 + 1:end);

Fmax = zeros(1,length(bestPerm)/2);

%Defining the differents band-edge frequencies

Fpass = Fp;

Fs = Fsin;

if (multistage_method == 2)
    Fstop = Fsin/2;
end

% Choose how many samples to plot
nbr_samples = 10000;

% Measuring the delay introduced by each stage:
delay = zeros(1,length(bestPerm)/2 + 1);

% Filtering through the different kind of filters:

% Parks-McClellan
if filter_choice == 1
    
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
            Fmax(1,i) = Fs*FL(1,i);

            if (multistage_method == 1)
                
                Fcutoff = (Fmax(1,i)/2) * min(1/FL(1,i),1/FM(1,i));

                f = [Fpass Fcutoff];
            
                [Order,fo,ao,w] = firpmord(f,A,dev,Fmax(1,i)); 
                
                    
            else
                
                pass_bands = Fpass/(Fmax(1,i)/2);
                
                if FL(1,i) > FM(1,i)
                    stop_bands = (Fs - Fstop)/(Fmax(1,i)/2);
                else
                    stop_bands = (Fs*(FL(1,i)/FM(1,i)) - Fstop)/(Fmax(1,i)/2);
                end


                f = [pass_bands stop_bands];


                [Order,fo,ao,w] = firpmord(f,A,dev);


            end
            
  
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
                
                disp('----------------------- Delay Coefficients ------------------------')
                X = ['a = ', num2str(a), ' and b = ', num2str(b)];
                disp(X)
                disp('-------------------------------------------------------------------')
                
                % Get polyphase components and compensate gain
                ek = myPolyphase(FL(1,i)*firpm(Order,fo,ao,w),FL(1,i),FM(1,i)); %,w
                
                %Have to filter xin through each branch 
                sumBranch = 0.0;
                
                % Delay introduced by the filter
                
                % Group Delay
                gd = grpdelay(FL(1,i)*firpm(Order,fo,ao),1);
           
                % The delay for each stage is computed as follow:
    % (Delay introduced by the previous stage (in samples)/Mi + gd/LiMi)*Li
                delay(1,i+1) = (delay(1,i)/FM(1,i) + gd(1)/(FL(1,i)*FM(1,i)))*FL(1,i);
                
                                
                if i == 1 
                    signal = input_signal; %true only for the first input
                end
                
                %Plots
                subplot(7,1,1)
                plot((0:1/(Fmax(1,i)/FL(1,i)):(nbr_samples-1)/(Fmax(1,i)/FL(1,i))),signal(1:nbr_samples))
                title(['Input Signal for Stage ', num2str(i)])
                
                for k = (FL(1,i)*FM(1,i)):-1:1 %Starting from the LM-1 branch
                    %Creating the other branches before summation
                    delayedBy_a = delayseq(signal,(k-1)*a); 
                    %Plots
                    subplot(7,1,2)
                    plot((0:1/(Fmax(1,i)/FL(1,i)):(nbr_samples-1)/(Fmax(1,i)/FL(1,i))),delayedBy_a(1:nbr_samples))
                    title(['Input Signal Delayed by ', num2str(k-1),' times a (branch nbr ' , num2str(k), ')'])
                    %
                    %
                    downsamp = downsample(delayedBy_a,FM(1,i));
                    
                    subplot(7,1,3)
                    plot((0:1/(Fmax(1,i)/(FL(1,i)*FM(1,i))):(nbr_samples-1)/(Fmax(1,i)/(FL(1,i)*FM(1,i)))),downsamp(1:nbr_samples))
                    title('Delayed Signal Downsampled ')
                    %
                    %
                    filter_polyphase = filter(dfilt.dffir(ek(k,:)),downsamp);
                    
                    subplot(7,1,4)
                    plot((0:1/(Fmax(1,i)/(FL(1,i)*FM(1,i))):(nbr_samples-1)/(Fmax(1,i)/(FL(1,i)*FM(1,i)))),filter_polyphase(1:nbr_samples))
                    title('Downsampled Signal filtered ')
                    %
                    %
                    upsamp = upsample(filter_polyphase,FL(1,i)); 
                    
                    subplot(7,1,5)
                    plot((0:1/(Fmax(1,i)/FM(1,i)):(nbr_samples-1)/(Fmax(1,i)/FM(1,i))),upsamp(1:nbr_samples))
                    title('Filtered Signal Upsampled ')
                    %
                    %
                    sumBranch = sumBranch + upsamp;

                    if k > 1
                        sumBranch = delayseq(sumBranch,-b);
                    end
                    
                    subplot(7,1,6)
                    plot((0:1/(Fmax(1,i)/FM(1,i)):(nbr_samples-1)/(Fmax(1,i)/FM(1,i))),sumBranch(1:nbr_samples))
                    title('Sum of signals (different branches) with advanced samples (b)')
                    
                end
                

                signal = sumBranch;
                
                subplot(7,1,7)
                plot((0:1/(Fmax(1,i)/FM(1,i)):(nbr_samples-1)/(Fmax(1,i)/FM(1,i))),signal(1:nbr_samples))
                title('Output Signal ')
                
                
            else
                
                           
                % Creating the filter
                Filter = dfilt.dffir(FL(1,i)*firpm(Order,fo,ao));%,w
                fvtool(Filter);
                
                % Delay introduced by the filter
                % Group Delay
                gd = grpdelay(FL(1,i)*firpm(Order,fo,ao),1);
                % Since it's a FIR filter, the gd is constant
                delay(1,i+1) = gd(1)/Fmax(1,i); %Given in seconds
                
                %Saving file to use it in a C code file
                %Move it to the appropriate directory
%                 file = fopen(['PM_filter',num2str(i),'.txt'],'w');
%                 fprintf(file,'%.20e\n',vpa(Filter.Numerator));
%                 fclose(file);
%                 set(Filter,'arithmetic','fixed');
%                 fcfwrite(Filter,['Pm_filter', num2str(i), '.fcf']);
               
                if i == 1 
                    signal = input_signal; %true only for the first input
                end
                
                %Plots
                subplot(4,1,1)
                plot((0:1/(Fmax(1,i)/FL(1,i)):(nbr_samples-1)/(Fmax(1,i)/FL(1,i))),signal(1:nbr_samples))
                title(['Input Signal for Stage ', num2str(i)])
                
                %Upsampling
                signal = upsample(signal,FL(1,i));
                
                subplot(4,1,2)
                plot((0:1/Fmax(1,i):(nbr_samples-1)/Fmax(1,i)),signal(1:nbr_samples))
                title('Signal Upsampled')
                %Filtering
                input_filtered = filter(Filter,signal);
                
                subplot(4,1,3)
                plot((0:1/(Fmax(1,i)):(nbr_samples-1)/(Fmax(1,i))),input_filtered(1:nbr_samples))
                title('Signal Filtered')
                %Downsampling
                signal = downsample(input_filtered,FM(1,i));
                
                subplot(4,1,4)
                plot((0:1/(Fmax(1,i)/FM(1,i)):(nbr_samples-1)/(Fmax(1,i)/FM(1,i))),signal(1:nbr_samples))
                title('Signal Downsampled')
                
            end    
            
         
            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * FL(1,i))/FM(1,i);
            
            
    end         

    %signal = signal;%/max(signal); %_To prevent data from clipping when writing file
    
    
% Elliptic    
elseif filter_choice == 2
    
    % ---------------------------------------------------------------------
    %                          Filters' parameters
    % ---------------------------------------------------------------------
    
    %Need to adapt the ripple in the passband
    Rp = Rp/length(FL);   
        
    % Ask if we want to use polyphase
    ispolyphase = input('1 if you want to use polyphase decomposition, 0 otherwise [1]: ');
    
    % Group delay for non-polyphase implementation
%     gp_direct = zeros(512,length(FL));
%     w_direct = zeros(512,length(FL));
%     f_direct = zeros(512,length(FL));
     
    for i = 1:length(FL)

%           % Frequency bands
            % Frequency bands
            Fmax(1,i) = Fs*FL(1,i);

            if (multistage_method == 1)
                
                Fcutoff = (Fmax(1,i)/2) * min(1/FL(1,i),1/FM(1,i));

            
                % Getting the order of the filters
                [Order,Wp] = ellipord(Fpass/(Fmax(1,i)/2),Fcutoff/(Fmax(1,i)/2),Rp,Rs);
                [z_ellip,p_ellip,k_ellip] = ellip(Order,Rp,Rs,Wp);
                [b_ellip,a_ellip] = ellip(Order,Rp,Rs,Wp);
                
                    
            else
                
                pass_bands = Fpass/(Fmax(1,i)/2);
                
                if FL(1,i) > FM(1,i)
                    stop_bands = (Fs - Fstop)/(Fmax(1,i)/2);
                else
                    stop_bands = (Fs*(FL(1,i)/FM(1,i)) - Fstop)/(Fmax(1,i)/2);
                end



                % Getting the order of the filters
                [Order,Wp] = ellipord(pass_bands,stop_bands,Rp,Rs);
                [z_ellip,p_ellip,k_ellip] = ellip(Order,Rp,Rs,Wp);
                [b_ellip,a_ellip] = ellip(Order,Rp,Rs,Wp);

                
                
            end
            
            
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

%                     elseif rem(Nl,2) ~= 0
%                         error('Nl has to be even')
                    end

                else

                    disp(['Np = ', num2str(Np)])    
                        
                    Nl = input('Choose even number of pole Nl < Np to be assigned to coefficient L: ');
                    
                    if isempty(Nl)
                          error('Choose Nl');

                    elseif Nl > Np
                        error('Nl has to be smaller than Np')

%                     elseif rem(Nl,2) ~= 0
%                         error('Nl has to be even')
                    end

                end

                Nm = Np - Nl;
    
                if i == 1 
                    signal = input_signal; %true only for the first input
                    
                    % Need to set the group delay and related frequencies
                    % to 0 before the first stage
                    gd_russell = [];
                    f_russell = [];
                    
                end
                
                [signal, flag, gd_russell, f_russell] = russell(z_ellip,p_ellip,FL(1,i)*k_ellip,...
                    FL(1,i),FM(1,i),b_ellip,Nl,Nm,signal,Fmax(1,i),gd_russell, f_russell);

                
                %In case the above technique is not feasible 
                if flag == 1
                    
                    disp('Could not create Russell filter')
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
                fvtool(Filter);
                
                % Group delay for each stage
                [gd_direct(:,i), w_direct(:,i)] = grpdelay(Filter); 
                % Have to relate the stage frequencies to the output one
                cumprodM = cumprod(FM(1,i:end));
                cumprodL = cumprod(FL(1,i+1:end));
                if isempty(cumprodL)
                    cumprodL = 1;
                elseif isempty(cumprodM)    
                    cumprodM = 1;
                end 
                f_direct(:,i) = w_direct(:,i).*((Fmax(1,i)/(2*pi))*cumprodL(end)/cumprodM(end));
                
                % Have to compute the time-scaled delay with the right
                % sampling frequency;
                gd_direct(:,i) = gd_direct(:,i)./Fmax(1,i);
                
                if i == 1 
                    signal = input_signal; %true only for the first input
                end
                
                %Plots
                subplot(4,1,1)
                plot((0:1/(Fmax(1,i)/FL(1,i)):(nbr_samples-1)/(Fmax(1,i)/FL(1,i))),signal(1:nbr_samples))
                title(['Input Signal for Stage ', num2str(i)])
                
                %fvtool(Filter)
                
                %Upsampling
                signal = upsample(signal,FL(1,i));
                
                subplot(4,1,2)
                plot((0:1/(Fmax(1,i)):(nbr_samples-1)/(Fmax(1,i))),signal(1:nbr_samples))
                title('Signal Upsampled')
                %Filtering
                input_filtered = filter(Filter,signal);
                
                subplot(4,1,3)
                plot((0:1/(Fmax(1,i)):(nbr_samples-1)/(Fmax(1,i))),input_filtered(1:nbr_samples))
                title('Signal Filtered')
                %Downsampling
                signal = downsample(input_filtered,FM(1,i));
                
                subplot(4,1,4)
                plot((0:1/(Fmax(1,i)/FM(1,i)):(nbr_samples-1)/(Fmax(1,i)/FM(1,i))),signal(1:nbr_samples))
                title('Signal Downsampled')
%                 output = signal;
%                 return
                
            end    
            
         
            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * FL(1,i))/FM(1,i);
            
            
    end
    
    

else
        
        %Combination: First filter is an IIR Elliptic filter and the rest is Schuessler
        %filters
        
        %Need to adapt the ripple in the passband
        Rp = Rp/length(FM); 
        
              
        %Design of the Schuessler filter
        %We want the same specs i.e. same Rp and Rs, thus we need to adapt
        %the specs of the PM filter first
        
        %Need of specific variables for the design

        a = [1 0];

        %Ripples
        Delta1 = 10^(Rp/20); %Reduce the passband ripple by factor of length(FM) so that passband
        %ripple in the cascade of the length(FM) filters doesn't exeed Delta1
        Delta2 = 10^(-Rs/20);           


        dev = [(Delta1 - 1)/(Delta1 + 1) Delta2]; %abs(Delta1 - 1)

%----------------------------- First Filter -------------------------------
           
            %Frequency bands
            Fmax(1,1) = Fs*FL(1,1);

            
            if (multistage_method == 1)
                
                Fpassband = Fpass;
                Fcutoff = (Fmax(1,1))/2 * min(1/FL(1,1),1/FM(1,1));
            
                %By Matlab estimation
                [Order,Wp] = ellipord(Fpassband/(Fmax(1,1)/2),Fcutoff/(Fmax(1,1)/2),Rp,Rs);

                [z_ellip,p_ellip,k_ellip] = ellip(Order,Rp,Rs,Wp);
                %Creating filter
                Filter =  dfilt.df2sos(zp2sos(z_ellip,p_ellip,FL(1,1)*k_ellip));
                
                    
            else
                
                pass_bands = Fpass/(Fmax(1,1)/2);
                
                if FL(1,1) > FM(1,1)
                    stop_bands = (Fs - Fstop)/(Fmax(1,1)/2);
                else
                    stop_bands = (Fs*(FL(1,1)/FM(1,1)) - Fstop)/(Fmax(1,1)/2);
                end



                %By Matlab estimation
                [Order,Wp] = ellipord(pass_bands,stop_bands,Rp,Rs);

                [z_ellip,p_ellip,k_ellip] = ellip(Order,Rp,Rs,Wp);
                %Creating filter
                Filter =  dfilt.df2sos(zp2sos(z_ellip,p_ellip,FL(1,1)*k_ellip));

                
                
            end
            
            
            
            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * FL(1,1))/FM(1,1);
            
            %Filtering through first filter
            signal = input_signal;
            
            %Upsampling
            signal = upsample(signal,FL(1,1));
            %Filtering
            input_filtered = filter(Filter,signal);
            %Downsampling
            signal = downsample(input_filtered,FM(1,1));
            
            
%----------------------------- Other Filter -------------------------------            
    

       

        for i = 2:length(FM) %Don't forget that the first filter is an elliptic
            
            %Frequency bands
            Fmax(1,i) = Fs*FL(1,i);
            
            if (multistage_method == 1)
                
                %Defining the limit frequencies for the design
                Fcutoff = (Fmax(1,i)/2) * min(1/FL(1,i),1/FM(1,i));

                f = [Fpass Fcutoff];
            
                [Order,fo,ao,w] = firpmord(f,a,dev,Fmax(1,i)); 
                
                    
            else
                
                pass_bands = Fpass/(Fmax(1,i)/2);
                
                if FL(1,i) > FM(1,i)
                    stop_bands = (Fs - Fstop)/(Fmax(1,i)/2);
                else
                    stop_bands = (Fs*(FL(1,i)/FM(1,i)) - Fstop)/(Fmax(1,i)/2);
                end


                f = [pass_bands stop_bands];

                [Order,fo,ao,w] = firpmord(f,a,dev);

            end
            
            
            
            
            %If we want to use the spectral factorisation, the filter has to
            %have a nonnegative zerophase reponse, which implies, has to have
            %an even order
        
            if rem(Order,2)~= 0
                Order = Order + 1;
            end
            
            b = firpm(Order,fo,ao,w); 
            
            %fvtool(b)
            
            %Now, we create the Schuessler filters
            
            H_schuessler = schuessler(b,Delta2);
            
            if H_schuessler == 0  
                %Upsampling
                signal = upsample(signal,FL(1,i));
                %Filtering
                input_filtered = filter(b,1,signal);
                %Downsampling
                signal = downsample(input_filtered,FM(1,i));
                disp('Cannot create Schuessler filter');
                continue
            end    
            
            
            %Upsampling
            signal = upsample(signal,FL(1,i));
            %Filtering
            input_filtered = filter(H_schuessler,1,signal);
            %Downsampling
            signal = downsample(input_filtered,FM(1,i));
            


        end


end


disp('-------------------------- Overall Delay --------------------------')

if filter_choice == 1 && ispolyphase == 0
    overall_delay = sum(delay);
    
    disp('Overall delay between input and output signal, for both channel, is: ');
    disp([num2str(overall_delay), ' s']);
    disp('-------------------------------------------------------------------')
elseif filter_choice == 1 && ispolyphase == 1
    overall_delay = delay(1,end)/Fsout;
    
    disp('Overall delay between input and output signal, for both channel, is: ');
    disp([num2str(overall_delay), ' s']);
    disp('-------------------------------------------------------------------')
elseif filter_choice == 2 && ispolyphase == 0
    overall_delay = 0;
    
    disp(['Overall delay between input and output signal, for both channel, is'...
        ' displayed by the group delay plot']);
    disp('-------------------------------------------------------------------')
    
    % Ploting the overall delay
    figure
    gd_overall = sum(gd_direct,2);
    f_overall = f_direct(:,1); % all the columns are the same
    plot(f_overall,gd_overall);
    title('Overall Delay using Elliptic filters with Polyphase (Russell)')
    ylabel('Delay (in s)')
    xlabel('Frequency (in Hz)')
    grid on;
    
elseif filter_choice == 2 && ispolyphase == 1
    overall_delay = 0;
    
    disp(['Overall delay between input and output signal, for both channel, is'...
        ' displayed by the group delay plot']);
    disp('-------------------------------------------------------------------')
    
    % Ploting the overall delay
    figure
    plot(f_russell,gd_russell);
    title('Overall Delay using Elliptic filters without Polyphase')
    ylabel('Delay (in s)')
    xlabel('Frequency (in Hz)')
    grid on;
end



disp('------------------------ Hardware Complexity -----------------------')
disp(['The maximum intermediary frequency is: ', num2str(max(Fmax)), ' Hz']);
disp('-------------------------------------------------------------------')


output = signal;

end