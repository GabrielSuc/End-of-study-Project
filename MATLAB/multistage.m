function [output, Fint_max] = multistage(Fsin,Fp,Rp,Rs,input_signal,bestPerm,filter_choice, multistage_method)

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
                
                %Compensate gain
                ek = myPolyphase(FL(1,i)*firpm(Order,fo,ao,w),1,FL(1,i),FM(1,i),'2'); %,w
                
                

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
                    
%                     if size(sumBranch,1) > 1
%                         if size(upsamp,1) > size(sumBranch,1)
%                             sumBranch = cat(1,sumBranch, zeros(size(upsamp,1)-size(sumBranch,1),2)); 
%                         elseif size(upsamp,1) < size(sumBranch,1)
%                             upsamp = cat(1,upsamp, zeros(size(sumBranch,1)-size(upsamp,1),2));
%                         end
%                     end    
                    sumBranch = sumBranch + upsamp;

                    if k > 1
                        sumBranch = delayseq(sumBranch,-b);
                    end
                end
                

                signal = sumBranch;
                
                
                
            else
                
                           
                % Creating the filter
                Filter = dfilt.dffir(FL(1,i)*firpm(Order,fo,ao));%,w
                fvtool(Filter);
                %Saving file to use it in a C code file
                %Move it to the appropriate directory
                file = fopen(['PM_filter',num2str(i),'.txt'],'w');
                fprintf(file,'%.20e\n',vpa(Filter.Numerator));
                fclose(file);
                set(Filter,'arithmetic','fixed');
                fcfwrite(Filter,['Pm_filter', num2str(i), '.fcf']);
               
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
elseif filter_choice == 2
    
    % ---------------------------------------------------------------------
    %                          Filters' parameters
    % ---------------------------------------------------------------------
    
    %Need to adapt the ripple in the passband
    Rp = Rp/length(FL);   
        
    % Ask if we want to use polyphase
    ispolyphase = input('1 if you want to use polyphase decomposition, 0 otherwise [1]: ');
     
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
                end
                
                [signal, flag] = russell(z_ellip,p_ellip,FL(1,i)*k_ellip,FL(1,i),FM(1,i),b_ellip,Nl,Nm,signal);

                
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
            Fmax(1,i) = Fs*FL(1,1);

            
            if (multistage_method == 1)
                
                Fpassband = Fpass;
                Fcutoff = (Fmax(1,i))/2 * min(1/FL(1,1),1/FM(1,1));
            
                %By Matlab estimation
                [Order,Wp] = ellipord(Fpassband/(Fmax(1,i)/2),Fcutoff/(Fmax(1,i)/2),Rp,Rs);

                [z_ellip,p_ellip,k_ellip] = ellip(Order,Rp,Rs,Wp);
                %Creating filter
                Filter =  dfilt.df2sos(zp2sos(z_ellip,p_ellip,FL(1,1)*k_ellip));
                
                    
            else
                
                pass_bands = Fpass/(Fmax(1,i)/2);
                
                if FL(1,1) > FM(1,1)
                    stop_bands = (Fs - Fstop)/(Fmax(1,i)/2);
                else
                    stop_bands = (Fs*(FL(1,1)/FM(1,1)) - Fstop)/(Fmax(1,i)/2);
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




output = signal;

Fint_max = max(Fmax);


end