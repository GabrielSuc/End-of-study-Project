function [bestPerm,filter_choice,multistage_method] = multi_stage(L,M,Fsin,Fsout,Fp,Rp,Rs)

% This function returns the best combination of stages for different kind of filter: 
% Parks-McClellan, Elliptic, Schuessler

% Multistage decomposition used in the open-source project Smarc ParisTech
% Method based on the paper "Optimum FIR Digital Filter Implementations for Decimation, 
% Interpolation, and Narrow-Band Filtering" by Crochiere and Rabiner.

%Get the list of factors L & M for the different stages
[FM,FL] = getListStages(L,M);

%First of all, we need to find which combinaison of L and M is the best
%We create array of all the possible permutations
permsFM = perms(FM);
permsFL = perms(FL);

%High probability that there are duplicates in those previous lists
%Let's count them

if length(unique(permsFM,'rows')) == length(unique(permsFL,'rows'))
    permsFM = unique(permsFM,'rows');
    permsFL = unique(permsFL,'rows');
end    
    

if length(permsFL) == length(permsFM)

%Ask for the type of filter we want to be returned
filter_choice = input('1 to return a Parks-McClellan multistage filter, 2 for Elliptic, 3 for combination:');

    if isempty(filter_choice)
          disp('Choose a filter');
    end
    
%Decide which method to use to design the multistage filters
multistage_method = input('Choose multistage method to design the filters: [1] for classic method, [2] for Crochiere & Rabiner s method: ');
        
     if isempty(multistage_method)
          disp('Choose a multistage method');
    end   

%Global permutation
permutation = zeros(length(permsFL)*length(permsFL),length(FL) + length(FM));

%Orders of the filters accross each permutations
Order = zeros(length(permsFL)*length(permsFL),length(FL));

%Total number of MPOS
MPOS = zeros(length(permsFL)*length(permsFL),1);

%Flag if a Schuessler FIlter is created
isSchuessler = zeros(length(permsFL)*length(permsFL),length(FM) - 1);

%We now are going to perform the below algorithm over every permutations of 
%FL and Fm and find which ones gives the best result in terms of computational 
%efficiency


for l = 1:length(permsFL)
    for m = 1:length(permsFM)
        
        %Storing the different permutations tested
        permutation((l-1)*length(permsFL) + m, 1:length(FL)) = permsFL(l,:);
        permutation((l-1)*length(permsFL) + m, (length(FL) + 1):end) = permsFM(m,:);

%         Defining the differents band-edge frequencies
        Fpass = Fp;
        
        Fs = Fsin;
        
        if (multistage_method == 2)
            Fstop = Fsin/2;
        end    
%--------------------------------------------------------------------------
% ------------------------- Parks-McClellan -------------------------------
%--------------------------------------------------------------------------


        if filter_choice ==1 
        
        %Need of specific variables for the design
        a = [1 0];
        f = zeros(length(FL),length(a));

        %Ripples
        Delta1 = 10^(Rp/20); 
        Delta2 = 10^(-Rs/20);           

        %Reduce the passband ripple by factor of length(FL) so that passband
        %ripple in the cascade of the length(FL) filters doesn't exeed Delta1
        dev = [(Delta1 - 1)/(length(FL)*(Delta1 + 1)) Delta2]; 
        % -------------------------- Filters' parameters --------------------------

        %Band-edge frequencies
        if (multistage_method == 1)
            [Fpassband, Fcutoff] = deal(zeros(1,length(FM)));
        else    
            [pass_bands, stop_bands] = deal(zeros(1,length(FM)));
        end    

        %Number of MPOS at each stage
        Rmpos =  zeros(length(FM),1);
        
        wrong_stage = 0;

        for i = 1:length(FM)

            %Frequency bands
            Fmax = Fs*permsFL(l,i);
            
            if wrong_stage == 1
                break
            end    
            
            if (multistage_method == 1)
                
                if Fsin < Fsout
                    if (Fs * permsFL(l,i))/permsFM(m,i) < Fsin 
                        wrong_stage = 1;
                        continue
                    end

                elseif Fsout < Fsin

                    if (Fs * permsFL(l,i))/permsFM(m,i) < Fsout 
                        wrong_stage = 1;
                        continue
                    end

                end 
                
                
                Fpassband(i) = Fpass;
                Fcutoff(i) = (Fmax/2) * min(1/permsFL(l,i),1/permsFM(m,i));

                f(i,:) = [Fpassband(i) Fcutoff(i)];

                Order((l-1)*length(permsFL) + m,i) = firpmord(f(i,:),a,dev,Fmax);
                
                
                    
            else
                
                pass_bands(i) = Fpass/(Fmax/2);

                if permsFL(l,i) > permsFM(m,i)
                    stop_bands(i) = (Fs - Fstop)/(Fmax/2);
                else
                    stop_bands(i) = (Fs*(permsFL(l,i)/permsFM(m,i)) - Fstop)/(Fmax/2);
                end   

%                 Some of combinations can lead to negative stopband frequencies
%                 Need to skip them
                if stop_bands(i) < 0 || (stop_bands(i) < pass_bands(i))
                    continue 
                end
                
                f(i,:) = [pass_bands(i) stop_bands(i)];
                
                Order((l-1)*length(permsFL) + m,i) = firpmord(f(i,:),a,dev);
                
            end     

            %Number of MPOS
            %Defined as polyphase implementation
            cumprodM = cumprod(permsFM(m,i+1:end));
            cumprodL = cumprod(permsFL(l,i+1:end));
            if isempty(cumprodL) && isempty(cumprodM)
                cumprodL = 1;
                cumprodM = 1;
            end    
            Rmpos(i) =  ((Order((l-1)*length(permsFL) + m,i) + 1)/(permsFL(l,i))) * (cumprodL(end)/cumprodM(end)); 
            
            
            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * permsFL(l,i))/permsFM(m,i);
                                                                                 

        end
        
        
        %Computing the total number of MPOS for each permutation
        
        %We first need to remove the lines where the Order have 0
        if (multistage_method == 1)
            
            if  wrong_stage == 1 
           
            MPOS((l-1)*length(permsFL) + m) = 0;
            permutation((l-1)*length(permsFL) + m) = 0;
            
            else
            
            B = 0;
            for i = 1:length(FL)
                B = B + Rmpos(i);
            end

            MPOS((l-1)*length(permsFL) + m) = B;

            end
            
        else
            
            if find(stop_bands < 0 | (stop_bands < pass_bands))
                
                MPOS((l-1)*length(permsFL) + m) = 0;
                permutation((l-1)*length(permsFL) + m) = 0;
            
            else
            
                B = 0;
                for i = 1:length(FL)
                    B = B + Rmpos(i);
                end

                MPOS((l-1)*length(permsFL) + m) = B;

            end
            
        end
             
%--------------------------------------------------------------------------
% ----------------------------- Elliptic ----------------------------------
%--------------------------------------------------------------------------


        elseif filter_choice == 2 

        %Need to adapt the ripple in the passband
        Rp = Rp/length(FM);     

        % -------------------------- Filters' parameters --------------------------

        %Band-edge frequencies
        if (multistage_method == 1)
            [Fpassband, Fcutoff] = deal(zeros(1,length(FM)));
        else    
            [pass_bands, stop_bands] = deal(zeros(1,length(FM)));
        end 
        
        %Number of MPOS at each stage
        Rmpos =  zeros(length(FM),1);
        
        wrong_stage = 0;
      
        for i = 1:length(FM)
            
            
            %Frequency bands
            Fmax = Fs*permsFL(l,i);
            
            if wrong_stage == 1
                break
            end  
            
            if (multistage_method == 1)
                
                if Fsin < Fsout
                    if (Fs * permsFL(l,i))/permsFM(m,i) < Fsin 
                        wrong_stage = 1;
                        continue
                    end

                elseif Fsout < Fsin

                    if (Fs * permsFL(l,i))/permsFM(m,i) < Fsout 
                        wrong_stage = 1;
                        continue
                    end

                end 
                
                
                Fpassband(i) = Fpass;
                Fcutoff(i) = (Fmax/2) * min(1/permsFL(l,i),1/permsFM(m,i));


                Order((l-1)*length(permsFL) + m,i) = ellipord(Fpassband(i)/(Fmax/2),Fcutoff(i)/(Fmax/2),Rp,Rs);
                [z,p,k] = ellip(Order((l-1)*length(permsFL) + m,i),Rp,Rs,Fpassband(i)/(Fmax/2));
            
                Nzi = length(z);
                Npi = length(p);
                
                
                    
            else
                
                pass_bands(i) = Fpass/(Fmax/2);

                if permsFL(l,i) > permsFM(m,i)
                    stop_bands(i) = (Fs - Fstop)/(Fmax/2);
                else
                    stop_bands(i) = (Fs*(permsFL(l,i)/permsFM(m,i)) - Fstop)/(Fmax/2);
                end   

%                 Some of combinations can lead to negative stopband frequencies
%                 Need to skip them
                if stop_bands(i) < 0 || (stop_bands(i) < pass_bands(i))
                    continue 
                end
              
                Order((l-1)*length(permsFL) + m,i) = ellipord(pass_bands(i),stop_bands(i),Rp,Rs);
                [z,p,k] = ellip(Order((l-1)*length(permsFL) + m,i),Rp,Rs,pass_bands(i));
            
                Nzi = length(z);
                Npi = length(p);
                
            end        

            
            %Number of MPOS
            %Defined as polyphase implementation
            cumprodM = cumprod(permsFM(m,i+1:end));
            cumprodL = cumprod(permsFL(l,i+1:end));
            if isempty(cumprodL) && isempty(cumprodM)
                cumprodL = 1;
                cumprodM = 1;
            end    
            Rmpos(i) =  ((Nzi + 1)/permsFL(l,i) + ((permsFL(l,i) + permsFM(m,i) - 1)/permsFL(l,i)) * Npi)...
                * (cumprodL(end)/cumprodM(end)); 
            

            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * permsFL(l,i))/permsFM(m,i);

        end
        
        
        %Computing the total number of MADS for each permutation
        
        %We first need to remove the lines where the Order have 0
        
        if (multistage_method == 1)
            
            if  wrong_stage == 1 
           
                MPOS((l-1)*length(permsFL) + m) = 0;
                permutation((l-1)*length(permsFL) + m) = 0;
            
            else
            
                B = 0;
                for i = 1:length(FL)
                    B = B + Rmpos(i);
                end

                MPOS((l-1)*length(permsFL) + m) = B;

            end
            
        else
            
            if find(stop_bands < 0 | (stop_bands < pass_bands))
                
                MPOS((l-1)*length(permsFL) + m) = 0;
                permutation((l-1)*length(permsFL) + m) = 0;
            
            else
            
                B = 0;
                for i = 1:length(FL)
                    B = B + Rmpos(i);
                end

                MPOS((l-1)*length(permsFL) + m) = B;

            end
            
        end
        
        
        if filter_choice == 2
            Rp = Rp*length(FM); %Otherwise, it keeps dividing Rp at each for-loop
        end
        
        
        
        
        
        
%--------------------------------------------------------------------------
% ----------------------- Elliptic + Schuessler ---------------------------
%--------------------------------------------------------------------------        
        
        
        
        
        
        
        elseif filter_choice == 3 %Any combination
        
        %First filter is an IIR Elliptic filter and the rest is Schuessler
        %filters
        
        
        %Need to adapt the ripple in the passband
        Rp = Rp/length(FM); 
        
              
        %Design of the Schuessler filter
        %We want the same specs i.e. same Rp and Rs, thus we need to adapt
        %the specs of the PM filter first
        
        %Need of specific variables for the design

        a = [1 0];
        f = zeros(length(FL),length(a));

        %Ripples
        Delta1 = 10^(Rp/20); %Reduce the passband ripple by factor of length(FM) so that passband
        %ripple in the cascade of the length(FM) filters doesn't exeed Delta1
        Delta2 = 10^(-Rs/20);           


        dev = [(Delta1 - 1)/(Delta1 + 1) Delta2]; %abs(Delta1 - 1)
        % -------------------------- Filters' parameters --------------------------

        %Band-edge frequencies
        if (multistage_method == 1)
            [Fpassband, Fcutoff] = deal(zeros(1,length(FM)));
        else    
            [pass_bands, stop_bands] = deal(zeros(1,length(FM)));
        end  

        %Number of MPOS at each stage
        Rmpos =  zeros(length(FM),1);

        wrong_stage = 0;
        

%----------------------------- First Filter -------------------------------
           
        %Frequency bands
        Fmax = Fs*permsFL(l,1);  
        
        if (multistage_method == 1)

            if Fsin < Fsout
                if (Fs * permsFL(l,1))/permsFM(m,1) < Fsin 
                    wrong_stage = 1;
                    Rp = Rp*length(FM); %Otherwise, it keeps dividing Rp at each for loop
                    continue
                end

            elseif Fsout < Fsin

                if (Fs * permsFL(l,1))/permsFM(m,1) < Fsout 
                    wrong_stage = 1;
                    Rp = Rp*length(FM); %Otherwise, it keeps dividing Rp at each for loop
                    continue
                end

            end 


            Fpassband(1) = Fpass;
            Fcutoff(1) = (Fmax/2) * min(1/permsFL(l,1),1/permsFM(m,1));


            Order((l-1)*length(permsFL) + m,1) = ellipord(Fpassband(1)/(Fmax/2),Fcutoff(1)/(Fmax/2),Rp,Rs);
            [z1,p1,k1] = ellip(Order((l-1)*length(permsFL) + m,1),Rp,Rs,Fpassband(1)/(Fmax/2));

            Nz1 = length(z1);
            Np1 = length(p1);



        else

            pass_bands(1) = Fpass/(Fmax/2);

            if permsFL(l,1) > permsFM(m,1)
                stop_bands(1) = (Fs - Fstop)/(Fmax/2);
            else
                stop_bands(1) = (Fs*(permsFL(l,1)/permsFM(m,1)) - Fstop)/(Fmax/2);
            end   

%                 Some of combinations can lead to negative stopband frequencies
%                 Need to skip them
            if stop_bands(1) < 0 || (stop_bands(1) < pass_bands(1))
                Rp = Rp*length(FM); %Otherwise, it keeps dividing Rp at each for loop
                MPOS((l-1)*length(permsFL) + m) = 0;
                permutation((l-1)*length(permsFL) + m) = 0;
                continue 
            end

            Order((l-1)*length(permsFL) + m,1) = ellipord(pass_bands(1),stop_bands(1),Rp,Rs);
            [z1,p1,k1] = ellip(Order((l-1)*length(permsFL) + m,1),Rp,Rs,pass_bands(1));

            Nz1 = length(z1);
            Np1 = length(p1);

        end   

        %MPOS
        %Defined as polyphase implementation
        cumprodM = cumprod(permsFM(m,2:end));
        cumprodL = cumprod(permsFL(l,2:end));

        Rmpos(1) =  ((Nz1 + 1)/permsFL(l,1) + ((permsFL(l,1) + permsFM(m,1) - 1)/permsFL(l,1)) * Np1)...
            * (cumprodL(end)/cumprodM(end)); 

        %Need to adapt the input frequency after passing through each stage
        Fs = (Fs * permsFL(l,1))/permsFM(m,1);
            
  
%----------------------------- Other Filter -------------------------------            
    

       

        for i = 2:length(FM) %Don't forget that the first filter is an elliptic
            
            %Frequency bands
            Fmax = Fs*permsFL(l,i);
            
            if wrong_stage == 1
                break
            end  
            
            if (multistage_method == 1)
                
                if Fsin < Fsout
                    if (Fs * permsFL(l,i))/permsFM(m,i) < Fsin 
                        wrong_stage = 1;
                        continue
                    end

                elseif Fsout < Fsin

                    if (Fs * permsFL(l,i))/permsFM(m,i) < Fsout 
                        wrong_stage = 1;
                        continue
                    end

                end 
                
                
                Fpassband(i) = Fpass;
                Fcutoff(i) = (Fmax/2) * min(1/permsFL(l,i),1/permsFM(m,i));

                f(i,:) = [Fpassband(i) Fcutoff(i)];

                Order((l-1)*length(permsFL) + m,i) = firpmord(f(i,:),a,dev,Fmax);
                
                %If we want to use the spectral factorisation, the filter has to
                %have a nonnegative zerophase reponse, which implies, has to have
                %an even order

                if rem(Order((l-1)*length(permsFL) + m,i),2)~= 0
                    Order((l-1)*length(permsFL) + m,i) = Order((l-1)*length(permsFL) + m,i) + 1;
                end

                fo = [0 Fpassband(i)/(Fmax/2) Fcutoff(i)/(Fmax/2) 1];
                ao = [1 1 0 0];

                b = firpm(Order((l-1)*length(permsFL) + m,i),fo,ao); 

                %fvtool(b)

                %Now, we create the Schuessler filters

                H_schuessler = schuessler(b,Delta2);

                if H_schuessler == 0  
                    continue
                end    

                Order((l-1)*length(permsFL) + m,i) = length(H_schuessler) - 1;
                
                
                    
            else
                
                pass_bands(i) = Fpass/(Fmax/2);

                if permsFL(l,i) > permsFM(m,i)
                    stop_bands(i) = (Fs - Fstop)/(Fmax/2);
                else
                    stop_bands(i) = (Fs*(permsFL(l,i)/permsFM(m,i)) - Fstop)/(Fmax/2);
                end   

%                 Some of combinations can lead to negative stopband frequencies
%                 Need to skip them
                if stop_bands(i) < 0 || (stop_bands(i) < pass_bands(i))
                    break
                end
                
                f(i,:) = [pass_bands(i) stop_bands(i)];
                
                Order((l-1)*length(permsFL) + m,i) = firpmord(f(i,:),a,dev);
                
                %If we want to use the spectral factorisation, the filter has to
                %have a nonnegative zerophase reponse, which implies, has to have
                %an even order

                if rem(Order((l-1)*length(permsFL) + m,i),2)~= 0
                    Order((l-1)*length(permsFL) + m,i) = Order((l-1)*length(permsFL) + m,i) + 1;
                end

                fo = [0 pass_bands(i) stop_bands(i) 1];
                ao = [1 1 0 0];

                b = firpm(Order((l-1)*length(permsFL) + m,i),fo,ao); 

                %fvtool(b)

                %Now, we create the Schuessler filters

                H_schuessler = schuessler(b,Delta2);

                if H_schuessler == 0  
                    continue
                end    

                Order((l-1)*length(permsFL) + m,i) = length(H_schuessler) - 1;

            end     
             

            
             %Number of MPOS
            %Defined as polyphase implementation
            cumprodM = cumprod(permsFM(m,i+1:end));
            cumprodL = cumprod(permsFL(l,i+1:end));
            if isempty(cumprodL) && isempty(cumprodM)
                cumprodL = 1;
                cumprodM = 1;
            end    
            Rmpos(i) =  (((Order((l-1)*length(permsFL) + m,i) + 1))/(permsFL(l,i))) * (cumprodL(end)/cumprodM(end)); 
           
            Fs = (Fs * permsFL(l,i))/permsFM(m,i);
            
            %Need to create a flag to inform that a Schuessler form has
            %been used
            isSchuessler((l-1)*length(permsFL) + m,i) = 1;
            
            
        end
        
    
        
 

        %Computing the total number of MADS for each permutation
        
        %We first need to remove the lines where the Order has 0
        
        
        if (multistage_method == 1)
            
            if  wrong_stage == 1 
           
                MPOS((l-1)*length(permsFL) + m) = 0;
                permutation((l-1)*length(permsFL) + m) = 0;
            
            else
            
                B = 0;
                for i = 1:length(FL)
                    B = B + Rmpos(i);
                end

                MPOS((l-1)*length(permsFL) + m) = B;

            end
            
        else
            
            if find(stop_bands < 0 | (stop_bands < pass_bands))
                
                MPOS((l-1)*length(permsFL) + m) = 0;
                permutation((l-1)*length(permsFL) + m) = 0;
            
            else
            
                B = 0;
                for i = 1:length(FL)
                    B = B + Rmpos(i);
                end

                MPOS((l-1)*length(permsFL) + m) = B;

            end
            
        end
        

        if filter_choice == 3
            Rp = Rp*length(FM); %Otherwise, it keeps dividing Rp at each for loop
        end    
        

        end
    
    end
    
end

end

%Some of the orders can be 0 (negative stopbands), we need to remove them
Order(~all(Order,2),:) = [];
MPOS(~all(MPOS,2)) = [];
permutation(~all(permutation,2),:) = [];
isSchuessler(~any(isSchuessler,2),:) = [];

index = find(MPOS == min(MPOS)); %Can be several value (different stages leading to same efficiency)


for i = 1:length(index)
disp('-------------------------------------------------------------------')
X = ['The most efficient combination of ',num2str(length(FL)),'-stage Combination filter is of order: ',...
    num2str(sum(Order(index(i),:)))];
disp(X)
disp(' ')
if ~isempty(isSchuessler)
    if sum(isSchuessler(index(i),:)) ~= 0
        X = ['The combination consists of a first Elliptic filter, ', ...
            num2str(sum(isSchuessler(index(i),:))),' Schuessler Filters and ', ...
            num2str((length(FM)-1) - sum(isSchuessler(index(i),:))), ' Parks-McClellan filters' ];
        disp(X)
        else
        X = ['The combination consists of a first Elliptic filter and ', ...
            num2str(length(FM) - 1),' Parks-McClellan Filters' ];  
        disp(X)
    end
end
disp('')
X = ['With stages: L = ', num2str(permutation(index(i),1:length(FL))), ' and M = ', num2str(permutation(index(i),(length(FL)+1):end))];
disp(X)
disp(' ')
X = ['With a number of MPOS ', num2str(MPOS(index(i)))];
disp(X)
disp('-------------------------------------------------------------------')

bestPerm = permutation(index(i),:);
multistage_method = multistage_method;
end



end