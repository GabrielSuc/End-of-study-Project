function [bestPerm,manual] = multi_stage(L,M,Fsin,Fsout,Fp,Rp,Rs)

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
manual = input('1 to return a Parks-McClellan multistage filter, 2 for Elliptic, 3 for combination:');

    if isempty(manual)
          disp('Choose a filter');
    end
    


%Global permutation
permutation = zeros(length(permsFL)*length(permsFL),length(FL) + length(FM));

%Orders of the filters accross each permutations
Order = zeros(length(permsFL)*length(permsFL),length(FL));

%Total number of MADS accross each permutations
MADS = zeros(length(permsFL)*length(permsFL),1);

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
        
%--------------------------------------------------------------------------
% ------------------------- Parks-McClellan -------------------------------
%--------------------------------------------------------------------------


        if manual ==1 


        

        %Need of specific variables for the design
        a = [1 0];
        f = zeros(length(FL),length(a));

        %Ripples
        Delta1 = 10^(Rp/20); %Reduce the passband ripple by factor of length(FM) so that passband
        %ripple in the cascade of the length(FM) filters doesn't exeed Delta1
        Delta2 = 10^(-Rs/20);           


        dev = [(Delta1 - 1)/(length(FL)*(Delta1 + 1)) Delta2]; %abs(Delta1 - 1)
        % -------------------------- Filters' parameters --------------------------

        %Band-edge frequencies
        [Fpassband, Fcutoff] = deal(zeros(1,length(FM)));
        
        %Number of MPOS at each stage
        Rmpos =  zeros(length(FM),1);

        %Final Filter
        H = cell(1,length(FM));
        
        wrong_stage = 0;

        for i = 1:length(FM)

%             %Frequency bands
            Fmax = Fs*permsFL(l,i);
%             pass_bands(i) = Fpass/Fmax;
%             if permsFL(l,i) > permsFM(m,i)
%                 stop_bands(i) = (Fs - Fstop)/Fmax;
%             else
%                 stop_bands(i) = (Fs*(permsFL(l,i)/permsFM(m,i)) - Fstop)/Fmax;
%             end   
% 
%             %Some of combinations can lead to negative stopband frequencies
%             %Need to pass them
%             if stop_bands(i) < 0 || (stop_bands(i) < pass_bands(i))
%                 continue 
%             end    
%             
            
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
%             
%             
            f(i,:) = [Fpassband(i) Fcutoff(i)];

%             f(i,:) = [pass_bands(i) stop_bands(i)];

            [Order((l-1)*length(permsFL) + m,i),fo,ao,w] = firpmord(f(i,:),a,dev,Fmax);
            
            %Defining the limit frequencies for the design
%             f(i,:) = [pass_bands(i) stop_bands(i)];
            %fo(i,:) = [0 pass_bands(i) stop_bands(i) 1];

            %Getting the order of the filters
            %By Remez estimation
%             Order((l-1)*length(permsFL) + m,i) = ceil(remlpord(pass_bands(i),stop_bands(i),dev(1),dev(2)));
            %By Matlab estimation
%             [Order((l-1)*length(permsFL) + m,i),fo,ao,w] = firpmord(f(i,:),a,dev);


            %Filter length must be must be 2*K*Mi +1 so that the delay is
            %integer 
%             k = ceil((Order((l-1)*length(permsFL) + m,i) - 1)/(2*permsFM(m,i)));
%             while (2*k*permsFM(m,i) + 1 < Order((l-1)*length(permsFL) + m,i))
%                 k = k + 1;
%             end
%             Order((l-1)*length(permsFL) + m,i) = 2*k*permsFM(m,i);
            
            %Number of multiplies ans adds(MADS) which must be performed by sec
            %Rmads(i) = ((Order((l-1)*length(permsFL) + m,i)+1)*Fs)/permsFM(m,i);
            
            
            %Number of MPOS
            %Defined as polyphase implementation
            cumprodM = cumprod(permsFM(m,i+1:end));
            cumprodL = cumprod(permsFL(l,i+1:end));
            if isempty(cumprodL) && isempty(cumprodM)
                cumprodL = 1;
                cumprodM = 1;
            end    
            Rmpos(i) =  (((Order((l-1)*length(permsFL) + m,i) + 1)*permsFM(m,i))/permsFL(l,i)) * (cumprodL(end)/cumprodM(end)); 
            
           
           
            
            %Creating the P-M filters
            %H{i} = dfilt.dfsymfir(firpm(Order(k,i),fo,ao,w)); %Folded implementation, reduces multiples by 2
            %fvtool(H{i})
            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * permsFL(l,i))/permsFM(m,i);
                                                                                 

        end
        
        
        %Computing the total number of MADS for each permutation
        
        %We first need to remove the lines where the Order have 0
        
        if wrong_stage == 1 %find(stop_bands < 0 | (stop_bands < pass_bands)) %
            %MADS((l-1)*length(permsFL) + m) = 0;
            MPOS((l-1)*length(permsFL) + m) = 0;
            permutation((l-1)*length(permsFL) + m) = 0;
            
        else
            
        [A,B] = deal(0);
        for i = 1:length(FL)
            %A = A + Rmads(i); 
            B = B + Rmpos(i);
        end

        %MADS((l-1)*length(permsFL) + m) = A;
        MPOS((l-1)*length(permsFL) + m) = B;
        
        end
        
             
%--------------------------------------------------------------------------
% ----------------------------- Elliptic ----------------------------------
%--------------------------------------------------------------------------


        elseif manual == 2 
            
        %coef = cell(1,3);

        %Need to adapt the ripple in the passband
        Rp = Rp/length(FM);     

        % -------------------------- Filters' parameters --------------------------

        %Band-edge frequencies
        [Fpassband,Fcutoff] = deal(zeros(1,length(FM))); %Or FL, same length 
        
        %Number of MPOS at each stage
        Rmpos =  zeros(length(FM),1);
        
        wrong_stage = 0;
      
        for i = 1:length(FM)
            
            
            %Frequency bands
            Fmax = Fs*permsFL(l,i);
            
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
            Fcutoff(i) = (Fmax)/2 * min(1/permsFL(l,i),1/permsFM(m,i));
            
            
            f(i,:) = [Fpassband(i) Fcutoff(i)];

            Order((l-1)*length(permsFL) + m,i) = ellipord(Fpassband(i)/(Fmax/2),Fcutoff(i)/(Fmax/2),Rp,Rs);
            
            %Number of MPOS
            %Defined as polyphase implementation
            cumprodM = cumprod(permsFM(m,i+1:end));
            cumprodL = cumprod(permsFL(l,i+1:end));
            if isempty(cumprodL) && isempty(cumprodM)
                cumprodL = 1;
                cumprodM = 1;
            end    
            Rmpos(i) =  (((Order((l-1)*length(permsFL) + m,i) + 1)*permsFM(m,i))/permsFL(l,i)) * (cumprodL(end)/cumprodM(end)); 
            

            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * permsFL(l,i))/permsFM(m,i);

        end
        
        
        %Computing the total number of MADS for each permutation
        
        %We first need to remove the lines where the Order have 0
        
        if wrong_stage == 1
            %MADS((l-1)*length(permsFL) + m) = 0;
            MPOS((l-1)*length(permsFL) + m) = 0;
            permutation((l-1)*length(permsFL) + m) = 0;
            
        else
            
        [A,B] = deal(0);
        for i = 1:length(FL)
            %A = A + Rmads(i); 
            B = B + Rmpos(i);
        end

        %MADS((l-1)*length(permsFL) + m) = A;
        MPOS((l-1)*length(permsFL) + m) = B;
        
        end
        
        
        if manual == 2
            Rp = Rp*length(FM); %Otherwise, it keeps dividing Rp at each for loop
        end
        
        
        
        
        
        
%--------------------------------------------------------------------------
% ----------------------- Elliptic + Schuessler ---------------------------
%--------------------------------------------------------------------------        
        
        
        
        
        
        
        elseif manual == 3 %Any combination
        
        %First filter is an IIR Elliptic filter and the rest is Schuessler
        %filters
        
        %Elliptic Filter
        %coef = cell(1,3);
        
        %Need to adapt the ripple in the passband
        Rp = Rp/length(FM); %Faut convertir en lineaire
        
              
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
        [Fpassband,Fcutoff] = deal(zeros(1,length(FM))); %Or FL, same length 

        %Rmads = zeros(length(FM),1);
        
        %Number of MPOS at each stage
        Rmpos =  zeros(length(FM),1);

        wrong_stage = 0;
        
        
%         [coef{1,1},coef{1,2},coef{1,3}] = ellip(Order(1),Rp,Rs,Wp); %[z,p,k]
%         H{1} = dfilt.df2sos(zp2sos(coef{1,1},coef{1,2},coef{1,3})); %Avoid error of quantization
%         fvtool(H{1})
        


%----------------------------- First Filter -------------------------------
           
            %Frequency bands
            Fmax = Fs*permsFL(l,1);

            
            %Frequency bands
            if Fsin < Fsout

                if (Fs * permsFL(l,1))/permsFM(m,1) < Fsin 
                    wrong_stage = 1;
                    Rp = Rp*length(FM);
                    MPOS((l-1)*length(permsFL) + m) = 0;
                    permutation((l-1)*length(permsFL) + m) = 0;
                    continue
                end
            
            elseif Fsout < Fsin
                
                if (Fs * permsFL(l,1))/permsFM(m,1) < Fsout 
                    wrong_stage = 1;
                    Rp = Rp*length(FM);
                    MPOS((l-1)*length(permsFL) + m) = 0;
                    permutation((l-1)*length(permsFL) + m) = 0;
                    continue
                end
             
            end   

            Fpassband(1) = Fpass;
            Fcutoff(1) = (Fmax)/2 * min(1/permsFL(l,1),1/permsFM(m,1));
            
            %By Matlab estimation
            Order((l-1)*length(permsFL) + m,1) = ellipord(Fpassband(1)/(Fmax/2),Fcutoff(1)/(Fmax/2),Rp,Rs);




            
            %MPOS
            %Defined as polyphase implementation
            cumprodM = cumprod(permsFM(m,2:end));
            cumprodL = cumprod(permsFL(l,2:end));
                
            Rmpos(1) =  (((Order((l-1)*length(permsFL) + m,1) + 1)*permsFM(m,1))/permsFL(l,1)) * (cumprodL(end)/cumprodM(end));

            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * permsFL(l,1))/permsFM(m,1);
            
  
%----------------------------- Other Filter -------------------------------            
    

       

        for i = 2:length(FM) %Don't forget that the first filter is an elliptic
            
            %Frequency bands
            Fmax = Fs*permsFL(l,i);
            
            %Frequency bands
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

            %Defining the limit frequencies for the design
            Fpassband(i) = Fpass;
            Fcutoff(i) = (Fmax)/2 * min(1/permsFL(l,i),1/permsFM(m,i));
            
            f(i,:) = [Fpassband(i) Fcutoff(i)];
            

            
            [Order((l-1)*length(permsFL) + m,i),fo,ao,w] = firpmord(f(i,:),a,dev,Fmax);
            

            %Number of MPOS
            %Defined as polyphase implementation
            cumprodM = cumprod(permsFM(m,i+1:end));
            cumprodL = cumprod(permsFL(l,i+1:end));
            if isempty(cumprodL) && isempty(cumprodM)
                cumprodL = 1;
                cumprodM = 1;
            end    
            Rmpos(i) =  (((Order((l-1)*length(permsFL) + m,i) + 1)*permsFM(m,i))/permsFL(l,i)) * (cumprodL(end)/cumprodM(end));

            %Creating the P-M filters
            %H{i} = dfilt.dfsymfir(firpm(Order(k,i),fo,ao,w)); %Folded implementation, reduces multiples by 2
            %fvtool(H{i})
            %Need to adapt the input frequency after passing through each stage
            Fs = (Fs * permsFL(l,i))/permsFM(m,i);
            
            
            
            %If we want to use the spectral factorisation, the filter has to
            %have a nonnegative zerophase reponse, which implies, has to have
            %an even order
        
            if rem(Order((l-1)*length(permsFL) + m,i),2)~= 0
                Order((l-1)*length(permsFL) + m,i) = Order((l-1)*length(permsFL) + m,i) + 1;
            end
            
            b = firpm(Order((l-1)*length(permsFL) + m,i),fo,ao,w); 
            
            %fvtool(b)
            
            %Now, we create the Schuessler filters
            
            H_schuessler = schuessler(b,Delta2);
            
            if H_schuessler == 0  
                continue
            end    
            
            Order((l-1)*length(permsFL) + m,i) = length(H_schuessler) - 1;
            
            
            %Need to create a flag to inform that a Schuessler form has
            %been used
            isSchuessler((l-1)*length(permsFL) + m,i) = 1;
            
            %H{i} = dfilt.dffirt(schuessler(b,Delta2));
            
            %fvtool(H{i}) 
        end
        
    
        
 

        %Computing the total number of MADS for each permutation
        
        %We first need to remove the lines where the Order have 0
        
        if wrong_stage == 1
            %MADS((l-1)*length(permsFL) + m) = 0;
            MPOS((l-1)*length(permsFL) + m) = 0;
            permutation((l-1)*length(permsFL) + m) = 0;
            
        else
            
        [A,B] = deal(0);
        for i = 1:length(FL)
            %A = A + Rmads(i); 
            B = B + Rmpos(i);
        end

        %MADS((l-1)*length(permsFL) + m) = A;
        MPOS((l-1)*length(permsFL) + m) = B;
        
        end
        
        
        
        if manual == 3
            Rp = Rp*length(FM); %Otherwise, it keeps dividing Rp at each for loop
        end    
        
        
        
        end
    
    end
    
end

end

%Some of the orders can be 0 (negative stopbands), we need to remove them
Order(~all(Order,2),:) = [];
%MADS(~all(MADS,2)) = [];
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
end



end