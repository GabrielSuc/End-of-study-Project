function H_multi = myMultistage(L,M,Fx,Fy,Fsc,Fp,Rp,Rs,Wc,Fc)

%Multistage decomposition based on Proakis's method.
%According to the highest number of decomposition possible (M or L), 
%decompose a desired LP filter in multiple stages.   



%Ask for the type of filter we want to be returned

manual = input('1 to return a Parks-McClellan multistage filter, 2 for Elliptic, 3 for combination:');

    if isempty(manual)
          disp('Choose a filter');
    end
    
    
%Divide L & M in prime factors and sort them as decreasing numbers
% FM = sort(factor(M),'descend');
% FL = sort(factor(L),'descend');

[FM,FL] = getListStages(L,M);
 
FM = [3,7,7];

if length(FM) > length(FL) %Determine highest number of stages we'll do
   

%Have to inverse Fx and Fy considering the overall system as a decimator  
%or an interpolator
if M > L
     [Fy Fx] = deal(Fx,Fy);
end    

%The first factor can be too small(abs(F1-Fsc) < Fp), therefore we have to 
%modify the list of factors to overcome this issue
% while abs((Fx/FM(1))-Fsc)< Fp
%     
%     FM(1) = FM(1)*FM(end);
%     FM(end) = [];
%     
% end   
    
    Fs = zeros(length(FM)+1,1);
    FSC = zeros(length(FM),1);
    
    % ------------------------ Filters' parameters ------------------------
    
    %Sampling and cutoff frequencies of each of the stage
    Fs(1) = L*Fx;
    %Orders
    [Order1,Order2] = deal(zeros(length(FM),1));  
    %Final Filter
    H = cell(1,length(FM));

    for i = 2:(length(FM)+1)
        Fs(i) = Fs(i-1)/FM(i-1);  %(Fs(i-1)*FL(i-1))/(FM(i-1)); 
        FSC(i-1) = abs(Fs(i)-Fsc); %Don't know which one is bigger, depends on the width of the passband                                         
    end
    
    if manual == 1 % Parks-McClellan
    
        %Need of specific variables for the design
        a = [1 0];
        f = zeros(length(FM),length(a));
        
        
        %Ripples
        Delta1 = (10^((Rp/length(FM))/20)); %Reduce the passband ripple by factor of length(FM) so that passband
        Delta2 = (10^(-Rs/20));            %ripple in the cascade of the length(FM) filters doesn't exeed Delta1

          
        dev = [abs(Delta1 - 1)/(Delta1 + 1) Delta2]; 
        
        for i = 1:length(FM)
            %Defining the limit frequencies for the design
            f(i,:) = [Fp FSC(i)];
            Dinf = -(0.00266*((log10(Delta1))^2) + 0.5941*log10(Delta1))
            Order2(i) = 
            %Creating the P-M filters
            [Order1(i),fo,ao,w] = firpmord(f(i,:),a,dev,Fs(i));
            H{i} = dfilt.dffirt(firpm(Order1(i),fo,ao,w)); %One implementation among others
            fvtool(H{i})
        end  
        
        %Cascading them
        H_multi = [Order1;Order2];%dfilt.cascade(H{1,1:end});
        return

    elseif manual == 2 %Elliptic
        
        Wsc = zeros(length(FM),1);
        Wp = Wc*(Fp/Fc);%/pi;
        coef = cell(1,3);

        %Need to adapt the ripple in the passband
        Rp = Rp/length(FM); 
        
        for i = 1:length(FM)
            % Discrete domain equivalent
            Wsc(i) = Wc*(FSC(i)/Fc);%/pi;
            %Order
            Order(i) = ellipord(Wp,Wsc(i),Rp,Rs);
            %Creating the Elliptic Filters
            [coef{1,1},coef{1,2},coef{1,3}] = ellip(Order(i),Rp,Rs,Wp); %[z,p,k]
            H{i} = dfilt.df2sos(zp2sos(coef{1,1},coef{1,2},coef{1,3})); %Avoid error of quantization
            fvtool(H{i})
        end    
        
        %Cascading them
        H_multi = dfilt.cascade(H{1,1:end});
        return
        
    elseif manual == 3 %Any combination
        
        %First filter is an IIR Elliptic filter and the rest is Schuessler
        %filters
        
        %Elliptic Filter
        Wsc = Wc*(FSC(1)/Fc)/pi;
        Wp = Wc*(Fp/Fc)/pi;
        coef = cell(1,3);
        
        %Need to adapt the ripple in the passband
        Rp = Rp/length(FM); %Faut convertir en lineaire
        
        Order(1) =  ellipord(Wp,Wsc,Rp,Rs);
        
        [coef{1,1},coef{1,2},coef{1,3}] = ellip(Order(1),Rp,Rs,Wp); %[z,p,k]
        H{1} = dfilt.df2sos(zp2sos(coef{1,1},coef{1,2},coef{1,3})); %Avoid error of quantization
        fvtool(H{1})
        
        %Design of the Schuessler filter
        %We want the same specs i.e. same Rp and Rs, thus we need to adapt
        %the specs of the PM filter first
        
        Delta2 = (10^-(Rs/20))^2/(2 - (10^-(Rs/20))^2);
        %Delta1 = (1 + Delta2)*((10^((length(FM)*Rp)/20) + 1)^2 - 1);        
        Delta1 = 10^(1/20);
        
        dev = [(Delta1 - 1)/(Delta1 + 1) Delta2];

        a = [1 0];
        f = zeros(length(FM) - 1,length(a));

        for i = 2:length(FM) %Don't forget that the first filter is an elliptic

            %Defining the limit frequencies for the design
            f(i-1,:) = [(Fp*Wc)/pi (FSC(i)*Wc)/pi];
            %Creating the P-M filters
            [Order(i),fo,ao,w] = firpmord(f(i-1,:),a,dev,Fs(i));
            
            %If we want to use the spectral factorisation, the filter has to
            %have a nonnegative zerophase reponse, which implies, has to have
            %an even order
        
            if rem(Order(i),2)~= 0
                Order(i) = Order(i) + 1;
            end
            
            b = firpm(Order(i),fo,ao,w); 
            
            fvtool(b)
            
            %Now, we create the Schuessler filters
            H{i} = dfilt.dffirt(schuessler(b,Delta2));
            
            fvtool(H{i})
        end  
        
        %Cascading them
        H_multi = dfilt.cascade(H{1,1:end});
        return
        
        
    else
         disp('You did not choose an appropriate filter');
         
    end
    
   
        
   
   
   
elseif length(FL) > length(FM) %Determine highest number of stages we'll do
    
    %Determine highest number of stages we'll do
   

    %Have to inverse Fx and Fy considering the overall system as a decimator  
    %or an interpolator
    if M > L
    [Fy Fx] = deal(Fx,Fy);
    end    


    Fs = zeros(length(FL)+1,1);
    FSC = zeros(length(FL),1);
    
    % ------------------------ Filters' parameters ------------------------
    
    %Sampling and cutoff frequencies of each of the stage
    Fs(end) = Fx; %Fj in Proakis
    %Orders
    Order = zeros(length(FL),1);  
    %Final Filter                              
    H = cell(1,length(FL));
    
    %Since we are going backwards this time, we can flip the list FL
    FL = fliplr(FL);

    for i = length(FL):-1:1
        Fs(i) =  FL(i)*Fs(i+1);
        FSC(i) = Fs(i)-Fsc; %Don't know which one is bigger, depends on the width of the passband                                         
    end
    
    if manual == 1 % Parks-McClellan
    
        %Need of specific variables for the design
        a = [1 0];
        f = zeros(length(FL),length(a));
        
        
        %Ripples
        Delta1 = (10^((Rp/length(FL))/20)); %Reduce the passband ripple by factor of length(FM) so that passband
        Delta2 = (10^(-Rs/20));            %ripple in the cascade of the length(FM) filters doesn't exeed Delta1

          
        dev = [abs(Delta1 - 1)/(Delta1 + 1) Delta2]; 
        
        for i = length(FL):-1:1
            %Defining the limit frequencies for the design
            f(i,:) = [(Fp*Wc) (FSC(i)*Wc)];
            %Creating the P-M filters
            [Order(i),fo,ao,w] = firpmord(f(i,:),a,dev,Fs(i));
            H{i} = (firpm(Order(i),fo,ao,w)); %One implementation among others dfilt.dffirt
            %fvtool(H{i})
                   
        end  
        
        %Cascading them
        H_multi = Order;%dfilt.cascade(H{1,1:end});
        return
    
        
        elseif manual == 2 %Elliptic
        
        Wsc = zeros(length(FL),1);
        Wp = Wc*(Fp/Fc)/pi;
        coef = cell(1,3);

        %Need to adapt the ripple in the passband
        Rp = Rp/length(FL); 
        
        for i = length(FL):-1:1
            % Discrete domain equivalent
            Wsc(i) = Wc*(FSC(i)/Fc)/pi;
            %Order
            Order(i) = ellipord(Wp,Wsc(i),Rp,Rs);
            %Creating the Elliptic Filters
            [coef{1,1},coef{1,2},coef{1,3}] = ellip(Order(i),Rp,Rs,Wp); %[z,p,k]
            H{i} = dfilt.df2sos(zp2sos(coef{1,1},coef{1,2},coef{1,3})); %Avoid error of quantization
            fvtool(H{i})
        end    
        
        %Cascading them
        H_multi = dfilt.cascade(H{1,1:end});
        return

    else
            disp('You did not choose an appropriate filter');
         
    end
    
end    



end