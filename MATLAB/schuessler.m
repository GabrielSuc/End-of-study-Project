function H_schuess = schuessler(b,delta2)
    
%This function aims to create a minimum-phase FIR filter according to
%Schuessler's method. It turns an FIR filter (Parks-McClellan) into this
%new type of filter.

%STEP ONE
%Creating the Parks-McClellan with the required specs i.e. those to produce
%the schuessler filter with the expected spec Rp, Rs

%Done directly in myMultistage.m


% STEP TWO
%Step two is to adjust the filter b so that the filter can  be  
%spectrally  factored. 

%'lifting' the amplitude response

h1 = b + delta2;

  
% STEP THREE
%We can obtain a minimum-phase filter by spectrally factoring  the  filter
%h1

%However, h1 has to have a nonnegative zero-phase reponse, odd length
%can resolve the pb

if min(zerophase(h1)) < -1e-3
    H_schuess = 0;
    return
end
   

h2 = firminphase(h1);


% STEP FOUR
%Finally,  we  need  only  scale h3(n) so  that  it  has  an equi-ripple
%behavior in the pass-band

H_schuess = h2./sqrt(1+delta2);

end