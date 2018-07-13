function N = remlpord(f1,f2,d1,d2)
%See smarc open-source -> remez_lp.c

% Computation of FIR filter order by the Remez formula. 
% Based on "Optimum FIR Digital Filter Implementations for Decimation, 
% Interpolation, and Narrow-Band Filtering" by Crochiere and Rabiner,
% Equation (9).

A1A3 = [5.309e-03, 7.114e-02, -4.761e-01];
A4A6 = [-2.660e-03, -5.941e-01, -4.278e-01];

d1 = log10(d1);
d2 = log10(d2);

lh = [ d1 * d1, d1, 1];
DD1 = [0,0,0];

for  i = 1:3
    DD1(i) = d2*A1A3(i) + A4A6(i);
end        
        
D = 0.0;
for  i = 1:3
    D = D + DD1(i)*lh(i);
end        

fK = 11.01217 + (d1 - d2) * 0.51244;
df = f2 - f1;


N = D / df - fK * df + 1;



end