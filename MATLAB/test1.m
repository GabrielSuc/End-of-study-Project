clear all
close all

fftsize = 1024;
genNew = 0;

if (genNew==1)
    x = rand(1024,1);
    fid2 = fopen('/home/gabriel/Documents/End-of-study-Project/MATLAB/PM_filter1.txt','wb');
    fwrite(fid2,single(x),'single');
    fclose(fid2);
else
    fid2 = fopen('/home/gabriel/Documents/End-of-study-Project/MATLAB/PM_filter1.raw','rb');
    x = fread(fid2,'single');
    fclose(fid2);
end

y = fft(x); 

fid3=fopen('/home/gabriel/Documents/End-of-study-Project/MATLAB/temp.raw','rb');
N=fftsize/2+1;
y_target = (zeros(N,1));
for (n=1:N)
    realpart=fread(fid3,1,'single');
    imagpart=fread(fid3,1,'single');
    y_target(n) = realpart + 1i*imagpart;
   
end
fclose(fid3);


figure(1)
plot(20*log10(abs(real(y_target)-real(y(1:N)))))

figure(2)
plot(20*log10(abs(imag(y_target)-imag(y(1:N)))))

