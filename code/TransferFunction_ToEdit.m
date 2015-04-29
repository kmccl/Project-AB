fft_orig = fft(orig,NFFT)/NFFT;
amp_orig = (2*abs(fft_orig(1:NFFT/2+1))).*NFFT;
figure,plot(amp_orig)
figure,plot(amp_orig)
figure,plot(f,amp_orig)
fft_lauren = fft(lauren,NFFT)/NFFT;
amp_lauren = (2*abs(fft_lauren(1:NFFT/2+1))).*NFFT;
hold on, plot(amp_lauren,'g')
hold on, plot(f,amp_lauren,'g')
figure,plot(f,amp_orig)
hold on, plot(f,amp_lauren,'g')
amp_ratio = amp_lauren./amp_orig;
figure,plot(amp_ratio)
figure,plot(f,amp_ratio)
hold on,plot(f,db2amp(interp_audiogram))
close
figure,plot(f,amp_ratio)
hold on,plot(f,db2amp(interp_audiogram),'g')

H = fft(zero_lauren)./fft(zero_orig);
figure,plot(abs(H))
Txy3=tfestimate(zero_orig,zero_lauren);
figure,plot(Txy3)
figure,plot(abs(Txy3))
help tfestimate

H = fft(lauren)/fft(orig);
realH = abs(H);
figure,plot(amp_filt_full)
zero_lauren = [zeros(length(lauren)+length(orig),1); lauren; zeros(length(lauren)+length(orig),1)];

figure,plot(zero_lauren)
zero_orig = [zeros(length(lauren)+length(orig),1); orig; zeros(length(lauren)+length(orig),1)];
H2 = fft(zero_lauren)/fft(zero_orig);
realH2 = abs(H2);
figure,plot([realH realH2])
size(realH)
size(realH2)
figure,plot(realH2)
figure,plot(realH2(find(realH2~=0)))
find(realH~=0)