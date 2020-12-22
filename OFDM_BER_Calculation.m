%Ishita Pal (IXP180006)

clc;
clear all;
close all;
M = 4;  % Modulation 
k = log2(M);  % Bits/symbol
nFFT = 64;
nDSC = 52;  %number of subcarrier
nbitspersym = 104; %number of bits per symbol
N = 100; %number of symbols
l= nbitspersym*N;

EbNo = [2:2:8]; %bit to noise
EsNo = EbNo + 10*log10(nDSC/nFFT) + 10*log10(64/80); %symbol to noise

for i = 1:length(EbNo)
    
    %generating symbols
    si=2*(round(rand(1,l))-0.5); %In-phase symbol generation
    sq=2*(round(rand(1,l))-0.5); %Quadrature symbol generation
    s=si+j*sq;
    n=(1/sqrt(2*EbNo(i)))*(randn(1,l)+j*randn(1,l)); %Random noise generation
    r=s+n;
    grp = reshape(r,nbitspersym,N).'; %grouping the symbols
    
    %assigning modulated symbols to subcarrier
    M = [zeros(N,6) grp(:,[1:nbitspersym/2]) zeros(N,1) grp(:,[nbitspersym/2+1:nbitspersym]) zeros(N,5)];
    
    F = (nFFT/sqrt(nDSC))*ifft(fftshift(M.')).'; %taking FFT and normalizing the power
    
    F = [F(:,[49:64]) F]; %adding cyclic prefix
    
    F = reshape(F.',1,N*132);
    
    %adding white guassian noise
    A = 1/sqrt(2)*[randn(1,N*132) + j*randn(1,N*132)];
    y = sqrt(132/64)*F + 10^(-EsNo(i)/20)*A;
    
    y = reshape(y.',132,N).'; % formatting the received vector into symbols
    y = y(:,[17:80]); %removing cyclic prefix
    
    Y = (sqrt(nDSC)/nFFT)*fftshift(fft(y.')).';
    yMod = Y(:,[6+[1:nbitspersym/4] 7+[nbitspersym/4+1:nbitspersym/2] ]);
    si = sign(real(yMod)); %In-phase demodulation
    si = reshape(si,1,l/2).';
    sq = sign(imag(yMod)); %quadrature phase demodulation
    sq = reshape(sq,1,l/2).';
    BERi=(l-sum(si==si))/l; %In-phase BER calculation
    BERq=(l-sum(sq==sq))/l; %Quadrature BER calculation
    ber(i)= mean([BERi BERq])
end

%calculating BER
BER = ber/(N*nbitspersym); 
%theortical calculation
BERTheory = (1/2)*erfc(sqrt(10.^(EbNo/10)));    
%displaying the results
disp('BER:');
disp(BER);

%Plotting the graph
figure;
semilogy(EbNo,BERTheory,'-','Linewidth',2);
hold on
semilogy(EbNo,BER,':*','Linewidth',2);
grid on
legend('theortical', 'simulation');
xlabel('Eb/No  (dB)')
ylabel('Bit Error Rate')
title('Bit error probability curve for QPSK using OFDM')