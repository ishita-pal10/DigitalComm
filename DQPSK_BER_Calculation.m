%Ishita Pal (IXP180006)

clc;
close all;
clear all;
N = 100000; % number of bits or symbols
M = 4; 
k = log2(M); % number of bits per symbol
EbNo = [4:2:10]; 
EsNo = EbNo + 10*log10(k);
for i = 1:length(EbNo)
     
     s = rand(1,N)>0.5;  % generating random binary signals

     grp = reshape(s,2,N/2).';
     bintodec = ones(N/2,1)*2.^[k-1:-1:0];
     dec_s = sum(grp.*bintodec,2);

     % converting to gray coded symbols
     gray_s = bitxor(dec_s,floor(dec_s/2));
     phasegray = 2*gray_s.'+1;
  
     % generating differential modulated symbols
     diffPhase = filter([ 1 ],[1 -1],phasegray); % start with 0 phase
     dqpsk_s = exp(j*diffPhase*pi/4);

     % white gaussian noise, 0 mean
     a = 1/sqrt(2)*[randn(1,N/2) + j*randn(1,N/2)];
     dqpsk = dqpsk_s + 10^(-(EsNo(i))/20)*a; % additive white gaussian noise

     %demodulation
     estphase = angle(dqpsk);
     estdphase = filter([1 -1],1,estphase)*4/pi;
     quant_dphase = 2*floor(estdphase/2)+1; % quantizing

     % gray to binary
     quant_dphase(find(quant_dphase<0))=quant_dphase(find(quant_dphase<0))+8;
     bin_phase = floor(bitxor(quant_dphase,floor(quant_dphase/2))/2);
     estBit = (dec2bin(bin_phase.')).';
     estBit = str2num(estBit(1:end).').';

     TotError(i) = size(find([s - estBit]),2); %calculating total error

     % theoretical BER computation
     a = sqrt(2*10.^(EbNo(i)/10)*(1-sqrt(1/2)));
     b = sqrt(2*10.^(EbNo(i)/10)*(1+sqrt(1/2)));
     k_b = 0:10;
     temp = exp(-((a.^2+b.^2)/2)).*sum((a/b).^k_b.*besseli(k_b,a*b));
     dqpskTheoritical(i) = temp - 0.5*besseli(0,a*b)*exp(-((a.^2+b.^2)/2));
end

%calculating BER
Ber = TotError/N;

disp('Bit error rate');
disp(Ber);

%plotting the BER
figure;
semilogy(EbNo,dqpskTheoritical,'-','Linewidth',2);
hold on
semilogy(EbNo,Ber,':*','Linewidth',2);
grid on
legend('theoritical', 'simulation');
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
title('Bit error probability curve for pi/4 DQPSK');
