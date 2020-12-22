%Ishita Pal  (IXP180006)

clc;
clear all;
close all;
M = 16;  %modulation index
k = log2(M); 
symbols = [-3,-1,1,3];
EbNo = [6:2:12];
EsNo = k*EbNo;
N = 100000; %number of symbol
S = zeros(1,N);

for i = 1:length(EbNo)
     
    s = randsrc(1,N,symbols) + j*randsrc(1,N,symbols); 
    n = (1/sqrt(10))*s;  %normalizing 
    a = (1/sqrt(2))*[randn(1,N) + j*randn(1,N)];
    y = n + 10^(-EbNo(i)/20)*a;  %adding white guassian noise
    
    %demodulation
    y_x = real(y);
    y_y = imag(y);
    
    S_re(find(y_x< -2/sqrt(10))) = -3;
    S_re(find(y_x>= -2/sqrt(10) & y_x<0)) = -1;
    S_re(find(y_x >=0 & y_x<2/sqrt(10))) = 1;
    S_re(find(y_x >= 2/sqrt(10))) = 3;
    
    S_im(find(y_y>= 2/sqrt(10))) = 3;
    S_im(find(y_y<2/sqrt(10) & y_y>=0)) = 1;
    S_im(find(y_y<0 & y_y>=-2/sqrt(10))) = -1;
    S_im(find(y_y<-2/sqrt(10))) = -3;
    
    S = S_re + j*S_im;
    TotError(i) = size(find(s-S),2); %calculating total error
end     
    
BER = TotError/N; % calculating bit error rate
SER = k*BER; %calculation sysmbol error rate

%displaying the BER and SER
disp('Bit error rate');
disp(BER);
disp('symbol error rate');
disp(SER);
    
%calualting using theortitical formula
BERtheory = 3/2*erfc(sqrt(0.1*(10.^(EbNo/10))));
SERtheory = k * BERtheory;

%plotting BER graph
figure;
semilogy(EbNo,BERtheory,'-','Linewidth',2);
hold on;
semilogy(EbNo,BER,':*','Linewidth',2);
grid on;
legend('Theoritical', 'Simulation');
xlabel('Eb/No  (dB)')
ylabel('Bit Error Rate')
title('Bit error probability for 16-QAM modulation');

%plotting SER graph
figure;
semilogy(EsNo,SERtheory,'-','Linewidth',2);
hold on;
semilogy(EsNo,SER,':*','Linewidth',2);
grid on;
legend('Theoritical', 'Simulation');
xlabel('Es/No  (dB)')
ylabel('Symbol Error Rate')
title('Symbol error probability for 16-QAM modulation');

