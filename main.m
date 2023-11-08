close all
clear all
%%Parameters
len_data = 512;
pol=[0.1,0.42,0.9];
zer =[0.5];
SNR =10;
poles_ = [0.99*exp(pol(1)*pi*1i); 0.99*exp(-1*pol(1)*pi*1i);
          0.99*exp(pol(2)*pi*1i); 0.99*exp(-1*pol(2)*pi*1i)
          0.99*exp(pol(3)*pi*1i); 0.99*exp(-1*pol(3)*pi*1i)]; 
zeros_ = [ 0.97*exp(zer(1)*pi*1i); 0.97*exp(-1*zer(1)*pi*1i)
           ];
noise_flag = 1;
%% original spectrum
p = length(poles_);
q = length(zeros_);
b = poly(zeros_);
a = poly(poles_);
f = 0:0.001:0.5;
matrix_b = ones(1,length(f));
for i = 1:(length(b)-1)
    matrix_b = [matrix_b;matrix_b(i,:).*exp(-1i*2*pi*f)];
end
matrix_a = ones(1,length(f));
for i = 1:(length(a)-1)
    matrix_a = [matrix_a;matrix_a(i,:).*exp(-1i*2*pi*f)];
end
P = (abs(b*matrix_b).^2)./(abs(a*matrix_a).^2);
figure;
plot(f,10*log10(P));
xlabel(' Frequency(f/fs)');
ylabel('Magnitude (dB)');
title('PSD of Xn');
grid
%% spectrum estimation
for i=1:length(len_data)
Nn = randn(len_data(i),1);
Xn = filter(b,a,Nn);
    if noise_flag
    Ps = mean(Xn.^2);
    sigma_noise = sqrt(Ps/(10^(SNR/10)));
    noise = sigma_noise * randn(length(Xn),1);
    Xn = Xn + noise;
    end
p_algo = p+1;
q_algo = q+1;
[f, P_ARMA] = ARMA(Xn, p_algo, q_algo);
figure;
plot(f,10*log10(P));
xlabel(' Frequency(f/fs)');
ylabel('Magnitude(dB)');
hold on
plot(f,10*log10(real(P_ARMA)));
grid
xlabel(' Frequency(f/fs)');
ylabel('Magnitude(dB)');
p_str = num2str(p_algo/2);
q_str = num2str(q_algo/2);
len_str = num2str(len_data(i));
title_ = strcat(len_str,'samples');
title(title_);
legend('Original PSD','Estimated PSD');
end
