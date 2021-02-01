% Part 1

load('Part2.mat')   % Loading the given data
figure;             % Plotting r(t) against t (time)
plot(t,r);
xlabel("t (seconds)");
ylabel("r(t)");
title("Samples of r(t)")

R = FT(r);          % Taking the Fourier Transform of r(t)
figure;             % Plotting the magnitude of the Fourier Transform R(w)
plot(omega,abs(R));
xlabel('\omega (rad/s)');
ylabel("|R(\omega)|"),
title("Fourier Transform of r(t) = R(\omega)")

ListenToMyRadio(r,55000,t,omega);   % Listening and plotting the signal at fc = 55kHz


% Part 2

load('Part2.mat')   % Loading the given data

% Part 3

R = FT(r);                  % Taking the Fourier Transform of r(t)
H = ones(1,length(omega));  % Creating the rectangular function for Band-Pass Filtering
H(abs(abs(omega) - 2*pi*40000) > 2*pi*4100) = 0;  
W = R .* H;                 % Band-Pass Filtered Signal W(w)
w = IFT(W);                 % The Inverse Fouirer Transfomr of W(w)

figure;
plot(omega,abs(H));
title("Band-Pass Filter H(\omega)")
xlabel("\omega (rad/s)");
ylabel("|H(\omega)|")

% Part 4

figure;
plot(omega,abs(W));
title("Fourier Transform of w(t) = W(\omega) (Band-Pass Filtered R(\omega))")
xlabel("\omega (rad/s)");
ylabel("|W(\omega)|")

% Part 5

ListenToMyRadio(w,40000,t,omega); % Listening and Plotting the signal at fc = 40kHz

% Part 6

H1 = ones(1, length(omega));        % Creating the rectangular function for Low-Pass Filtering
H1(abs(omega) >= 2*pi*40000) = 0;
Q = W .* H1;                        % Low-Pass Filtered W(w)
q = IFT(Q);                         % The Inverse Fouirer Transfomr of Q(w)

figure;
plot(omega,abs(H1));
title("Low-Pass Filter H1(\omega)")
xlabel("\omega (rad/s)");
ylabel("|H1(\omega)|")

% Part 7

figure;
plot(omega,abs(Q));
title("Fourier Transform of q(t) = Q(\omega) (Low-Pass Filtered W(\omega))")
xlabel("\omega (rad/s)");
ylabel("|Q(\omega)|")

% Part 8

ListenToMyRadio(q,40000,t,omega);   % Listening and Plotting the signal at fc = 40kHz


function [] = ListenToMyRadio(r,fc,t,omega)
%LÝSTENTOMYRADÝO AM Receiver
%   r: Sampled signal in time domain
%   fc: Frequency of the radio channel (Hz)
%   t: Time
%   omega = Angular Frequency (radian/seconds)

d = r .* cos(2*pi*fc*t);    % Modulating the sampled signal with fc
D = FT(d);                  % Taking the Fourier Transform of the modulated signal
M = D;                      % Filtering the signal D in frequency domain
M(abs(omega)>=8200*pi) = 0;
m = IFT(M);                 % Taking the Inverse Fourier Transform of the filtered signal
soundsc(m,200000);          % Listening the signal

figure;                     % Plotting d(t), |D(w)|    
subplot(1,2,1);
plot(t,d);
title("Modulated Signal d(t)");
xlabel("t (seconds)");
ylabel("d(t)");

subplot(1,2,2);
plot(omega,abs(D));
title("Fourier Transform of d(t) = D(\omega)");
xlabel("\omega (rad/s)");
ylabel("|D(\omega)|");

figure;                     % Plotting |M(w)|, m(t)   
subplot(1,2,1);
plot(omega,abs(M));
title(sprintf("Filtered Signal M(\\omega) for fc = %d Hz", fc));
xlabel("\omega (rad/s)");
ylabel("|M(\omega)|");

subplot(1,2,2);
plot(t,m);
title(sprintf("Signal of Station with fc = %d Hz", fc));
xlabel("t (seconds)");
ylabel("m(t)");

end

function output=FT(input)
M=size(input,2);
t=exp(j*pi*(M-1)/M*[0:1:M-1]);

output=exp(-j*pi*(M-1)^2/(2*M))*t.*1/(M)^0.5.*fft(input.*t);
end

function output=IFT(input)
M=size(input,2);
t=exp(-j*pi*(M-1)/M*[0:1:M-1]);
output=real(exp(j*pi*(M-1)^2/(2*M))*t.*(M)^0.5.*ifft(input.*t));
end

