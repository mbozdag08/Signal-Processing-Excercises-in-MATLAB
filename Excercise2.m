%Part 1.1 - DTMF Transmitter
telnum = [0 5 3 1 8 4 5 0 3 4 9];
x1 = DTMFTRA(telnum);
%soundsc(x1, 8192);

%Part 1.2 - DTMF Receiver
idnum = [3 6 1 2]; %ID: 21702163
x2 = DTMFTRA(idnum);
%soundsc(x2, 8192);

X2 = FT(x2); %Fourier Transform of x2(t)
omega=linspace(-8192*pi,8192*pi,8193); %Given omega
omega=omega(1:8192);
%Ploting the dialed number 3612
figure;
plot(omega, abs(X2));xlabel('Omega (w)');
ylabel('X(jw)');
title('Fourier Transform of Reverse ID Digits ([3 6 1 2])')

t = [0 : 1/8191 : 1];
%Ploting each digit seperately to see the dialed numbers
for k = 1 : 4 
    recx2{k} = x2 .* rectangularPulse( (k-1) * 0.25 , k * 0.25, t );
    recX2{k} = FT(recx2{k});
    figure;
    plot(omega,abs(recX2{k}));
    xlabel('Omega (w)'); 
    ylabel(sprintf('|X_%d(w)|', k)); 
    title(sprintf('Magnitude of X_%d(w)', k));
end

%Part 2
[y, Fs] = audioread('ses.mp3'); %Reading the audio file
[P,Q] = rat(8192/Fs);
x = resample(y,P,Q); %Resampling the audio recording to satisfy the criteria
x = x(10001:91920,:); %Croping the signal to 81920 samples
x=transpose(x);

t = [0 : 1/8192 : 10-1/8192]; %Given values
M = 5;
Ai = [0.35 0.5 0.65 0.05 0.15];
ti = [0.50 1.25 2.5 3 2.75];

%soundsc(x); %Listening to the original signal
figure;
plot(t,x)
xlabel("Time (seconds)");
ylabel("x(t)")
title("Original Speech Signal");

%Adding the echoes
y = x; 
for k=1:M 
    t0 = ti(k) * 8192; 
    xi = x(:,[1 : (round((10 - ti(k)) * 8192) + 1) ]);
    y(t0:end) = y(t0:end) + Ai(k) * xi; 
end 
%soundsc(y); %Listening the echoed signal
figure;
plot(t,y)
xlabel("Time (seconds)");
ylabel("y(t)")
title("Speech Signal with Echoes");

Y = FT(y); %Fourier Transform of y(t)

omega=linspace(-8192*pi,8192*pi,81921); %Given omega
omega=omega(1:81920);

H = ones(size(omega)); %Frequency Response H(jw)
for k = 1:M 
     H = H + Ai(k) * exp(-1i * omega * ti(k));
end 

h = IFT(H); %Impulse Response h(t)
figure;
plot(t,h);
xlabel('Time (seconds)'); 
ylabel('h(t)'); 
title('Impulse Response');
figure;
plot(omega,real(H));
xlabel('w (rad/sec)'); 
ylabel('H(jw)'); 
title('Frequency Response');
Xe = Y ./ H; %Fourier Transform of xe(t), the Estimated Speech Signal
xe = IFT(Xe);
soundsc(xe); %Listening the estimation of the original signal

figure;
plot(t,xe);
xlabel('Time (seconds)'); 
ylabel('xe(t)'); 
title('Estimated Speech Signal');


function [x] = DTMFTRA(Number)
%DTMFTRA : Double Tone Multi Frequency Transmitter
%   Number: The pressed numbers in a DTMF Keypad, 1xN array
%   x: The output signal, x(t);

Ts = 1/8192;
t = 0 : Ts : 0.25-Ts;
x = zeros(1, length(t).*length(Number));

 for k = 1 : length(Number)
     
     if Number(k) == 1
         y = cos(2*pi*697*t) + cos(2*pi*1209*t);
     elseif Number(k) == 2
         y = cos(2*pi*697*t) + cos(2*pi*1336*t);
     elseif Number(k) == 3
         y = cos(2*pi*697*t) + cos(2*pi*1477*t);
     elseif Number(k) == 4
         y = cos(2*pi*770*t) + cos(2*pi*1209*t);
     elseif Number(k) == 5
         y = cos(2*pi*770*t) + cos(2*pi*1336*t);
     elseif Number(k) == 6
         y = cos(2*pi*770*t) + cos(2*pi*1477*t);
     elseif Number(k) == 7
         y = cos(2*pi*852*t) + cos(2*pi*1209*t);
     elseif Number(k) == 8
         y = cos(2*pi*852*t) + cos(2*pi*1336*t);
     elseif Number(k) == 9
         y = cos(2*pi*852*t) + cos(2*pi*1477*t);
     elseif Number(k) == 0
         y = cos(2*pi*941*t) + cos(2*pi*1336*t);
     end
     
     %Rearranging t and x
     t = k*0.25 : Ts : (k+1)*0.25 -Ts;
     if k == 1
         x = y;
     else
     x = [x , y];
     end
     
 end
 
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