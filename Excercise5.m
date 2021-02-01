% Part 1

D = 21702163;
D5 = rem(21702163,5);    % D5 = 3
M = 4 + D5;              % The lenght of b[k]x[n-k] terms

a = 0;                   % Coefficients of y[n-l] terms
k = 0 : M-1;             % Array for k, index of x[n-k]
x = [1, zeros(1,10)];    % x[n] = Unit Impulse
b = exp(-k);             % Coefficients of x[n-k] terms
Ny = length(x);          % y[n] only consists of x[n] components
hn = DTLTI(a, b, x, Ny); % Calculating the Impulse Response
n = 0 : 10;              % Index for the plot
figure;                  % Plotting the Impulse Response
stem(n, hn);             
xlabel("n");
ylabel("h[n]"),
title("Impulse Response of the Filter");


% Part 2

omega = -pi:0.001:pi;       % Omega and H(e^j.omega) 
hOmega = (1 - exp(-8.* (1 + 1i*omega))) ./ (1 - exp(-1 - 1i*omega));
maghOmega = abs(hOmega);    % Magnitude of H(e^j.omega)
figure;                     % Plotting the DTFT of the Impulse Response
plot(omega,maghOmega);      
xlabel("\Omega (rad/s)")
xticks([-pi -3*pi/4 -pi/2 -pi/4 0 pi/4 pi/2  3*pi/4 pi]);
xticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
ylabel("|H(e^{j\Omega})|")
title("Magnitude of the Calculated DTFT of the Impulse Response")


% Part 3

fs = 2000;              % Sampling Frequency = 2kHz
tVals = [1 10 1000];    % Given intervals for t
n = 1;
figure;

for t0 = tVals
    t = 0 : 1/fs : t0;                   % t of the Chirp Signal
    k = 1000/t0;                         % k for the Chirp Signal (finst = kt)
    xch = cos(2*pi*(k/2)*t.^2);          % The Chirp Signal
    ych = DTLTI(a, b, xch, length(xch)); % Output of the Chirp Signal
    
    subplot(3,2,n);
    plot(linspace(-pi,pi,length(xch)), abs(xch));
    title(sprintf('The Chirp Signal for 0 < t < %d', t0));
    ylabel('|X|(e^{j\Omega})');
    xlabel('\Omega (rad/s)');
    xlim([-pi pi]);
    xticks([-pi -pi/2 0 pi/2 pi]);
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    n = n + 1;
    
    subplot(3,2,n);
    plot(linspace(0,pi,length(ych)), abs(ych));
    title(sprintf('Output of the Chirp Signal for 0 < t < %d', t0));
    ylabel('|Y|(e^{j\Omega})');
    xlabel('\Omega (rad/s)');
    xlim([0 pi]);
    xticks([0 pi/4 pi/2  3*pi/4 pi]);
    xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'});
    n = n + 1;
end

% Part 4

z1 = 0.640 + 1j*0.768;
p1 = 0.924 + 1j*0.370;
p2 = 0.507 + 1j*0.845;
omega = -pi:0.01:pi;
H = (exp(1j*omega) - z1) ./ ((exp(1j*omega) - p1).*(exp(1j*omega) - p2));
absH = abs(H);

figure;
plot(omega,absH);
title("DTFT of the Filter");
ylabel("|H(e^j\Omega)|")
xlabel("\Omega");
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels(["-\pi", "-\pi/2", "0", "\pi/2", "\pi"]);


% Part 5

a = [(p1+p2) -p1*p2];   % Found in Part 4b analytically
b = [0 1 -z1];

fs = 1000;              % Sampling Frequency = 1kHz
tVals = [1 10 1000];    % Given intervals for t
n = 1;
figure;

for t0 = tVals
    t = 0 : 1/fs : t0;                           % t of the Chirp Signal
    k = 1000/t0;                                 % k for the Chirp Signal (finst = kt)
    xch = exp(2j*pi*(-500*t + (k/2)*t.^2));      % The Chirp Signal
    ych = DTLTI(a, b, xch, length(xch));         % Output of the Chirp Signal
    
    subplot(3,2,n);
    plot(linspace(-pi,pi,length(ych)), abs(ych));
    title(sprintf('Magnitude of the Output for 0 < t < %d', t0));
    ylabel('|Y|(e^{j\Omega})');
    xlabel('\Omega (rad/s)');
    xlim([-pi pi]);
    xticks([-pi -pi/2 0 pi/2 pi]);
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    n = n + 1;
    
    subplot(3,2,n);
    plot(linspace(-pi,pi,length(ych)), angle(ych));
    title(sprintf('Phase of Output for 0 < t < %d', t0));
    ylabel('<Y(e^{j\Omega})');
    xlabel('\Omega (rad/s)');
    xlim([-pi pi]);
    xticks([-pi -pi/2 0 pi/2 pi]);
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    n = n + 1;
end


function [y] = DTLTI(a,b,x,Ny)
%DTLTI 
%   a: Coefficients of y[n-l], size 1xN
%   b: Coefficients of x[n-k], size 1xM
%   x: Input signal, size 1xNx
%   Ny: Size of the output signal y[n]

y = zeros(1,Ny); % Memory allocation
N = length(a);   % Length of coefficients of y[n-l]
M = length(b)-1;   % Length of coefficients of x[n-k]

for n = 1 : Ny
    for l = 1 : N       % Adding y[n-l] terms with coefficients a[l]
        if (n-l) >= 1
            y(n) = y(n) + a(l)*y(n-l); 
        end
    end
    for k = 0 : M       % Adding x[n-k] terms with coefficients b[k]
        if (n-k) >= 1
            y(n) = y(n) + b(k+1)*x(n-k); 
        end
    end
end
end



