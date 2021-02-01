% Part 1

ID = 21702163;      % My Student ID
dur = rem(ID, 5);   % dur = 3
Ts = dur/3;         % Ts = 1
t = -dur/2 : Ts/1000 : dur/2 - Ts/1000; % t = time, interval given as Ts/1000

pZ = generateInterp(0,Ts,dur); % Zero-Order Hold
figure;
plot(t, pZ);
title("Zero-Order Hold Interpolation Function");
xlabel("t (seconds)");
ylabel("pZ(t)");

pL = generateInterp(1,Ts,dur); % Linear
figure;
plot(t, pL);
title("Linear Interpolation Function");
xlabel("t (seconds)");
ylabel("pL(t)");

pI = generateInterp(2,Ts,dur); % Ideal Bandlimited (sinc)
figure;
plot(t, pI);
title("Ideal Bandlimited Interpolation Function");
xlabel("t (seconds)");
ylabel("pI(t)");

% Part 2

a = randi([2,6]);       % a = random integer between 2 and 6
Ts = 1/(25*a);          % Ts given as 1/25a seconds
t = -5 : Ts : 5;        % g(t) for -5<t<5 with interval = Ts
n = -5/Ts : 1 : 5/Ts;   % g(nTs) for -5/Ts<n<5/Ts with interval = 1
g = zeros(size(t));     % Memory allocation for g(t)
g( t>= -1 & t < 0 ) = -1;      % g(t) = -1 for -1 <= t < 0
g( t> 0 & t <= 1 ) = 1;        % g(t) = 1 for 0 < t <= 1

figure;
stem(n,g);              % Stem plot of discrete g(t) which is equal to g(nTs)
title("Stem Plot of Discrete g(t) = g(nTs)")
xlabel("n")
ylabel("g(nTs)")

gRZ = DtoA(0, Ts, 10, g);   % Reconstructed g(t) with Zero-Order Hold Interpolation
gRL = DtoA(1, Ts, 10, g);   % Reconstructed g(t) with Linear Interpolation
gRI = DtoA(2, Ts, 10, g);   % Reconstructed g(t) with Ideal Bandlimited Interpolation

figure;
plot(linspace(-5,5,length(gRZ)), gRZ);
title("Zero-Order Hold Reconstruction of g(nTs)");
xlabel("t (seconds)");
ylabel("gR(t)")

figure;
plot(linspace(-5,5,length(gRL)), gRL);
title("Linear Reconstruction of g(nTs)");
xlabel("t (seconds)");
ylabel("gR(t)")

figure;
plot(linspace(-5,5,length(gRI)), gRI);
title("Ideal Bandlimited Reconstruction of g(nTs)");
xlabel("t (seconds)");
ylabel("gR(t)")

clear;

% Part 3

ID = 21702163;     % My Student ID
D5 = rem(ID, 5);   % D5 = 3

TsVals = [0.005*(D5+1), 0.25 + 0.01*D5, 0.18 + 0.005*(D5+1), 0.099, 0.15, 0.05, 0.01]; % Values of Ts

for Ts = TsVals
    
    t = -2 : Ts/1000 : 2;  % t for continuous indexing
    n = -2/Ts : 1 : 2/Ts;  % n for discrete indexing
    
    xt = 0.25*cos(2*pi*3*t + pi/8) + 0.4*cos(2*pi*5*t - 1.2) + 0.9*cos(2*pi*t + pi/4);           % Continuous x(t)
    xN = 0.25*cos(2*pi*3*n*Ts + pi/8) + 0.4*cos(2*pi*5*n*Ts - 1.2) + 0.9*cos(2*pi*n*Ts + pi/4);  % Discrete x[n]
    
    figure;
    plot(t, xt, 'b');
    hold on;
    stemPlot = stem(n*Ts, xN, 'r');
    hold off;
    xlabel("t (seconds)")
    ylabel('amplitude'); 
    title(sprintf('Plot of x[n] and x(t) with T_s=%g', Ts)); 
    legend('x(t)', 'x[n]', 'Location', 'northeast');
    
    figure;
    xRZ = DtoA(0, Ts, 4, xN);
    plot(linspace(-2,2,length(xRZ)), xRZ);
    xlabel("t (seconds)")
    ylabel('xR(t)'); 
    title(sprintf('Reconstruction of x(t) from x[n], Zero-Order Hold Interpolation, T_s=%g', Ts)); 
    
    figure;
    xRL = DtoA(1, Ts, 4, xN);
    plot(linspace(-2,2,length(xRL)), xRL);
    xlabel("t (seconds)")
    ylabel('xR(t)'); 
    title(sprintf('Reconstruction of x(t) from x[n], Lineer Interpolation, T_s=%g', Ts));
    
    figure;
    xRI = DtoA(2, Ts, 4, xN);
    plot(linspace(-2,2,length(xRL)), xRI);
    xlabel("t (seconds)")
    ylabel('xR(t)'); 
    title(sprintf('Reconstruction of x(t) from x[n], Ideal Bandlimited Interpolation, T_s=%g', Ts));
    
end


function [xR] = DtoA(type,Ts,dur,Xn)
%DTOA Digital to Analog conversion
%   type: Interpolation Type (zero-order hold = 0, linear = 1, ideal bandlimited (sinc) = 2) 
%   Ts: Sampling rate
%   dur: Duration of the signal, symmetric wrt vertical axis (dur = 2, pulse will be between -1 and 1)
%   Xn: Samples of the signal x(t), size 1 x N, Xn(n) = x((n-1)Ts)

p = generateInterp(type,Ts,dur); % Generating the interpolation function
N = length(Xn);                  % Length of the samples Xn
Np = length(p);                  % Length of the interpolation function p
xR = zeros(1, N * 1000 + Np);    % Memory Allocation for xR, the reconstructed signal

for i = 1 : N
    xR( (i-1)*1000 + 1 : (i-1)*1000 + Np ) = xR( (i-1)*1000 + 1 : (i-1)*1000 + Np ) + p*Xn(i);
end
    xR = xR( 500*N + 1 : end - 500*N );
end

function [p] = generateInterp(type,Ts,dur)
%GENERATEINTERP Generates a function p(t) for interpolation 
%   type: Interpolation Type (zero-order hold = 0, linear = 1, ideal bandlimited (sinc) = 2)
%   Ts: Sampling rate
%   dur: Duration of the signal, symmetric wrt vertical axis (dur = 2, pulse will be between -1 and 1)

t = -dur/2 : Ts/1000 : dur/2 - Ts/1000; % t = time, interval given as Ts/1000
p = zeros(1, length(t));
if type == 0;
    p(t >= -Ts/2 & t < Ts/2) = 1;
elseif type ==1;
    p(t >= -Ts & t <= Ts) = 1 - abs(t(t >= -Ts & t <= Ts)) / Ts;
else
    p(t == 0) = 1;
    p = sin(pi * t/Ts) ./ (pi * t/Ts);
end

end

