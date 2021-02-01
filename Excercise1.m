%Part 1
t=[0:0.001:1];
n = mod(21702163, 41);
A = 3*rand(1,n) + 3j*rand(1,n);
omega = pi*rand(1,n);
xs = SUMCS(t,A,omega);
rxs = real(xs);
imxs = imag(xs);
amp = abs(xs);
phi = angle(xs);

figure;
plot(t,rxs);title('Real part of xs(t)');
xlabel('t');
ylabel('Re(xs(t))');

figure;
plot(t,imxs);
title('Imaginary part of xs(t)');
xlabel('t');
ylabel('Im(xs(t))');

figure;
plot(t,amp);
title('Amplitude of xs(t)');
xlabel('t');
ylabel('Amplitude');

figure;
plot(t,phi);
title('Phase angle of xs(t)');
xlabel('t');
ylabel('phi');


%Part 2
D11 = mod(21702163, 11); % D11 = 10
D4 = mod(21702163, 4); %D4 = 3
t = [-5:0.001:5];
T = 2; W = 1; K = 23 + D11;
xt = FSWave(t, K, T, W);

figure;
plot(t, real(xt));
xlabel('t');
ylabel('Re(x(t))');
title(sprintf('Real part of x(t) for K=%d', K));

figure;
plot(t, imag(xt));
xlabel('t'); 
ylabel('Im(x(t))'); 
title(sprintf('Imaginary part of x(t) for K=%d', K));

Kvals = [1+D4 7+D4 16+D4 100+D4 200+D4];

for K = Kvals % plot real part of x(t) for given K values
    figure; 
    xt = FSWave(t, K, T, W); 
    plot(t, real(xt)); 
    xlabel('t'); 
    ylabel('x(t)'); 
    title(sprintf('Real part of x(t) for K=%d', K)); 
end


%Part 3
t = [-5:0.001:5];
T = 2; W = 1; K = 23 + D11;

ytA = FSWave4a(t, K, T, W);

figure;
plot(t, real(ytA));
xlabel('t');
ylabel('Re(y(t))');
title(sprintf('Part 4a: Real part of y(t) for K=%d', K));

figure;
plot(t, imag(ytA));
xlabel('t'); 
ylabel('Im(y(t))'); 
title(sprintf('Part 4a: Imaginary part of y(t) for K=%d', K));

ytB = FSWave4b(t, K, T, W);

figure;
plot(t, real(ytB));
xlabel('t');
ylabel('Re(y(t))');
title(sprintf('Part 4b: Real part of y(t) for K=%d', K));

figure;
plot(t, imag(ytB));
xlabel('t'); 
ylabel('Im(y(t))'); 
title(sprintf('Part 4b: Imaginary part of y(t) for K=%d', K));

ytC = FSWave4c(t, K, T, W);

figure;
plot(t, real(ytC));
xlabel('t');
ylabel('Re(y(t))');
title(sprintf('Part 4c: Real part of y(t) for K=%d', K));

figure;
plot(t, imag(ytC));
xlabel('t'); 
ylabel('Im(y(t))'); 
title(sprintf('Part 4c: Imaginary part of y(t) for K=%d', K));
ytD = FSWave4d(t, K, T, W);

figure;
plot(t, real(ytD));
xlabel('t');
ylabel('Re(y(t))');
title(sprintf('Part 4d: Real part of y(t) for K=%d', K));

figure;
plot(t, imag(ytD));
xlabel('t'); 
ylabel('Im(y(t))'); 
title(sprintf('Part 4d: Imaginary part of y(t) for K=%d', K));

function [xs] = SUMCS(t, A, omega)
%SUMCS: The superposition of M complex exponentials.
%   t: 1xN vector that contains the time instants over which xs is computed.
%   A: 1xM complex-valued vector. i th element is Ai.
%   omega: 1xM vector. i th element is wi.
xs = zeros(1, length(t));
for i = 1:length(A)
    xs = xs + A(i).*exp(1j.*(omega(i)./2).*t);
end
end

function [xt] = FSWave(t, K, T, W)
%FSWAVE Fourier Sythesis (Series Expansion)
%   t: Time grid
%   K: Fourier series range (from -K to +K)
%   T: Fundamental Period 
%   W: W<T, variable of the function
xk = zeros(1, 2*K + 1);
omega = zeros(1, 2*K +1);
for k = -K:K
    if k == 0
        xk(k + K + 1) = (W/T)*(1/2 - (W^2)/24);
    else
        xk(k + K + 1) = ( (1/(pi*k)) - ((W^2)/(4*pi*k)) + ...
            ((2*T^2)/((pi*k)^3)) ) * sin((k*pi*W)/(2*T)) - ...
            ( W*T/(pi^2*k^2) ) * cos((k*pi*W)/(2*T));
    end
    omega(k + K + 1) = (2*pi*k)/T;
end
xt = SUMCS(t, xk, omega);
end

function [xt] = FSWave4a(t, K, T, W)
%FSWAVE Fourier Sythesis (Series Expansion)
%   t: Time grid
%   K: Fourier series range (from -K to +K)
%   T: Fundamental Period 
%   W: W<T, variable of the function
xk = zeros(1, 2*K + 1);
omega = zeros(1, 2*K +1);
for k = -K:K
    if k == 0
        xk(k + K + 1) = (W/T)*(1/2 - (W^2)/24);
    else
        xk(k + K + 1) = ( (1/(pi*(-k))) - ((W^2)/(4*pi*(-k))) + ...
            ((2*T^2)/((pi*(-k))^3)) ) * sin(((-k)*pi*W)/(2*T)) - ...
            ( W*T/(pi^2*(-k)^2) ) * cos(((-k)*pi*W)/(2*T));
    end
    omega(k + K + 1) = (2*pi*k)/T;
end
xt = SUMCS(t, xk, omega);
end

function [xt] = FSWave4b(t, K, T, W)
%FSWAVE Fourier Sythesis (Series Expansion)
%   t: Time grid
%   K: Fourier series range (from -K to +K)
%   T: Fundamental Period 
%   W: W<T, variable of the function
xk = zeros(1, 2*K + 1);
omega = zeros(1, 2*K +1);
for k = -K:K
    if k == 0
        xk(k + K + 1) = (W/T)*(1/2 - (W^2)/24);
    else
        xk(k + K + 1) = ( ( (1/(pi*k)) - ((W^2)/(4*pi*k)) + ...
            ((2*T^2)/((pi*k)^3)) ) * sin((k*pi*W)/(2*T)) - ...
            ( W*T/(pi^2*k^2) ) * cos((k*pi*W)/(2*T)) ) * (exp((-1j*2*pi*k*0.65)/T));
    end
    omega(k + K + 1) = (2*pi*k)/T;
end
xt = SUMCS(t, xk, omega);
end

function [xt] = FSWave4c(t, K, T, W)
%FSWAVE Fourier Sythesis (Series Expansion)
%   t: Time grid
%   K: Fourier series range (from -K to +K)
%   T: Fundamental Period 
%   W: W<T, variable of the function
xk = zeros(1, 2*K + 1);
omega = zeros(1, 2*K +1);
for k = -K:K
    if k == 0
        xk(k + K + 1) = (W/T)*(1/2 - (W^2)/24);
    else
        xk(k + K + 1) = ( ( (1/(pi*k)) - ((W^2)/(4*pi*k)) + ...
            ((2*T^2)/((pi*k)^3)) ) * sin((k*pi*W)/(2*T)) - ...
            ( W*T/(pi^2*k^2) ) * cos((k*pi*W)/(2*T)) ) * ((1j*k*2*pi)/T);
    end
    omega(k + K + 1) = (2*pi*k)/T;
end
xt = SUMCS(t, xk, omega);
end

function [xt] = FSWave4d(t, K, T, W)
%FSWAVE Fourier Sythesis (Series Expansion)
%   t: Time grid
%   K: Fourier series range (from -K to +K)
%   T: Fundamental Period 
%   W: W<T, variable of the function
xk = zeros(1, 2*K + 1);
omega = zeros(1, 2*K +1);
for k = -K:K
    if k < 0
        g = -(K + 1 + k);
        xk(k + K + 1) = ( (1/(pi*g)) - ((W^2)/(4*pi*g)) + ...
            ((2*T^2)/((pi*g)^3)) ) * sin((g*pi*W)/(2*T)) - ...
            ( W*T/(pi^2*g^2) ) * cos((g*pi*W)/(2*T));
    elseif k == 0
        xk(k + K + 1) = (W/T)*(1/2 - (W^2)/24);
    elseif k > 0
        h = K + 1 - k;
        xk(k + K + 1) = ( (1/(pi*h)) - ((W^2)/(4*pi*h)) + ...
            ((2*T^2)/((pi*h)^3)) ) * sin((h*pi*W)/(2*T)) - ...
            ( W*T/(pi^2*h^2) ) * cos((h*pi*W)/(2*T));
    end
    omega(k + K + 1) = (2*pi*k)/T;
end
xt = SUMCS(t, xk, omega);
end


