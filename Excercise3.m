%Part 1

A = ReadMyImage('Part5.bmp');
DisplayMyImage(A);
clear;

%Part 2

A = zeros(3,3);
A(1:2 ,1) = [1,2];
A(3:3,1:1) = [3];
A(1:1,2:3) = [4 7];
B = A(1:2,2:3);
B = A(3:3,1:2);
clear;

x = [8,1,6 ; 3,5,7 ; 4,9,2]; %Input
h = [1,3 ; 4,2];             %Impulse Response
y = DSLSI2D(h,x)             %Output: y = x ** h
clear;

%Part 3

D = 21702163;       %My Student ID
D17 = rem(D,17);    %Modulo 17 of my ID = 14
Mh = 20 + D17;      %Height of h
Nh = 20 + D17;      %Width of h

subplot(2,2,1);     %Creating a subplot
x = ReadMyImage('Part4.bmp');
imshow(x,[]);       %Reading and Displaying the Original Image
title("Original Image");

for i = 1 : 3       %2D convoluiton with 3 different h for 3 different B values
    
    if i == 1;
        B = 0.8;
    elseif i == 2;
        B = 0.5;
    else
        B = 0.2;
    end
    
    for m = 1 : Mh  %Creating the impulse response h according to B
        for n = 1 : Nh
            h(m, n) = sinc(B * (m - (Mh-1)/2)) .* sinc(B * (n - (Nh-1)/2));
        end
    end
    
    y(:,:,i) = DSLSI2D(h,x); %Saving the output for differnet B values
    
    subplot(2,2,i+1);
    im = y(:,:,i);
    imshow(im,[]);
    title("Denoised Image with B = " + string(B));
    
end
clear;

%Part 4

x = ReadMyImage("Part5.bmp"); %Reading the image

h1 = [1 -1];        %Definig h1 and h2, finding y1 = h1 ** x and y2 = h2 ** x
y1 = DSLSI2D(h1,x);
h2 = [1; -1];
y2 = DSLSI2D(h2,x);

figure;
subplot(1,2,1);
imshow(y1,[]);
title("y1");

s1 = y1 .* y1;
subplot(1,2,2);
imshow(s1,[]);
title("s1");

figure;
subplot(1,2,1);
imshow(y2,[]);
title("y2");

s2 = y2 .* y2;
subplot(1,2,2);
imshow(s2,[]);
title("s2");

h3 = 0.5 .* h1 + 0.5 .* h2; %Calculating h3 and finding y3 = h3 ** x
y3 = DSLSI2D(h3,x);

figure;
subplot(1,2,1);
imshow(y3,[]);
title("y3");

s3 = y3 .* y3;
subplot(1,2,2);
imshow(s3,[]);
title("s3");

figure;     %Showing all edges along with the original image
subplot(2,2,1);
imshow(x);
title("Original Image");
subplot(2,2,2);
imshow(s1);
title("Vertical Edges");
subplot(2,2,3);
imshow(s2);
title("Horizontal Edges");
subplot(2,2,4);
imshow(s3);
title("All Edges");

%Part 5

x = ReadMyImage("Part6x.bmp"); %Reading the image
figure;
DisplayMyImage(x);
title("National Soccer Team");

h6 = ReadMyImage("Part6h.bmp"); %Reading the image
figure;
DisplayMyImage(h6);
title("Inverted Face");

y = DSLSI2D(h6,x);
yabs = abs(y);
yabs3 = yabs .^3;
yabs5 = yabs .^5;
figure;
DisplayMyImage(yabs);
title("Output")

figure;
DisplayMyImage(yabs3);
title("3rd Power of Output: |y|^3")

figure;
DisplayMyImage(yabs5);
title("5th Power of Output: |y|^5") 



function [y] = DSLSI2D(h,x)
%DSLSI2D is the 2D convolution of h and x: x ** h = y
%   h: 2D Impulse Response of size Mh x Nh
%   x: 2D Input of size Mx x Nx
%   y: 2D Output of size (Mx + Mh - 1) x (Nx + Nh - 1) = My x Ny

[Mh, Nh] = size(h); %Taking the dimensions of x
[Mx, Nx] = size(x); %Taking the dimensions of h
y = zeros(Mx+Mh-1, Nx+Nh-1); %Creating an emty array y of size My x Ny

for k = 0 : Mh - 1
    for l = 0 : Nh - 1
        y(k+1 : k+Mx, l+1 : l+Nx) = y(k+1 : k+Mx, l+1 : l+Nx) + ( h(k+1, l+1) * x );
    end
end
end

function [x] = ReadMyImage(string)
x=double((rgb2gray(imread(string))));
x=x-min(min(x));
x=x/max(max(x));
x=x-0.5;
end

function []=DisplayMyImage(Image)
Image=Image-min(min(Image));
figure;
imshow(uint8(255*Image/max(max(abs(Image)))));
end