function    C = multiresolution_FT( A, n, pfilt, s, b, sigma )
% description: This code implements MFT for 2D signals (image)
% In this code none-overlapping windowed Fourier transform is considered
% Inputs 
% A: m*n matrix
% n: number of Laplacian pyramid level
% pfilt: pyramid filter name (see PFILTERS)
% s: size of the window (s*s matrix)
% b: b=0: boxing window for the windowed Fourier transform
% b=1: gaussian window for the windowed Fourier transform
% b=2: laplacian of gaussian window for the windowed fourier transform
% sigma: standard deviation of the gaussian function
% output
% C: MFT coefficients; cell array.

p = cell(1,n+1);
y = lpd(A,pfilt,n);
for i=1:n+1
    p{i} = size(y{i});
end
%resizing if needed (making square matrix)
for i=1:n+1
    l = p{i};
    A_1 = y{i}; 
    if  mod(l(1),s)~= 0
    A_1(l(1)+1:l(1)+s-mod(l(1),s),:) = 0;
    end
    if  mod(l(2),s)~= 0
    A_1(:,l(2)+1:l(2)+s-mod(l(2),s)) = 0;
    end
    g{i} = size(A_1);
    %making s*s Blocks
    A_2{i} = cell((g{i}(1)/s),(g{i}(2)/s));
    for x=0:g{i}(1)/s-1
        for j=0:g{i}(2)/s-1
        A_2{i}{x+1,j+1}=A_1(s*x+1:s*(x+1),s*j+1:s*(j+1));
        end
    end
end

switch b
    case 0
        for i=1:n+1
            B{i} = A_2{i};
        end
    case 1
        filt = fspecial('gaussian', s, sigma);
        for i=1:n+1
            for x=0:g{i}(1)/s-1
                for j=0:g{i}(2)/s-1
                    B{i}{x+1,j+1} = A_2{i}{x+1,j+1}.*filt;
                end
            end
        end
    case 2
        filt = fspecial('log', s, sigma);
        for i=1:n+1
            for x=0:g{i}(1)/s-1
                for j=0:g{i}(2)/s-1
                    B{i}{x+1,j+1} = A_2{i}{x+1,j+1}.*filt;
                end
            end
        end
end

for i=1:n+1
    c{i} = size(B{i});
end

for i=1:n+1
    z = c{i};
    for h=1:z(1)
        for a=1:z(2)
            C{i}{h,a} = fft2(B{i}{h,a});
        end
    end
end
end

