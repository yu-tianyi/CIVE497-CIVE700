clc, clear, close all

%--------Exercise 001--------

for i = 1:10
    for j = 1:10
    M3(i,j) = 1;
    end
end

M3;

for k = 1:5
    for l = 1:5
        if k == l
            M5(k,l)=1;
        else
            M5(k,l)=0;
        end
    end
end

M5;

%--------Exercise 002--------

V = [1, 1, 1, 1, 1];
E2 = diag(V,1);

%--------Exxercise 003--------

M1 = reshape(1:64, 8, 8);

M2 = transpose(M1);

M3 = reshape(1:64, 8, 8);

M3(1:4,1:4) = 1;

M4 = reshape(1:64, 8, 8);

M4(4:5, 4:5) = 0;

M5 = reshape(1:64, 8, 8);
M_5 = repelem(eye(4, 4), 2, 2);
M5(M_5 == 1) = 1;

M6 = reshape(1:64, 8, 8);
M_6 = rot90(M_5);
M6(M_6 == 1) = 1;

M7 = reshape(1:64, 8, 8);
M_7 = repelem(eye(4, 4), 2, 2);
M7(M_7 == 1) = 0;

M8 = reshape(1:64, 8, 8);
M_8 = repelem(eye(4, 4), 2, 2);
M8(M_8 == 0) = 100

M9 = reshape(1:64, 8, 8);
M_9 = repelem(eye(4, 4), 2, 2);
M9(M_9 == 0) = 0