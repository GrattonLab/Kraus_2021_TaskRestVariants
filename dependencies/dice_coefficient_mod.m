function QS = dice_coefficient_mod(A,B)
%allpos = A+B;

C = sum(A & B);
QS = 2*C/(sum(A)+sum(B));

%C_vals = intersect(A,B);
%C = 0;
%for c = 2%1:length(C_vals)
%    C = C + nnz((A==C_vals(c)) & (B==C_vals(c)));
%end

%QS = 2*C/(length(A)+length(B));
