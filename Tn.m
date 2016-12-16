clear
close all
import_file
figure
B=cell2mat(B);
B=str2double(B);
C=A(2:end,1);
T=reshape(C,A(1),B);
Tn1=flip(T);
contourf(Tn1)

import2
figure
D=cell2mat(VarName2);
D=str2double(D);
E=VarName1(2:end,1);
Y=reshape(E,VarName1(1),D);
Yn=flip(Y);
contourf(Yn)

figure
diff_T=Tn1-Yn;
total_diff=sum(sum(diff_T));
contourf(diff_T)
