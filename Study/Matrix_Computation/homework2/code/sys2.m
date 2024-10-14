%test file
q3=sqrt(3);
q2=sqrt(2);

A=[(q2*x)/2 + 1/2,   3^(1/2)/2 - (3^(1/2)*q2*x)/6
;     -(q2*x)/2,               (3^(1/2)*q2*x)/6;
(q2*x)/2 - 1/2, - 3^(1/2)/2 - (3^(1/2)*q2*x)/6];



Ans=[1;sqrt(3)]


tic;

S=[       (3*x^2)/2 + 1/2, -(3^(1/2)*(x^2 - 1))/2;
-(3^(1/2)*(x^2 - 1))/2,            x^2/2 + 3/2]; %S=A^T*A
h=A'*b;
d1=Cholesky(S,h);
norm(d1-Ans)/norm(Ans)
d1
toc;

tic;

d1=Cholesky(S,A'*b);
norm(d1-Ans)/norm(Ans)
d1
toc;


tic;

d2=QR_Givens(A,b);
norm(d2-Ans)/norm(Ans)
d2
toc;