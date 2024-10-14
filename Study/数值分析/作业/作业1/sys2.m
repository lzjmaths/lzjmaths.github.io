%test the Hilbert matrix, i.e T2
for i=1:40
    h(i)=0;
    for j=1:40
        H(i,j)=1/(i+j-1);
        h(i)=h(i)+H(i,j);
    end
end
b=h';
A=H;
tic;
fprintf("Gauss elemination method\n")
d1=Gauss(A,b);
norm(d1)
d1
toc;
tic;
fprintf("Gauss chosen elemination method\n")
d2=Gauss_choose(A,b);
norm(d2)
d2
toc;
tic;
fprintf("Cholesky method\n")
d3=Cholesky(A,b);
norm(d3)
d3
toc;