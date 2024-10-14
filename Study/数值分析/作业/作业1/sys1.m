%, compute Ax=b for A fixed and x random.
 for i=2:99
     D(i,i)=10;
     D(i,i-1)=1;
     D(i,i+1)=1;
 end
 D(1,1)=10;
 D(1,2)=1;
 D(100,100)=10;
 D(100,99)=1;
 A=D;
x=rand(100,1);
b=A*x;
tic;
fprintf("Gauss elemination method\n")
s1=Gauss(A,b)-x;
norm(s1)/norm(x)
toc;
tic;
fprintf("Gauss chosen elemination method\n")
s2=Gauss_choose(A,b)-x;
norm(s2)/norm(x)
toc;
tic;
fprintf("Cholesky method\n")
s3=Cholesky(A,b);
s3=s3-x;
norm(s3)/norm(x)
toc;