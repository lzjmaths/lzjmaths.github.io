function X=QR_Givens(A,b);
[row,col]=size(A);

now=[A,b];
for i=1:col
    for j=i+1:row
        if now(j,i)==0;
            continue
        end
        m=sqrt(now(i,i)^2+now(j,i)^2);
        c=now(i,i)/m;
        s=now(j,i)/m;
        tmp=now(i,:);
        now(i,:)=now(i,:)*c+s*now(j,:);
        now(j,:)=-tmp*s+c*now(j,:);
    end
    now(i,:)=now(i,:)/now(i,i);
end
X=solution(now);
end

