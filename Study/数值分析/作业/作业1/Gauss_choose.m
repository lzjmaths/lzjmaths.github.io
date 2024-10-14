%Gauss elimination method by choosing row vector
function X=Gauss_choose(A,b);
    now=[A,b];
    [row,col]=size(A);
    for i=1:row
        for j=i+1:row
            if abs(now(i,i))<abs(now(i,i))
                temp=now(i,:);
                now(i,:)=now(j,:);
                now(j,:)=temp;
            end
        end
        
        if now(i,i)~=0
            for j=i+1:row
                t=now(j,i)/now(i,i);
                now(j,:)=now(j,:)-t*now(i,:);
            end
            now(i,:)=now(i,:)/now(i,i);
        
        else
        fprintf("There has to be some problem");
        end
    X=solution(now)';
end