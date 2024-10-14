%Gauss elimination method
function[X]=Gauss(A,b);
    now=[A,b];
    [row,col]=size(A);
    for i=1:row
        if now(i,i)~=0
            for j=i+1:row
                t=now(j,i)/now(i,i);
                now(j,:)=now(j,:)-t*now(i,:);
            end
            now(i,:)=now(i,:)/now(i,i);
        
        else
        fprintf("There has to be some problem");
        end
    end
    [X]=solution(now)';


end