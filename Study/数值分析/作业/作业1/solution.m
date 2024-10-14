%solve upper triangular matrix
function [X]=solution(T);
    [row,col]=size(T);
    for i=row:-1:1
        X(i)=T(i,col);
        for j=i+1:col-1
            X(i)=X(i)-X(j)*T(i,j);
        end
    end
end