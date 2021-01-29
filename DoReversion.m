% function q=DoReversion(p,C,D) % jay matrix C va D ra dar f avaz mikonad
% q=p;
% q(1: C, :) = p(C+1:C+D, :);
% q(C+1: end, :) = p(1:C, :);
% %q=p([1:i1-1 i1+1:i2 i1 i2+1:end]);
% end

function q=DoReversion(p,i1,i2)
    q=p;
    if i1<i2
        q(i1:i2,:)=p(i2:-1:i1,:);
    else
        q(i1:-1:i2,:)=p(i2:i1,:);
    end
end