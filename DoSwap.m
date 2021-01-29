function q=DoSwap(p,i1,i2)

q=p;
q(i1,:)=p(i2,:);
q(i2,:)=p(i1,:);
end
