function q=DoInsertion(p,i1,i2)% satr i1 ra be paeen satr i2 ezaf kon
if i1<i2
    
    q=[p(1:i1-1,:); p(i1+1:i2,:); p(i1,:); p(i2+1:end,:)];
else
    q=[p(1:i2,:) ;p(i1,:); p(i2+1:i1-1,:); p(i1+1:end,:)];
end
      