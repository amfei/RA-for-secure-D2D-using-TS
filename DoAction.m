function q=DoAction(p,a) % a has 3 elements

    switch a(1) % first element of a (action) shows the type of action 
        case 1
            % Swap
            q=DoSwap(p,a(2),a(3)); 
            
        case 2
           q=DoInsertion(p,a(2),a(3));
       case 3           
             q=DoReversion(p,a(2),a(3));
    end

end