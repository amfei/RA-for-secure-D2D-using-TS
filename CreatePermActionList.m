function ActionList=CreatePermActionList(C) % n number of element that we are going to calculate its permutation

nSwapC=((C*(C-1)/2)+1)^2;
% nReversionC=C*(C-1)/2;
% nInsertionC=C^2;
% nSwapD=D*(D-1)/2;
% nReversionD=D*(D-1)/2;
% nInsertionD=D^2;

%nAction=nSwapC+nReversionC+nInsertionC+nSwapD+nReversionD+nInsertionD;
%nAction=nSwapC+nInsertionC+nSwapD+nInsertionD;
nAction=nSwapC
ActionList=cell(nAction,1);

c=0;

% Add SWAP

for i=1:C-1
    for j=i+1:C
            c=c+1;
            ActionList{c}=[1 i j]
    end
end
% for i=C+1:C+D-1
%     for j=i+1:C+D
%         c=c+1;
%         ActionList{c}=[1 i j];
%     end
% end





% % Add Insertion
% for i=1:C
%     for j=1:C
%         if i<j
%         if abs(i-j)>1  % faghat kenar ham nabad bashad
%             c=c+1;
%             ActionList{c}=[2 i j];
%         end
%         else
%            if abs(i-j)>2  % faghat kenar ham nabad bashad
%                c=c+1;
%               ActionList{c}=[2 i j];
%           end
%         end
%     end
% end
% 
% for i=C+1:C+D
%     for j=C+1:C+D
%         if i<j
%         if abs(i-j)>1  % faghat kenar ham nabad bashad
%             c=c+1;
%             ActionList{c}=[2 i j];
%         end
%                 else
%                     if abs(i-j)>2  % faghat kenar ham nabad bashad
%                         c=c+1;
%                         ActionList{c}=[2 i j];
%                     end
%         end
%     end
% end
% 
% 
% % Add REVERSION
% for i=1:C-1
%     for j=i+1:C
%         if abs(i-j)>2  % kenar ham ya yek adad ekhtelaf dashte bashd reversen yeki ast va niyazi nist
%             c=c+1;
%             ActionList{c}=[3 i j];
%         end
%     end
% end
% 
% for i=C+1:C+D-1
%     for j=i+1:C+D
%         if abs(i-j)>2  % kenar ham ya yek adad ekhtelaf dashte bashd reversen yeki ast va niyazi nist
%             c=c+1;
%             ActionList{c}=[3 i j];
%         end
%     end
% end
% 
% 

ActionList=ActionList(1:c);

end