%% Create Random solution that meet resource allocation constraints
function f=createRandomSolution(model)
C=model.C;
D=model.D;
N_RB=model.N_RB;

x_ik=zeros(C,N_RB);
x_jk=zeros(D,N_RB);
Rc=randperm(C);
Rd=randperm(D);

for i=1:C
x_ik(Rc(i), i)=1;
end
for j=1:D
x_jk(Rd(j), j)=1;
end

f1= x_ik(:,1:N_RB);
f2=x_jk(:,1:N_RB);
f=[f1;f2];