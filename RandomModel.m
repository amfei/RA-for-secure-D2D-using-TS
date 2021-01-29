

function model = RandomModel(C,D,N_RB,Pc,Pd,CR,dd,Rcmin,Rdmin,sigma,alpha,beta)
%cell R and BS
xBS = 0 ;
yBS= 0;
%circle
t = linspace(0, 2*pi, 100);
xr = CR*cos(t) + xBS;
yr = CR*sin(t) + yBS;
% Eav in cell
E=2*pi*rand(1,1);
rE  = CR*rand(1,1);
xE =rE.*cos(E) ;
yE = rE.*sin(E);
%CUs in circle cell
centerC = [0 ,0];
anglec = 2*pi*rand(1,C);
rc = CR*sqrt(rand(1,C));
xc = rc.*cos(anglec)+ centerC(1) ;
yc = rc.*sin(anglec)+ centerC(2) ;
%DU transmiter in R
centerD = [0 ,0];
angledt = 2*pi*rand(1,D);
rdt = CR*sqrt(rand(1,D));
xdt = rdt.*cos(angledt)+ centerD(1);
ydt = rdt.*sin(angledt)+ centerD(2);

for dt=1:D
    centerDr=[xdt(dt) ydt(dt)];
    angledr= 2*pi*rand;
    xdr(dt) = dd(dt).*cos(angledr)+ centerDr(1);
    ydr(dt) = dd(dt).*sin(angledr)+ centerDr(2);

end


%% plot
% plot(xBS,yBS,'gh','LineWidth',2,'MarkerSize',19,'MarkerEdgeColor','r','MarkerFaceColor',[1 0 1])
% text((xBS),(yBS-30),'BS');
% hold on
% plot(xr,yr);
% hold on
% plot(xE,yE,'r^','LineWidth',2,'MarkerSize',10,'MarkerEdgeColor','r','MarkerFaceColor',[1 1 0]);
% text((xE),(yE-20),'Eve');
% hold on
% plot(xc,yc,'gs','LineWidth',2,'MarkerSize',8,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
% hold on
% plot(xdt,ydt,'*r')
% for dt=1:D
%     text((xdt(dt)-6),(ydt(dt)-6),num2str(dt));
% end
% hold on
% for dt=1:D
% plot(xdr,ydr,'*b')
% hold on
% plot([xdr(dt) xdt(dt)],[ydr(dt) ydt(dt)],'or-')
% end
% for cc=1:C
%     text((xc(cc)+8),(yc(cc)+8),num2str(cc));
% end
% % 

%% Distance
% CUs  to BS
diB=zeros(1,C);
for b=1:C
    diB(b)=sqrt((xc(b)).^2 + (yc(b)).^2);
end

%CUs to Eav
die= zeros(1,C);
for v=1:C
    die(v)=  sqrt((xc(v)-xE).^2 + (yc(v)-yE).^2); % 1xN
end


%  Dt to BS
djtB=zeros(1,D);
for g=1:D
    djtB(g)=sqrt(  (xdt(g))^2+(ydt(g)^2 ));
end

%  Dt to evesdropper
djt_e= zeros(1,D);
for pp=1:D
    djt_e(pp)=  sqrt(   (xdt(pp)-xE).^2 + (ydt(pp)-yE).^2); % 1xN
end

% CUs  to Dr
dijr=zeros(C,D);
for z=1:C
    for w=1:D
        dijr(z,w)=sqrt(   (xc(z)-xdr(w)).^2 + (yc(z)-ydr(w)).^2);
    end
end


%% pathloss iB
%PL_iB=40*(1-4*0.02)*log10(diB/1000)-18*log10(20+21*log10(1000*fc))+80;% h_BS=20 meter
%PL_iB2=36.7.*log10(diB)+22.7+26.*log10(fc);
PL_iB=35*log10(diB)+31.5;
PL_iB_lin=10.^(-PL_iB/10);% convert to linear
% Rayleigh channel coeffiect
%Rc_iB = ((abs(1/sqrt(2)*(randn(C,N_RB) + 1j*randn(C,N_RB)))).^2 );
h_iB = sqrt(0.5)*randn(C,N_RB) +sqrt(0.5)* 1j*randn(C,N_RB);
% iB channek gain
g_iB=zeros(C,N_RB);
for iB=1:1:C
    %g_iB(iB,:)=nG*PL_iB_lin(iB)*Rc_iB(iB,:);
    g_iB(iB,:)=(PL_iB_lin(iB)*abs(h_iB(iB,:)));
end

%% pathloss ie
%PL_ie=40*(1-4*0.02)*log10(die/1000)-18*log10(20+21*log10(1000*fc))+80;
PL_ie=35*log10(die)+31.5;
%PL_ie=36.7.*log10(die)+22.7+26.*log10(fc);
PL_ie_lin=10.^(-PL_ie/10);% convert to linear
% Rayleigh channel coeffiect
%Rc_ie = ((abs(1/sqrt(2)*(randn(C,N_RB) + 1j*randn(C,N_RB)))).^2 );
h_ie = sqrt(0.5)*randn(C,N_RB) +sqrt(0.5)* 1j*randn(C,N_RB);
% ie channek gain
g_ie=zeros(C,N_RB);
for i=1:1:C
    %g_ie(i,:)=nG*PL_ie_lin(i)*Rc_ie(i,:);
    g_ie(i,:)=(PL_ie_lin(i)*abs(h_ie(i,:)));
end

%% pathloss jt-B
%PL_jtB=36.7.*log10(djtB)+22.7+26.*log10(fc);
%PL_jtB=40*(1-4*0.02)*log10(djtB/1000)-18*log10(20+21*log10(1000*fc))+80;
PL_jtB=35*log10(djtB)+31.5;
PL_jtB_lin=10.^(-PL_jtB/10);% convert to linear

% Rayleigh channel coeffiect
%Rc_jtB = ((abs(1/sqrt(2)*(randn(D,N_RB) + 1j*randn(D,N_RB)))).^2 );
h_jtB = sqrt(0.5)*randn(D,N_RB) +sqrt(0.5)* 1j*randn(D,N_RB);
% ie channek gain
g_jtB=zeros(D,N_RB);
for jt=1:1:D
    g_jtB(jt,:)=(PL_jtB_lin(jt)*abs(h_jtB(jt,:)));
end

%% pathloss jt_e
%PL_jt_e=40*(1-4*0.02)*log10(djt_e/1000)-18*log10(20+21*log10(1000*fc))+80;
%PL_jt_e=36.7.*log10(djt_e)+22.7+26.*log10(fc);
PL_jt_e=35*log10(djt_e)+31.5;
PL_jt_e_lin=10.^(-PL_jt_e/10);% convert to linear
% Rayleigh channel coeffiect
%Rc_jt_e = ((abs(1/sqrt(2)*(randn(D,N_RB) + 1j*randn(D,N_RB)))).^2 );
h_jt_e = sqrt(0.5)*randn(D,N_RB) + sqrt(0.5)*1j*randn(D,N_RB);
% ie channek gain
g_jt_e=zeros(D,N_RB);
for jte=1:1:D
    g_jt_e(jte,:)=(PL_jt_e_lin(jte)*abs(h_jt_e(jte,:)));
end


%% pathloss D2D
%PL_D2D=36.7.*log10(dd)+22.7+26.*log10(fc);
PL_D2D=40*log10(dd)+31.5;
PL_D2D_lin=10.^(-PL_D2D/10);% convert to linear
%PL_D2D_lin=dd^-4
% Rayleigh channel coeffiect
%Rc_D2D = ((abs(1/sqrt(2)*(randn(D,N_RB) + 1j*randn(D,N_RB)))).^2 );
h_D2D = sqrt(0.5)*randn(D,N_RB) + sqrt(0.5)*1j*randn(D,N_RB);

% ie channek gain
g_J=zeros(D,N_RB);
for J=1:D
    g_J(J,:)=(PL_D2D_lin(D)*abs(h_D2D(J,:)));
end

%% pathloss i-jr
%PL_dijr=40*(1-4*0.02)*log10(dijr/1000)-18*log10(20+21*log10(1000*fc))+80;
%PL_dijr=36.7.*log10(dijr)+22.7+26.*log10(fc);
PL_dijr=35*log10(dijr)+31.5;
PL_dijr_lin=10.^(-PL_dijr/10);% convert to linear
% Rayleigh channel coeffiect
%Rc_dijr = ((abs(1/sqrt(2)*(randn(D,N_RB) + 1j*randn(D,N_RB)))).^2 );
% ijr channek gain
g_ijr=zeros(C,D,N_RB);

for ii=1:C
    for jr=1:D
        %Rc_dijr = ((abs(1/sqrt(2)*(randn(1,N_RB) + 1j*randn(1,N_RB)))).^2);
        h_dijr = sqrt(0.5)*randn(1,N_RB) +sqrt(0.5)* 1j*randn(1,N_RB);
        g_ijr(ii,jr,:)=(PL_dijr_lin(ii,jr).*abs(h_dijr)) ; %C*D*N_RB
    end
end

%%create random Assignment for initial solution of heuristic method
x_ik=zeros(C,N_RB);
x_jk=zeros(D,N_RB);
Rc=randperm(C);
Rd=randperm(D);

for xi=1:C
x_ik(Rc(xi), xi)=1;
end
for xj=1:D
x_jk(Rd(xj), xj)=1;
end
model.x_ik=x_ik;
model.x_jk=x_jk;
model.g_ijr=g_ijr;
model.g_iB=g_iB;
model.g_J=g_J;
model.g_jt_e=g_jt_e;
model.g_jtB=g_jtB;
model.g_ie =g_ie;
model.C=C;
model.D=D;
model.N_RB=N_RB;
model.Pc=Pc;
model.Pd=Pd;
model.sigma=sigma;
model.CR=CR;
model.dd=dd;
model.Rcmin=Rcmin;
model.Rdmin=Rdmin;
model.alpha=alpha;
model.beta=beta;
end