  function CostT=getSecCap(f,model)
% Secracy Capacity
%sigma=7.1659e-16;% 180khz*4e-21(-174dBm)
alpha=model.alpha;
beta=model.beta;
Rcmin=model.Rcmin;
Rdmin=model.Rdmin;
sigma=model.sigma;
C=model.C;
D=model.D;
N_RB=model.N_RB;
Pc=model.Pc;
Pd=model.Pd;
f1=f(1:C,:);
f2=f(C+1:C+D,:);
costc=zeros(1,N_RB);
  for k=1:N_RB
      for i=1:C
        for j=1:D
            if (f1(i,k) ~= 0) && (f2(j,k)~= 0)  %&& (log2(1+(Pc*g_iB)/(Pd*g_jtB+BOLTZ*293*BW/N_RB)))
            %g_ijr=model.g_ijr(i,j,k);
             g_iB=model.g_iB(i,k);
            %g_J=model.g_J(j,k);
            g_jt_e=model.g_jt_e(j,k);
            g_jtB=model.g_jtB(j,k);
            g_ie=model.g_ie(i,k);
           
           
            %costc(k)=max(0, log2(1+(Pc*g_iB)/(Pd*g_jtB+BOLTZ*293*BW/N_RB))-log2(1+(Pc*g_ie)  /(Pd*g_jt_e+BOLTZ*293*BW/N_RB)))
            ViolationCU=max(1-(log2(1+(Pc*g_iB)/(Pd*g_jtB+sigma)))/Rcmin,0);%(1-g/g0)
            %costc(k)=max(0, log2(1+(Pc*g_iB)/(Pd*g_jtB+BOLTZ*293*BW/N_RB))-log2(1+(Pc*g_ie)  /(Pd*g_jt_e+BOLTZ*293*BW/N_RB)))-alpha*ViolationCU;
           
             costc(k)= (max(log2(1+(Pc*g_iB)/(Pd*g_jtB+sigma))-log2(1+(Pc*g_ie)/(Pd*g_jt_e+sigma)),0)-alpha*ViolationCU);
       
            end
        end
      end  
  end
   CostTc=sum(costc);
   
  costd=zeros(1,N_RB);
for kk=1:N_RB
    for jj=1:D
        for ii=1:C
            if (f1(ii,kk) ~= 0) && (f2(jj,kk)~= 0)  %(log2(1+(Pd*g_J )/(Pc*g_ijr+BOLTZ*293*BW/N_RB)))
                g_ijr=model.g_ijr(ii,jj,kk);
                %g_iB=model.g_iB(i,k);
                g_J=model.g_J(jj,kk);
                g_jt_e=model.g_jt_e(jj,kk);
                %g_jtB=model.g_jtB(j,k);
                g_ie=model.g_ie(ii,kk);
                
                
                ViolationD2D=max(1-(log2(1+(Pd*g_J )/(Pc*g_ijr+sigma)))/Rdmin,0);
               %costd(kk)= max(0, log2(1+(Pd*g_J )/(Pc*g_ijr+BOLTZ*293*BW/N_RB))-log2(1+(Pd*g_jt_e)/(Pc*g_ie+BOLTZ*293*BW/N_RB)))
               % costd(kk)= max(0, log2(1+(Pd*g_J )/(Pc*g_ijr+BOLTZ*293*BW/N_RB))-log2(1+(Pd*g_jt_e)/(Pc*g_ie+BOLTZ*293*BW/N_RB)))-beta*ViolationD2D;
               costd(kk)=(max (log2(1+(Pd*g_J )/(Pc*g_ijr+sigma))-log2(1+(Pd*g_jt_e)/(Pc*g_ie+sigma)),0)-beta*ViolationD2D);
%ConstT.ViolationD2D=ViolationD2D;
            end
        end
    end
    
end
CostTd=sum(costd);
CostT=CostTc+CostTd;

  
  
 
  