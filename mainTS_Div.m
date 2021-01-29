clc
clear
close all
CR = 500;   %CellRadius
L=100;  %mont carlo loop
% BW=20e6; %BW
fc=2.3;
Pc=0.200 ;%
Pd=0.02 ;%

% nG=1.3955e+15 ;%  
% BOLTZ=1.3806488e-23;
sigma=7.1659e-16;% 180khz*4e-21(-174dBm)
alpha=1;
beta=10;
Rdmin=4;
Rcmin=1;


pc=0.200 ;%
gamma = 0.5;    % discount factor  % TODO : we need learning rate schedule
alpha = 0.5;    % learning rate    % TODO : we need exploration rate schedule
epsilon = 0.6;  % exploration probability (1-epsilon = exploit / epsilon = explore)
pathloss_parameter = 3.5;
Rc=400; % level of distance between Cu and d2d for assignment
T0=1; % min sinr guaranteeing QoS cus


C=3;
K=C;   
D=3;
tic
    MaxIt=(C+D);
    BestCostLoops=zeros(MaxIt,L);
    for loop=1:L
        dd = randi(20,1,D);
         model = RandomModel(C,D,K,pc,Pd,CR,dd,Rcmin,Rdmin,sigma,alpha,beta);
         
        CostFunction=@(f) getSecCap(f,pc_opt,pd_opt,model);
        ActionList=CreatePermActionList(C,K,D) ;    % Action List
        nAction=numel(ActionList);              % Number of Actions
        %% Tabu Search Parameters
        % Maximum Number of Iterations
        TL=round(0.5*nAction); % Tabu Length

        %% Initialization
        % Create Empty Individual Structure
        empty_individual.Position=[];
        empty_individual.Cost=[];  
        % Create Initial Solution
        sol=empty_individual;
        sol.Position=createRandomSolution(model) ;   %S.  Initial_sol
        
        sol.Cost=CostFunction(sol.Position);
        % Initialize Best Solution Ever Found
        BestSol=sol;                                             %S**=S.  Global_sol <-- Initial_sol
        
        % Array to Hold Best Costs
        BestCost=zeros(MaxIt,1);  % f(S*)=[0 0 0 ..]'
        
        % Initialize Action Tabu Counters
        TC=zeros(nAction,1);
        TCin=zeros(nAction,1);
        
        %% Tabu Search Main Loop
        for it=1:MaxIt
            if it> 0.7*MaxIt % Diversification condition
                %divsol.Position=DoReversion(sol.Position,1,C+D); % S0d
                divsol.Position=createRandomSolution(model); % S0d
                divsol.Cost=CostFunction(divsol.Position) ;    %f(S0d)
                sol=divsol;
            end
            
            bestnewsol.Cost=sol.Cost    ;   %f(S*)=f(S.)     Optimal_local_sol <--initial_sol
            % Apply Actions
            for i=1:nAction
                if TC(i)==0
                    newsol.Position=DoAction(sol.Position,ActionList{i}); % S'  local_sol
                    newsol.Cost=CostFunction(newsol.Position) ;    %f(S')
                    newsol.ActionIndex=i;   %
                    if newsol.Cost>=bestnewsol.Cost; % f(S')>f(S*)    local_sol
                        bestnewsol=newsol; %(f&S)*=(f&S)'       Optimal_local_sol <-- local_sol
                    end
                end
            end
            
            % Update Local Solution
            sol=bestnewsol;% (f&S.)=(f&S*)        initial_sol<< Optimal_local_sol
            % Update Tabu List
            for i=1:nAction
                if i==bestnewsol.ActionIndex
                    TC(i)=TL   ;       % Add To Tabu List
                else
                    TC(i)=max(TC(i)-1,0);  % Reduce Tabu Counter
                end
            end
            
            % Update Best Solution Ever Found
            if sol.Cost>=BestSol.Cost                 % initial_sol > Global_sol
                BestSol=sol;                   %Global sol <- initial_sol
            end
            
            % Save Best Cost Ever Found
            BestCost(it)=BestSol.Cost;

            %% Show Iteration Information
            disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
        end
        
        BestCostLoops(:,loop)=BestCost;
    end
    [MaxBestCostLoops, I]=max(max(BestCostLoops));
    [MinBestCostLoops ,  J]=min(max(BestCostLoops));
    AveBestCostLoops=max((sum(BestCostLoops,2))./L)
    toc
    
     plot(1:MaxIt,BestCostLoops(:,I),1:MaxIt,BestCostLoops(:,J),1:MaxIt,AveBestCostLoops,'LineWidth',1.5);
     hold on
     grid on
     xlabel('Number of Iteration');
     ylabel('Total Secracy Capacity(bit/sec/Hz)');
    
    

%     Secracy(loop)= getSecCap2(S,model); % we have differnt cost in each iteration because the chanel are random 
% 
% MaxSec(N_RB)=max(Secracy);
% AveSec(N_RB)=sum(Secracy)/l;
% MinSec(N_RB)=min(Secracy);
% 
% MaxSec(MaxSec==0)=[];
% AveSec(AveSec==0)=[]
% MinSec(MinSec==0)=[]
    
    

% legend('Max. Proposed TS','Ave. Proposed TS','Min. Proposed TS','FontSize',12,'FontWeight','bold','Location','northwest')

% MaxBestCostLoops(MaxBestCostLoops==0) = [];
% AveBestCostLoops(AveBestCostLoops==0) = []
% MinBestCostLoops(MinBestCostLoops==0) = [];
%  figure;
% plot(10:10:N_RB,MaxBestCostLoops,'r-o',10:10:N_RB,AveBestCostLoops,'b-*',10:10:N_RB,MinBestCostLoops,'k-s','LineWidth',1.5,'MarkerSize',8);
%  legend('Max. Proposed TS','Ave. Proposed TS','Min. Proposed TS','FontSize',12,'FontWeight','bold','Location','northwest')
% grid on
% xlabel('Number of D2D pairs');
% ylabel('Total Secracy Capacity(bit/sec/Hz)');




