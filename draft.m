clc
clear
close all
CR = 500;   %CellRadius

BW=20e6; %BW
fc=2.3;
Pc=0.200 ;%
Pd=0.02 ;%
L=1;  %mont carlo loop
nG=1.3955e+15 ;%  
BOLTZ=1.3806488e-23;

C=6;
N_RB=C;


D=C;
    tic
    dd = 20*ones(1,D); % 20 meter d2d distance
    MaxIt=(C+D);
    BestCostLoops=zeros(MaxIt,L);
    for loop=1:L
        
        model = RandomModel(C,D,fc,N_RB,Pc,Pd,CR,dd,BW,nG,BOLTZ);
        CostFunction=@(f) getSecCap(f,model);
        ActionList=CreatePermActionList(C,D) ;    % Action List
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
        sol.Position=createRandomSolution(model)  ;  %S.  Initial_sol
        
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
                if i==bestnewsol.ActionIndex;
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
            %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
        end
        
        BestCostLoops(:,loop)=BestCost;
    end
    [MaxBestCostLoops, I]=max(max(BestCostLoops))
    [MinBestCostLoops ,  J]=min(max(BestCostLoops));
    AveBestCostLoops=max((sum(BestCostLoops,2))./L);

    

    
    toc





