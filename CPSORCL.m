function [Epistasis, MI, time] = CPSOAST(FileName,ParticleNumber,IterationNumber,Order,ReturnNumber)

%%
% [Epistasis, MI, time] = CPSOAST('dat100', 100, 50, 3, 5)
% [Epistasis, MI, time] = CPSOAST('Model1_Set1', 100, 50, 2, 5)
%%
tic;
filename=[FileName,'.mat'];
temp = load(filename);
pts = temp.pts;
class = temp.class;

%% Unimportant Parameters 不重要的参数
snpdata = size(pts,2);
% maximum particle position  最大粒子位置
xmax = snpdata;
% minimum particle position 最小粒子位置
xmin = 1;
% maximum particle velocity 最大粒子速度
vmax = xmax;
% minimum particle velocity 最小粒子速度
vmin = - xmax;
% the acceleration coefficients controlling 加速度系数控制
c1 = 2.05;
c2 = 2.05;
% inertia weight惯性权重
W=0.65;
%初始化SNP权重矩阵
SNPWeight = zeros(IterationNumber,snpdata);

%% Initialize Parameters 初始参数
% initialize particle position 初始化粒子位置
par_x = cell(1,ParticleNumber);
par_mut = zeros(1,ParticleNumber);
for i = 1:1:ParticleNumber
    par_x{i} = randperm(xmax,Order);
    par_mut(i) = MutualInformation(pts, class, par_x{i});
end

% initialize particle velocity  初始化粒子速度
par_v = cell(1, ParticleNumber);
for i = 1:1:ParticleNumber 
    par_v{i} = randi([vmin, vmax], 1, Order);
end

 [~, index] = sort(par_mut, 'DESCEND');
 Epistasis=par_x(index(1,1:ReturnNumber));

disp('Waiting...');

%%
for t = 1:1:IterationNumber
    %% 随机拓扑更新策略
    m=t+1;
    if m>=sqrt(IterationNumber)
        m=round(sqrt(IterationNumber));
    end
    %%
    A=zeros(ParticleNumber,m);
    count=zeros(ParticleNumber,1);
    
    %%
    par_index_x=cell(ParticleNumber,m);
    par_index_v=cell(ParticleNumber,m);
    par_index_mut=zeros(ParticleNumber,m);
    
     Achieve_x = cell(1,ParticleNumber);
     Achieve_v = cell(1,ParticleNumber);
     Leaderbest_x=cell(ParticleNumber,1);
     Leaderbest_v=cell(ParticleNumber,1);
     Leaderworst_x=cell(ParticleNumber,1);
     Leaderworst_v=cell(ParticleNumber,1);
     
    %% 判断每个粒子是否存在两个leader
    for i=1:1:ParticleNumber
        A(i,:)=randperm(ParticleNumber,m);
        while sum(ismember(A(i,:),i))
            A(i,:)=randperm(ParticleNumber,m);
        end
    
        for j=1:1:m
            if par_mut(A(i,j))>par_mut(i)
                count(i)=count(i)+1;
                par_index_x{i,j}=par_x{A(i,j)};
                par_index_v{i,j}=par_v{A(i,j)};
                par_index_mut(i,j)=par_mut(A(i,j));
            end
        end
 
%% 若存在 保留最优以及最差解的位置和速度
        if count(i)>=2
            [~, loc] = sort(par_index_mut(i,:), 'DESCEND');
            Leaderbest_x{i}=par_index_x{i,loc(1)};            
            Leaderbest_v{i}=par_index_v{i,loc(1)};
            Leaderworst_x{i}=par_index_x{i,loc(count(i))};
            Leaderworst_v{i}=par_index_v{i,loc(count(i))};
            
            %% update particle velocity  更新粒子速度
            for k = 1:1:Order
                Achieve_v{i}(k) = W *par_v{i}(k) + c1 * rand(1) * (Leaderbest_x{i}(k) - par_x{i}(k))+ c2 * rand(1) * (Leaderworst_x{i}(k)- par_x{i}(k));
                Achieve_v{i}(k) = round(Achieve_v{i}(k));
            end
            
           %% update particle position 更新粒子位置
           % 放在函数最大最小值等三维连续不断平滑坡面的效果一定不错。但是，放在SNP互作上，尤其是随机SNP之间相互独立的话，
           % 这就是一个纯粹随机的粒子更新。
           % 这里要注意，其实在真实SNP数据之间的话，效果会比目前的仿真数据上效果要好。原因是仿真数据通常不仿真SNP之间的生物
           % 生物特性，比如LD或者单体型快等。真实数据上有这些属性，使得其更像连续不断平滑坡面。
           % 所有如果在仿真数据上仿真上生物特性的话，效果也会好很好。
                for k = 1:1:Order
                    Achieve_x{i}(k) = par_x{i}(k) + Achieve_v{i}(k);
                    if Achieve_x{i}(k) < xmin || Achieve_x{i}(k) > xmax
                        Achieve_x{i}(k) = randperm(xmax,1);
                    end 
                end
                
                if numel(unique(Achieve_x{i}))~=Order
                    Achieve_x{i} = randperm(xmax,Order);
                end              
                
                Achieve_mut = MutualInformation(pts, class,Achieve_x{i});
                if Achieve_mut>par_mut(i)
                    par_v{i}=Achieve_v{i};
                    par_x{i}=Achieve_x{i};
                    par_mut(i)= Achieve_mut;
                end  
         end
    end   
    
    %% 翻转策略
    %% 记录每次迭代SNP权重 
    % 这是取当前所有含该SNP的MI的最大值为其当前权重，二选一
    for i=1:1:ParticleNumber
        for j=1:1:Order
            if SNPWeight(t,par_x{i}(j))<par_mut(i)
                SNPWeight(t,par_x{i}(j))=par_mut(i);
            end
        end
    end
    
    % 这是取当前所有含该SNP的MI的累加值为其当前权重，二选一
    % 考虑选择这个的原因是，权重变化速度显著加快。既然前面更新粒子本质上趋向于随机。那么这里的更新粒子，就应该更趋向于
    % 向优。对于多次取到的SNP，或者SNP组合MI值显著高的SNP，通过累加其权重显著高。那么其翻转的可能性就变小。在粒子里
    % 不变的可能性变大。
    % for i=1:1:ParticleNumber
    %     for j=1:1:Order
    %         SNPWeight(t,par_x{i}(j))=SNPWeight(t,par_x{i}(j))+par_mut(i);
    %     end
    % end
    
    %% 归一化SNP权重
    SUM_SNPWeight=sum(SNPWeight);
    
    minValue = min(SUM_SNPWeight);
    maxValue = max(SUM_SNPWeight);
    SNP_W_Normal = (SUM_SNPWeight - minValue) / (maxValue - minValue);
 
    %% 计算每个SNP的翻转概率并更新粒子位置
    Probability_SNP=(1-SNP_W_Normal)*(1-((t-1)/IterationNumber));
    
    for i = 1:1:snpdata
        if rand < Probability_SNP(1,i)
            for j=1:1:ParticleNumber
                for k=1:1:Order
                    if par_x{j}(k)==i
                        par_x{j}(k) = snpdata+1-par_x{j}(k) ;
                    end
                end               
            end
        end
    end
    
    for i = 1:1:ParticleNumber
        if numel(unique(par_x{i}))~=Order
            par_x{i} = randperm(xmax,Order);
        end
        par_mut(i) = MutualInformation(pts, class, par_x{i});
    end

    %%
    [~, index] = sort(par_mut, 'DESCEND');
    Epistasis_temp=par_x(index(1,1:ReturnNumber));
    
    [Epistasis, ~] = Integration(pts, class, Epistasis,  Epistasis_temp, ReturnNumber);
    
end    

%% 遍历
avervalue = mean(SNP_W_Normal);
stdvalue = std(SNP_W_Normal);
SNPindex = 1:snpdata;
SNPs= SNPindex(SNP_W_Normal>(avervalue+stdvalue));

if size(SNPs,2)>100
    SNPWeights= SNP_W_Normal(SNP_W_Normal>(avervalue+stdvalue));
    [~,index]=sort(SNPWeights,'descend');
   SNPs=SNPs(index(1:100));
end

%%
SNPs2=[];
for i=1:1:ReturnNumber
    SNPs2=union(SNPs2, Epistasis{i});
end
SNPs=union(SNPs,SNPs2');

%%
Candidate = nchoosek(SNPs, Order);
CandidateNum = size(Candidate,1);

Candidate_Epi = cell(1,CandidateNum);
Candidate_mut = zeros(1,CandidateNum);

for i=1:1:CandidateNum
    Candidate_Epi{i} = Candidate(i,:);
    Candidate_mut(i) = MutualInformation(pts, class, Candidate_Epi{i});
end

[~, index] = sort(Candidate_mut, 'DESCEND');
Epistasis_temp=Candidate_Epi(index(1,1:min(CandidateNum,ReturnNumber)));

[Epistasis, MI] = Integration(pts, class, Epistasis,  Epistasis_temp, ReturnNumber);

%%
time=toc;

SaveName=[FileName,'_',num2str(ParticleNumber),'_',num2str(IterationNumber),'_',num2str(Order),'_',num2str(ReturnNumber),'.mat'];

save(SaveName,'Epistasis','MI','time');


function MI=MutualInformation(pts,class,factor)
% Compute Mutual Information between class label and a SNP set.计算类标签和SNP集之间的互信息
% MI(X;C)=H(X)+H(C)-H(X,C), where C is the class label and X is a SNP set.
data=pts(:,factor);
MI=Func_MutualInfo(double(data)',double(class));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MI=Func_MutualInfo(data,labels)
% Jan.19 2007 By Guoqiang Yu
%--- INPUT
% data:  features x samples; note that the 样本特征
% value of zero indicates a missing data 值为0表示缺失数据
% labels: the label matrix of each sample;  每个样本的标签矩阵
%--- OUTPUT
% MI: Information provided by the features 特性提供的信息
%%
data=data-min(data(:));          %minimum data is 0
labels=labels-min(labels(:))+1; %minimum label is 1

Num_Label=max(labels(:));
Num_DataType=max(data(:))+1;
[Num_Feature,Num_Sample]=size(data);

%Entropy of the random variable Label 随机变量Label的熵
H_Label=hist(labels,Num_Label);
P_Label=H_Label/Num_Sample;
Entropy_Label=-P_Label*log2(P_Label'+eps);

%Special dealing with the case of small Num_Feature 特殊处理小Num_Feature的情况
if Num_Feature<9
    ZZ=Num_DataType.^(Num_Feature-1:-1:0)';
    Hist_Label=zeros(Num_DataType^Num_Feature,Num_Label);%  Hist_Label is p(c,f)
    tempIndex=ZZ'*data+1;
    for j=1:Num_Sample
        Hist_Label(tempIndex(j),labels(j))=Hist_Label(tempIndex(j),labels(j))+1; % calculate p(c,f) 
    end
    
    sumHist=sum(Hist_Label,2);   %calculate p(f)
    repHist=repmat(sumHist,1,Num_Label); 
    pHist_Label=Hist_Label./(repHist+eps);%p(c/f)=p(c,f)/p(f).
    InfoIncre=-sum((log2(pHist_Label+eps).*pHist_Label).*(repHist));
    MI=Entropy_Label-sum(InfoIncre)/Num_Sample;
    return;
end

%Larger Feature Number follows the following procedure （大特征）遵循以下
mm=1;
Hist_Label=zeros(Num_Label,Num_Sample);
Hist_Label(labels(1,1),mm)=1;
Hist_SNP=zeros(Num_Feature,Num_Sample);
Hist_SNP(:,mm)=data(:,1);

for j=2:Num_Sample
    tempData=data(:,j);
    Index=0;
    for k=1:mm
        if isequal(Hist_SNP(:,k),tempData)
            Index=k;
            break;
        end
    end
    if Index==0
        mm=mm+1;
        Hist_SNP(:,mm)=tempData;
        Hist_Label(labels(j,1),mm)=1;
    else
        Hist_Label(labels(j,1),Index)=Hist_Label(labels(j,1),Index)+1;
    end
end

M1=mm;
InfoIncre=0;
for s=1:M1
    tempNum=sum(Hist_Label(:,s));
    P=Hist_Label(:,s)/tempNum;
    InfoIncre=InfoIncre-P'*log2(P+eps)*tempNum;
end
MI=Entropy_Label-InfoIncre/Num_Sample;


function [Epistasis, MI] = Integration(pts, class, Epistasis,  Epistasis_temp, ReturnNumber)

%%
for i=1:1:size(Epistasis,2)
    Epistasis{i}=sort(Epistasis{i});
end

for i=1:1:size(Epistasis_temp,2)
    Epistasis_temp{i}=sort(Epistasis_temp{i});
end

% 将整数向量转换为字符串
A_str = cellfun(@num2str, Epistasis, 'UniformOutput', false);
B_str = cellfun(@num2str, Epistasis_temp, 'UniformOutput', false);

% 合并并去除重复元素
mergedCell = unique([A_str, B_str], 'stable');

number=size(mergedCell,2);
Epistasis_all=cell(1,number);
MI_all=zeros(1,number);

for i=1:1:number
    Epistasis_all{i}=str2num(mergedCell{i}); %#ok<*ST2NM>
    MI_all(i)= MutualInformation(pts, class, Epistasis_all{i});
end

 [Sort_par_mut, index] = sort(MI_all, 'DESCEND');
 MI=Sort_par_mut(1,1:ReturnNumber);
 Epistasis=Epistasis_all(index(1,1:ReturnNumber));

