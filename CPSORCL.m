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

%% Unimportant Parameters ����Ҫ�Ĳ���
snpdata = size(pts,2);
% maximum particle position  �������λ��
xmax = snpdata;
% minimum particle position ��С����λ��
xmin = 1;
% maximum particle velocity ��������ٶ�
vmax = xmax;
% minimum particle velocity ��С�����ٶ�
vmin = - xmax;
% the acceleration coefficients controlling ���ٶ�ϵ������
c1 = 2.05;
c2 = 2.05;
% inertia weight����Ȩ��
W=0.65;
%��ʼ��SNPȨ�ؾ���
SNPWeight = zeros(IterationNumber,snpdata);

%% Initialize Parameters ��ʼ����
% initialize particle position ��ʼ������λ��
par_x = cell(1,ParticleNumber);
par_mut = zeros(1,ParticleNumber);
for i = 1:1:ParticleNumber
    par_x{i} = randperm(xmax,Order);
    par_mut(i) = MutualInformation(pts, class, par_x{i});
end

% initialize particle velocity  ��ʼ�������ٶ�
par_v = cell(1, ParticleNumber);
for i = 1:1:ParticleNumber 
    par_v{i} = randi([vmin, vmax], 1, Order);
end

 [~, index] = sort(par_mut, 'DESCEND');
 Epistasis=par_x(index(1,1:ReturnNumber));

disp('Waiting...');

%%
for t = 1:1:IterationNumber
    %% ������˸��²���
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
     
    %% �ж�ÿ�������Ƿ��������leader
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
 
%% ������ ���������Լ������λ�ú��ٶ�
        if count(i)>=2
            [~, loc] = sort(par_index_mut(i,:), 'DESCEND');
            Leaderbest_x{i}=par_index_x{i,loc(1)};            
            Leaderbest_v{i}=par_index_v{i,loc(1)};
            Leaderworst_x{i}=par_index_x{i,loc(count(i))};
            Leaderworst_v{i}=par_index_v{i,loc(count(i))};
            
            %% update particle velocity  ���������ٶ�
            for k = 1:1:Order
                Achieve_v{i}(k) = W *par_v{i}(k) + c1 * rand(1) * (Leaderbest_x{i}(k) - par_x{i}(k))+ c2 * rand(1) * (Leaderworst_x{i}(k)- par_x{i}(k));
                Achieve_v{i}(k) = round(Achieve_v{i}(k));
            end
            
           %% update particle position ��������λ��
           % ���ں��������Сֵ����ά��������ƽ�������Ч��һ���������ǣ�����SNP�����ϣ����������SNP֮���໥�����Ļ���
           % �����һ��������������Ӹ��¡�
           % ����Ҫע�⣬��ʵ����ʵSNP����֮��Ļ���Ч�����Ŀǰ�ķ���������Ч��Ҫ�á�ԭ���Ƿ�������ͨ��������SNP֮�������
           % �������ԣ�����LD���ߵ����Ϳ�ȡ���ʵ����������Щ���ԣ�ʹ���������������ƽ�����档
           % ��������ڷ��������Ϸ������������ԵĻ���Ч��Ҳ��úܺá�
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
    
    %% ��ת����
    %% ��¼ÿ�ε���SNPȨ�� 
    % ����ȡ��ǰ���к���SNP��MI�����ֵΪ�䵱ǰȨ�أ���ѡһ
    for i=1:1:ParticleNumber
        for j=1:1:Order
            if SNPWeight(t,par_x{i}(j))<par_mut(i)
                SNPWeight(t,par_x{i}(j))=par_mut(i);
            end
        end
    end
    
    % ����ȡ��ǰ���к���SNP��MI���ۼ�ֵΪ�䵱ǰȨ�أ���ѡһ
    % ����ѡ�������ԭ���ǣ�Ȩ�ر仯�ٶ������ӿ졣��Ȼǰ��������ӱ������������������ô����ĸ������ӣ���Ӧ�ø�������
    % ���š����ڶ��ȡ����SNP������SNP���MIֵ�����ߵ�SNP��ͨ���ۼ���Ȩ�������ߡ���ô�䷭ת�Ŀ����Ծͱ�С����������
    % ����Ŀ����Ա��
    % for i=1:1:ParticleNumber
    %     for j=1:1:Order
    %         SNPWeight(t,par_x{i}(j))=SNPWeight(t,par_x{i}(j))+par_mut(i);
    %     end
    % end
    
    %% ��һ��SNPȨ��
    SUM_SNPWeight=sum(SNPWeight);
    
    minValue = min(SUM_SNPWeight);
    maxValue = max(SUM_SNPWeight);
    SNP_W_Normal = (SUM_SNPWeight - minValue) / (maxValue - minValue);
 
    %% ����ÿ��SNP�ķ�ת���ʲ���������λ��
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

%% ����
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
% Compute Mutual Information between class label and a SNP set.�������ǩ��SNP��֮��Ļ���Ϣ
% MI(X;C)=H(X)+H(C)-H(X,C), where C is the class label and X is a SNP set.
data=pts(:,factor);
MI=Func_MutualInfo(double(data)',double(class));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MI=Func_MutualInfo(data,labels)
% Jan.19 2007 By Guoqiang Yu
%--- INPUT
% data:  features x samples; note that the ��������
% value of zero indicates a missing data ֵΪ0��ʾȱʧ����
% labels: the label matrix of each sample;  ÿ�������ı�ǩ����
%--- OUTPUT
% MI: Information provided by the features �����ṩ����Ϣ
%%
data=data-min(data(:));          %minimum data is 0
labels=labels-min(labels(:))+1; %minimum label is 1

Num_Label=max(labels(:));
Num_DataType=max(data(:))+1;
[Num_Feature,Num_Sample]=size(data);

%Entropy of the random variable Label �������Label����
H_Label=hist(labels,Num_Label);
P_Label=H_Label/Num_Sample;
Entropy_Label=-P_Label*log2(P_Label'+eps);

%Special dealing with the case of small Num_Feature ���⴦��СNum_Feature�����
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

%Larger Feature Number follows the following procedure ������������ѭ����
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

% ����������ת��Ϊ�ַ���
A_str = cellfun(@num2str, Epistasis, 'UniformOutput', false);
B_str = cellfun(@num2str, Epistasis_temp, 'UniformOutput', false);

% �ϲ���ȥ���ظ�Ԫ��
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

