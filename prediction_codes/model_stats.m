function [stats,corrmat,null_distribution]=model_stats(task,pred,ctx,num_perms)

% %  Created for usage in: Tik, N, Gal,S, Bernstein-Eliav, M, Tavor,
% I. Towards a generalized AI framework for predicting task-evoked brain
% activity from resting-state connectivity (2022)

%Input arguments:
% task: task contrast map
% pred: predicted map
% ctx: cortical vertices for testing (remove low signal vertices)
% num_perms: number of permutaions to test the null hypothesis


% Outputs:
% stats: success measures by the following order: [diagonal_mean
% (accuracy);diagonal_median; off_diagonal_mean; off_diagonal_median; 
% diagonality_index(spcificity) ; accuracy permutation p-val; 
% specificity permutation p-val; 


% Compute correlation between actual and predicted maps
[corrmat,pval]=corr(normalise(task(ctx,:)),normalise(pred(ctx,:)));

num_subs=size(corrmat,1);

% find diagonal and off diagonal values
diagon=diag(corrmat);
other=corrmat(~eye(num_subs));% change number if necessary
diagonality_index=mean(diagon)-mean(other);

% calcualte statistics of diagon vs other
diagonal_mean=mean(diagon);
other_mean=mean(other);
diagonal_median=median(diagon);
other_median=median(other);
%diagonal_std=std(diagon);
%other_std=std(other);

% significance and diagonality index using a permutaion test
diagonality_index_perm=zeros(num_perms,1);
all_perm_means=zeros(num_perms,1);

for i=1:num_perms
    disp(i)
    rand_other=randperm(num_subs^2-num_subs,num_subs);
    rand_all=randperm(num_subs^2,num_subs);
    all_perm_means(i)=mean(corrmat(rand_all));
    other_perm_mean=mean(other(rand_other));
    diagonality_index_perm(i)=all_perm_means(i)-other_perm_mean;
end

%diagonality_index=mean(diagonality_index_perm);
all_means=[diagonal_mean; all_perm_means];
all_means_sorted=sort(all_means,'descend');
diagonal_score=find(all_means_sorted==diagonal_mean);
pPerm=diagonal_score/(num_perms+1);

all_DIs=[diagonality_index; diagonality_index_perm];
all_DIs_sorted=sort(all_DIs,'descend');
DI_score=find(all_DIs_sorted==diagonality_index);
pPerm_DI=DI_score/(num_perms+1);


% Output stats
stats=[diagonal_mean; diagonal_median; other_mean; other_median; diagonality_index ;pPerm ; pPerm_DI];

% Create and plot a histogram for the true and permuted diagonals

figure;
set(gcf,'color','w')
hold on
null_distribution=histogram(all_perm_means);
plot([diagonal_mean diagonal_mean],[0 max(null_distribution.BinCounts)],'color','r','linewidth',2)
set(gca,'FontSize',16)
hold off

figure;
set(gcf,'color','w')
hold on
null_distribution_DI=histogram(diagonality_index_perm,'Facecolor',[0 0.8 0]);
plot([diagonality_index diagonality_index],[0 max(null_distribution_DI.BinCounts)],'color','r','linewidth',2)
set(gca,'FontSize',16)
end