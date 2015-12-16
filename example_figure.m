% make simple examples of a few more or less correlated cells
% and see how it affects decoding
%
%
ratemax=10;
ntrials=100;

rng(3) % For reproducibility

nstims=5;
tuning=normpdf([1:nstims],3.0,1.1);
tuning=10*(tuning./max(tuning));

ncells=6;
plotcells=[3 4];

cvs=[.5 0  -0.5];

ll=linspace(0,ratemax,50);
%
fstrs={'b'};
cmap = @(stim) 0.5+([sin(stim),-sin(stim*.2),-sin(stim*1.2)]/2);

% simulate two cells with and without correlation

for stim=1:nstims
    mu{stim}=1+ones(1,ncells)*tuning(stim);
end;

figure (1); clf;
for cv_i=1:numel(cvs)
    
    cv= cvs(cv_i);
    
    sigma=diag(ones(ncells-1,1)*cv,-1)+diag(ones(ncells-1,1)*cv,1)+eye(ncells);
    sigma(ncells,1)=cv;    sigma(1,ncells)=cv;
    sigma_p=(sigma*sigma'); % make positive semidefinite
    
    
    sigma_p=sigma_p./mean(diag(sigma_p));
    %\[L,D] = ldl(sigma);
    
    sigma_p;
    
    x=linspace(0,2*pi,100)';
    circle=[sin(x),cos(x)]*2;
    
    subplot(numel(cvs),3,1+((cv_i-1)*3));
    hold on;
    
    
    for stim=[1 4]
        
        thiscolor=cmap(stim);
        
        rates{stim} = mvnrnd(mu{stim},sigma_p,ntrials);
        
        plot(rates{stim}(:,plotcells(1)),rates{stim}(:,plotcells(2)),[fstrs{1} '.'],'color',thiscolor);
        
        plot(ll,2*normpdf(ll,mu{stim}(plotcells(1)),sigma(plotcells(1),plotcells(1))),fstrs{1},'color',thiscolor);
        plot(2*normpdf(ll,mu{stim}(plotcells(2)),sigma(plotcells(2),plotcells(2))),ll,fstrs{1},'color',thiscolor);
        
        cc=circle*sigma_p(plotcells,plotcells);
        plot( cc(:,1)+mu{stim}(plotcells(1)),cc(:,2)+mu{stim}(plotcells(2)),[fstrs{1} '-'],'color',thiscolor);
        
        %h=histc(r(:,1),ll);
        %plot(ll,h,'b');
        
        xlim([0 ratemax]);
        ylim([0 ratemax]);
        ax=gca;
        ax.XTick=[0 10];
        ax.YTick=[0 10];
        daspect([1 1 1]);
        
    end;
    xlabel(['cell ',num2str(plotcells(1)),' rate']);
    ylabel(['cell ',num2str(plotcells(2)),' rate']);
    title([' Two example cells, cv = ',num2str(cv)]);
    
    subplot(numel(cvs),3,2+((cv_i-1)*3));
    hold on;
    %
    
    trialwidth=.3;
    for stim=1:nstims
        
        rates{stim} = mvnrnd(mu{stim},sigma_p,ntrials);
        
        thiscolor=cmap(stim);
        
        plot(stim+[-trialwidth trialwidth],[1 1]*mu{stim}(1),'color',thiscolor);
        
        xl=linspace(-trialwidth,trialwidth,ncells);
        
        plot(xl+stim,rates{stim}(1,:),'k.','MarkerSize',10);
        
        plot(stim+[-trialwidth trialwidth],[1 1]*mean(rates{stim}(1,:)) ,'k--');
        ax=gca;
        ax.YTick=[0 10];
    end;
    ylim([0 15]);
    ylabel('Firing rates across cells');
    xlabel('Trials/stimuli');
    title([' all ',num2str(ncells),' cells, across all ',num2str(nstims),' stimuli, cv = ',num2str(cv)]);
    %subplot(numel(cvs),3,3+((cv_i-1)*3));
    %hold on;
    
    drawnow;
end;


%% quantify how the correlation affects overall decoding per trial
cvs_dense=linspace(-.5,1,100);
meanerr=[];
ntrials_dense=10000;
for cv_i=1:numel(cvs_dense)
    
    cv= cvs_dense(cv_i);
    
    sigma=diag(ones(ncells-1,1)*cv,-1)+diag(ones(ncells-1,1)*cv,1)+eye(ncells);
    sigma(ncells,1)=cv;    sigma(1,ncells)=cv;
    sigma_p=(sigma*sigma'); % make positive semidefinite
    
    for stim=1:nstims
        rates{stim} = mvnrnd(mu{stim},sigma_p,ntrials_dense);
        
        sigma_p=sigma_p./mean(diag(sigma_p));
        meanerr(cv_i,stim)= mean(abs(mean((rates{stim} - repmat(mu{stim},ntrials_dense,1))')));
    end;
end;

subplot(numel(cvs),3,3);
plot(cvs_dense,mean(meanerr'));
ylim([0 .8]);
ax=gca;
ax.XTick=[-0.5:.2:1];
ax.YTick=[0:0.2:1];
grid on;
ylabel('mean per-trial error from mean');
xlabel('correlation');
title('cv vs. per trial decoding error');

