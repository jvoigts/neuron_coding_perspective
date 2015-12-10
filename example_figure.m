%/2
%
%
rng(10) % For reproducibility
ratemax=10;
ntrials=100;

nstims=5;
ncells=6;
plotcells=[4 5]; 

ll=linspace(0,ratemax,50);
%
fstrs={'b'};
% simulate two cells with and without correlation

for stim=1:nstims
    mu{stim}=1+ones(1,ncells)*stim;
end;

cv= .9;

sigma=diag(ones(ncells-1,1)*cv,-1)+diag(ones(ncells-1,1)*cv,1)+eye(ncells);
sigma_p=(sigma*sigma')/2; % make positive semidefinite

x=linspace(0,2*pi,100)';
circle=[sin(x),cos(x)]*2;

figure (1); clf;
subplot(2,2,1);
hold on;

for stim=[2 5]
    
    thiscolor=0.5+([sin(stim),-sin(stim*.6),-sin(stim*1.2)]/2);
    
    rates{stim} = mvnrnd(mu{stim},sigma_p,ntrials);
    
    plot(rates{stim}(:,plotcells(1)),rates{stim}(:,plotcells(2)),[fstrs{1} '.'],'color',thiscolor);
    
    plot(ll,normpdf(ll,mu{stim}(plotcells(1)),sigma(plotcells(1),plotcells(1))),fstrs{1},'color',thiscolor);
    plot(normpdf(ll,mu{stim}(plotcells(2)),sigma(plotcells(2),plotcells(2))),ll,fstrs{1},'color',thiscolor);
    
    cc=circle*sigma_p(plotcells,plotcells);
    plot( cc(:,1)+mu{stim}(plotcells(1)),cc(:,2)+mu{stim}(plotcells(2)),[fstrs{1} '-'],'color',thiscolor);
    
    %h=histc(r(:,1),ll);
    %plot(ll,h,'b');
    
    xlim([0 ratemax]);
    ylim([0 ratemax]);
    daspect([1 1 1]);
end;

subplot(2,2,2);
hold on;
%

trialwidth=.3;
for stim=1:nstims
    
    rates{stim} = mvnrnd(mu{stim},sigma_p,ntrials);
    
    thiscolor=0.5+([sin(stim),-sin(stim*.6),-sin(stim*1.2)]/2);
    
    plot(stim+[-trialwidth trialwidth],[1 1]*mu{stim}(1),'color',thiscolor);
    
    ll=linspace(-trialwidth,trialwidth,ncells);
    
    plot(ll+stim,rates{stim}(1,:),'k.','MarkerSize',10);
    
    plot(stim+[-trialwidth trialwidth],[1 1]*mean(rates{stim}(1,:)) ,'k--');
    
end;