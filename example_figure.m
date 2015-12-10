%/2
%
%
rng(10) % For reproducibility
ratemax=10;
ntrials=100;

nstims=5;
ncells=4;

ll=linspace(0,ratemax,50);
%
fstrs={'b'};
% simulate two cells with and without correlation

for stim=1:nstims
    mu{stim}=ones(ncells,1)*stim;
end;

cv= -0.3;

sigma =(ones(ncells)*cv) .* (1-eye(ncells)) +eye(ncells);

x=linspace(0,2*pi,100)';
circle=[sin(x),cos(x)]*2;

figure (1); clf;
subplot(2,2,1);
hold on;
for stim=1:nstims
    
    thiscolor=0.5+([sin(stim),-sin(stim*.6),-sin(stim*1.2)]/2);
    
    rates{stim} = mvnrnd(mu{stim},sigma,ntrials);
    
    plot(rates{stim}(:,1),rates{stim}(:,2),[fstrs{1} '.'],'color',thiscolor);
    
    plot(ll,normpdf(ll,mu{stim}(1),sigma(1,1)),fstrs{1},'color',thiscolor);
    plot(normpdf(ll,mu{stim}(2),sigma(2,2)),ll,fstrs{1},'color',thiscolor);
    
    cc=circle*sigma(1:2,1:2);
    plot( cc(:,1)+mu{stim}(1),cc(:,2)+mu{stim}(2),[fstrs{1} '-'],'color',thiscolor);
    
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
    
    thiscolor=0.5+([sin(stim),-sin(stim*.6),-sin(stim*1.2)]/2);
    
    plot(stim+[-trialwidth trialwidth],[1 1]*mu{stim}(1),'color',thiscolor);
    
    ll=linspace(-trialwidth,trialwidth,ncells);
    plot(ll+stim,rates{stim}(1,:),'k.');
    
end;