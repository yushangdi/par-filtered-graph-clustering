function [X,cmap,Lh]=DrawPMFG(Rpm,F,WeightFlag);


[d dt pred]=bfs(Rpm,1);

if sum(d==-1)>0;
    [ci sizes]=components(Rpm);
    [ss is]=sort(sizes);clear ss
    Wdj=Rpm;
    for n=1:[length(sizes)-1];
        v1=find(ci==is(n));v2=find(ci==is(n+1));
        %v1=v1(ceil(n*rand(1)));v2=v2(ceil(n*rand(1)));
        v1=v1(ceil(length(v1)*rand(1)));v2=v2(ceil(length(v2)*rand(1)));
        Wdj(v1,v2)=1;
        Wdj(v2,v1)=1;
        clear v1 v2
    end
    X=SortLayout(Wdj);
    clear Wdj
    
else
    X=SortLayout(Rpm);
end

% if min(F)<0;
%     nF=(F+abs(min(F)))/[max(F)-min(F)];
% else
%     nF=F/max(F);
% end

if WeightFlag==0;
    figure;
    gplot(Rpm,X,'k-');
else
    [i j]=find(triu(Rpm)~=0);
    figure;hold on;
    for n=1:length(i);
        if Rpm(i(n),j(n))>0;
            plot(X([i(n) j(n)],1),X([i(n) j(n)],2),'k-','LineWidth',Rpm(i(n),j(n))*3);
        else
            plot(X([i(n) j(n)],1),X([i(n) j(n)],2),'k:');
        end
    end
    hold off
end

cmap=colormap(hsv(length(unique(F))));
fvec=unique(F);
gcf;
hold on
for n=1:length(fvec);
    %vec=[nF(n) 0.5 0];
    vec=cmap(fvec(n),:);
    plot(X(F==fvec(n),1),X(F==fvec(n),2),'o','markersize',7,'markeredgecolor',vec,'markerfacecolor',vec);
    child_handles = findobj(gcf,'Type','line');
    Lh(n)=child_handles(1);
    clear vec child_handles
end
hold off

%%
function X=SortLayout(A);

A=double(A~=0);

N = size(A,1);
X = [];
X0  = zeros(N,2);
err = 1;
itt = 0;
while err > 1e-4  & itt < 100
    itt = itt +1;
    [X] = kamada_kawai_spring_layout(A,'progressive',X);
    err = (sum(sum((X-X0).^2,1),2)/N);
    X0 = X;
end

