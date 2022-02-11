function [Adjv,Tc]=BubbleCluster8(Rpm,Dpm,Hb,Mb,Mv,CliqList)
% Obtains non-discrete and discrete clusterings from the bubble topology of
% PMFG. 
% 
% Function call: [Adjv,T8]=BubbleCluster8(Rpm,Dpm,Hb,Mb,Mv,CliqList);
%
% Input
% 
% Rpm = N x N sparse weighted adjacency matrix of PMFG
% Dpm = N x N shortest path lengths matrix of PMFG
% Hb = Undirected bubble tree of PMFG
% Mb = Nc x Nb bubble membership matrix for 3-cliques. Mb(n,bi)=1 indicates that
% 3-clique n belongs to bi bubble. 
% Mv = N x Nb bubble membership matrix for vertices. 
% CliqList = Nc x 3 matrix of list of 3-cliques. Each row vector contains
% the list of vertices for a particular 3-clique. 
%
% Output
% 
% Adjv = N x Nk cluster membership matrix for vertices for non-discrete
% clustering via the bubble topology. Adjv(n,k)=1 indicates cluster
% membership of vertex n to kth non-discrete cluster.
% Tc = N x 1 cluster membership vector. Tc(n)=k indicates cluster
% membership of vertex n to kth discrete cluster.

[Hc,Sep]=DirectHb(Rpm,Hb,Mb,Mv,CliqList);%Assign directions on the bubble tree
N=size(Rpm,1);% Number of vertices in the PMFG
indx=find(Sep==1);% Look for the converging bubbles
Adjv=[];

if length(indx)>1;
    Adjv=sparse(size(Mv,1),length(indx));%Set the non-discrete cluster membership matrix 'Adjv' at default
    
    % Identify the non-discrete cluster membership of vertices by each
    % converging bubble
    for n=1:length(indx);
        [d dt p]=bfs(Hc',indx(n));
        [r c]=find(Mv(:,d~=-1)~=0);
        Adjv(unique(r),n)=1;
        clear d dt p r c
    end

    Tc=zeros(N,1);% Set the discrete cluster membership vector at default
    Bubv=Mv(:,indx);% Gather the list of vertices in the converging bubbles
    cv=find(sum(Bubv')'==1);% Identify vertices which belong to single converging bubbles
    uv=find(sum(Bubv')'>1);% Identify vertices which belong to more than one converging bubbles.
    Mdjv=sparse(N,length(indx));% Set the cluster membership matrix for vertices in the converging bubbles at default
    Mdjv(cv,:)=Bubv(cv,:);% Assign vertices which belong to single converging bubbles to the rightful clusters.
    
    % Assign converging bubble membership of vertices in `uv'
    for v=1:length(uv);
        v_cont=sum(bsxfun(@times,Rpm(:,uv(v)),Bubv))';% sum of edge weights linked to uv(v) in each converging bubble
        all_cont=3*(full(sum(Bubv))-2);% number of edges in converging bubble
        [mx imx]=max(v_cont(:)./all_cont(:));% computing chi(v,b_{alpha})
        Mdjv(uv(v),imx(1))=1;% Pick the most strongly associated converging bubble
    end
    
    [v ci]=find(Mdjv~=0);Tc(v)=ci;clear v ci% Assign discrete cluster memebership of vertices in the converging bubbles.
    
    Udjv=Dpm*(Mdjv*diag(1./sum(Mdjv~=0)));Udjv(Adjv==0)=inf;% Compute the distance between a vertex and the converging bubbles.
    [mn imn]=min(Udjv(sum(Mdjv')==0,:)');% Look for the closest converging bubble
    Tc(Tc==0)=imn;% Assign discrete cluster membership according to the distances to the converging bubbles
    
else
    Tc=ones(N,1); % if there is one converging bubble, all vertices belong to a single cluster
end
