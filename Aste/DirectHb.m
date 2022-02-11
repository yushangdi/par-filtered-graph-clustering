function [Hc,Sep]=DirectHb(Rpm,Hb,Mb,Mv,CliqList);
% Computes directions on each separating 3-clique of a maximal planar
% graph, hence computes Directed Bubble Hierarchical Tree (DBHT). 
% 
% Function call: Hc=DirectHb(Rpm,Hb,Mb,Mv,CliqList);
%
% Input
% Rpm = N x N sparse weighted adjacency matrix of PMFG
% Hb = Undirected bubble tree of PMFG
% Mb = Nc x Nb bubble membership matrix for 3-cliques. Mb(n,bi)=1 indicates that
% 3-clique n belongs to bi bubble. 
% Mv = N x Nb bubble membership matrix for vertices. 
% CliqList = Nc x 3 matrix of list of 3-cliques. Each row vector contains
% the list of vertices for a particular 3-clique. 
%
% Output
% Hc = Nb x Nb unweighted directed adjacency matrix of DBHT. Hc(i,j)=1
% indicates a directed edge from bubble i to bubble j. 


Hb=(Hb~=0);
[r,c]=find(triu(Hb)~=0);
CliqEdge=[];
for n=1:length(r);
    CliqEdge=[CliqEdge;r(n) c(n) find((Mb(:,r(n))~=0)&(Mb(:,c(n))~=0))]; % the clique that's in both bubble
end
clear r c

kb=sum(Hb~=0);
Hc=sparse(size(Mv,2),size(Mv,2));
for n=1:size(CliqEdge,1);
    Temp=Hb;
    Temp(CliqEdge(n,1),CliqEdge(n,2))=0;
    Temp(CliqEdge(n,2),CliqEdge(n,1))=0;
    [d dt p]=bfs(Temp,1);
    vo=CliqList(CliqEdge(n,3),:);
    bleft=CliqEdge(n,1:2);bleft=bleft(d(bleft)~=-1);
    bright=CliqEdge(n,1:2);bright=bright(d(bright)==-1);
    [vleft c]=find(Mv(:,(d~=-1))~=0);vleft=setdiff(vleft,vo);
    [vright c]=find(Mv(:,(d==-1))~=0);vright=setdiff(vright,vo);clear c
    left=sum(sum(Rpm(vo,vleft)));
    right=sum(sum(Rpm(vo,vright)));
    if left>right;
        Hc(bright,bleft)=left;
    else
        Hc(bleft,bright)=right;
    end
    clear vleft vright vo Temp bleft bright right left
end

Sep=double((sum(Hc')==0));
Sep((sum(Hc)==0)&(kb>1))=2;
