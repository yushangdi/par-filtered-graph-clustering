function [T8,Rpm,Adjv,Dpm,Mv]=DBHT(D,S);
% Perform DBHT clustering, a deterministic technique which only requires a
% similarity matrix, and respective dissimilarity matrix. It makes
% extensive use of graph-theoretic filtering technique called Planar
% Maximally Filtered Graph (PMFG).
%
% Function call: [T8,Rpm,Adjv,Dpm,Mv]=DBHT(D,S);
% 
% Input
%
% D = NxN dissimilarity matrix.
% S = NxN similarity matrix.
%
% Output
%
% T8 = Nx1 cluster membership vector. 
% Rpm = NxN adjacency matrix of PMFG. 
% Adjv = Bubble cluster membership matrix from BubbleCluster8. 
% Dpm = NxN shortest path length matrix of PMFG
% Mv = NxNb bubble membership matrix. Nv(n,bi)=1 indicates vertex n is a
% vertex of bubble bi. 
%
% Note: In order to calculate the respective DBHT hierarchy, use the
% function call Z=HierarchyConstruct4(Rpm,Dpm,T8,Adjv,Mv);
tic;
Rpm=doPMFG(S);
t = toc;
fprintf('PMFG time: %.5f  \n',t)

fprintf('Outputing PMFG....');
[r,c,val] = find(Rpm); 
fprintf('%f \n', sum(val));
% save(strcat(outputdirname, dataset ,"-P-0"),'r','c','val');

Apm=Rpm;Apm(Apm~=0)=D(Apm~=0);
tic;
Dpm=all_shortest_paths(Apm);
t = toc;
fprintf('shortest path time: %.5f  \n',t)
tic;
[H1,Hb,Mb,CliqList,Sb]=CliqHierarchyTree2(Rpm,'uniqueroot');
clear H1 Sb

Mb=Mb(1:size(CliqList,1),:);
Mv=[];
for n=1:size(Mb,2);
    vec=sparse(size(Rpm,1),1);
    vec(unique(CliqList((Mb(:,n)~=0),:)))=1;
    Mv=[Mv vec];
end
[Adjv,T8]=BubbleCluster8(Rpm,Dpm,Hb,Mb,Mv,CliqList);
t = toc;
fprintf('bubble time: %.5f  \n',t)