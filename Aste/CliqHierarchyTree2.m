%% ClqHierarchyTree2 looks for 3-cliques of a maximal planar graph, then 
% construct hierarchy of the cliques with the definition of 'inside' a
% clique to be a subgraph with smaller size, when the entire graph is
% made disjoint by removing the clique. Refer and cite to:
%
% Won-Min Song, T. Di Matteo, and Tomaso Aste, Nested hierarchies in planar
% graphs, Discrete Applied Mathematics, Volume 159, Issue 17, 28 October 2011, Pages 2135-2146.
%
% Function call: [H1,Hb,Mb,CliqList,Sb]=CliqHierarchyTree2(Apm,method1);
%
% Input
% 
% Apm = N x N Adjacency matrix of a maximal planar graph
%
% method = Choose between 'uniqueroot' and 'equalroot'. Assigns
%          connections between final root cliques. Uses Voronoi
%          tesselation between tiling triangles. 
%
% Output   

% H1 = Nc x Nc adjacency matrix for 3-clique hierarchical tree where Nc is the number of 3-cliques
% H2 = Nb x Nb adjacency matrix for bubble hierarchical tree where Nb is the number of bubbles
% Mb = Nc x Nb matrix bubble membership matrix. Mb(n,bi)=1 indicates that 3-clique n belongs to bi bubble.
% CliqList = Nc x 3 matrix. Each row vector lists three vertices consisting a 3-clique in the maximal planar graph.
% Sb = Nc x 1 vector. Sb(n)=1 indicates nth 3-clique is separating. 

function [H1,H2,Mb,CliqList,Sb]=CliqHierarchyTree2(Apm,method1);

N=size(Apm,1);
IndxTotal=1:N;

if issparse(Apm)~=1;
    A=sparse(Apm~=0);
else
    A=(Apm~=0);
end

[K3,E,clique]=clique3(A);
clear K3 E N3

Nc=size(clique,1);
M=sparse(N,Nc);
CliqList=clique;
clear clique

for n=1:Nc;
    
    cliq_vec=CliqList(n,:);
    [T,IndxNot]=FindDisjoint(A,cliq_vec);
    % 1，2: inside/outside of clique, 0: clique inds
    indx1=find(T==1);indx2=find(T==2);indx0=find(T==0);
    
    if length(indx1)>length(indx2);
        indx_s=[indx2(:);indx0];
        clear indx1 indx2;
    else
        indx_s=[indx1(:);indx0];
        clear indx1 indx2
    end
    
    if isempty(indx_s)==1;
        Sb(n)=0;
    else
        Sb(n)=length(indx_s)-3;
    end
    M(indx_s,n)=sparse(1);
    clear Indicator InsideCliq count T Temp cliq_vec IndxNot InsideCliq
end

Pred=BuildHierarchy(M);
Root=find(Pred==0);

for n=1:length(Root);
    Components{n}=find(M(:,Root(n))==1);
end
    
clear n


switch lower(method1)
    
    case 'uniqueroot'
        
        if length(Root)>1;
            Pred=[Pred(:);0];
            Pred(Root)=length(Pred);
        end
        H=sparse(Nc+1,Nc+1);
        for n=1:length(Pred);
            if Pred(n)~=0;
                H(n,Pred(n))=sparse(1);
            end
        end
        H=H+H';
        
    case 'equalroot'
        
        if length(Root)>1;
            RootCliq=CliqList(Root,:);
            Adj=AdjCliq(A,CliqList,Root);
        end
        H=sparse(Nc,Nc);
        for n=1:length(Pred);
            if Pred(n)~=0;
                H(n,Pred(n))=sparse(1);
            end
        end
        if isempty(Pred)~=1;
            H=H+H';H=H+Adj;
        else
            H=[];
        end
     
end
H1=H;

if isempty(H1)~=1;
    [H2,Mb]=BubbleHierarchy(Pred,Sb,A,CliqList);
else
    H2=[];Mb=[];
end

H2=double(H2~=0);

Mb=Mb(1:size(CliqList,1),:);
    
%%
function Pred=BuildHierarchy(M);

Pred=zeros(size(M,2),1);

for n=1:size(M,2);
    
    Children=find(M(:,n)==1);
    ChildrenSum=sum(M(Children,:)); %for each clique, how many Children in them
    Parents=find((ChildrenSum==length(Children))); %cliques that have all children
    Parents=Parents(find(Parents~=n));
    
    if isempty(Parents)~=1;
        ParentSum=sum(M(:,Parents));
        a=find(ParentSum==min(ParentSum));
        if length(a)==1;
            Pred(n)=Parents(a);
        else
            Pred=[];
            break
        end
    else
        Pred(n)=0;
    end
end

%%
function [T,IndxNot]=FindDisjoint(Adj,Cliq);

N=size(Adj,1);
Temp=Adj;
T=zeros(N,1);
IndxTotal=1:N;
IndxNot=find((IndxTotal~=Cliq(1))&(IndxTotal~=Cliq(2))&(IndxTotal~=Cliq(3)));
Temp(Cliq,:)=0;Temp(:,Cliq)=0;
[d dt pred] = bfs(Temp,IndxNot(1));
Indx1=find(d==-1);Indx2=find(d~=-1);
T(Indx1)=1;T(Indx2)=2;T(Cliq)=0;
clear Temp

%%
function Adj=AdjCliq(A,CliqList,CliqRoot);
Nc=size(CliqList,1);
N=size(A,1);
Adj=sparse(Nc,Nc);

Indicator=zeros(N,1);
for n=1:length(CliqRoot);
    Indicator(CliqList(CliqRoot(n),:))=1;
    Indi=[Indicator(CliqList(CliqRoot,1)) Indicator(CliqList(CliqRoot,2)) Indicator(CliqList(CliqRoot,3))];
    adjacent=CliqRoot(find(sum(Indi')==2));
    Adj(adjacent,n)=1;
end
Adj=Adj+Adj';
%%
function [H,Mb]=BubbleHierarchy(Pred,Sb,A,CliqList);

Nc=size(Pred,1);
Root=find(Pred==0);
CliqCount=zeros(Nc,1);
CliqCount(Root)=1;
Mb=[];
k=1;


if length(Root)>1;
    TempVec=sparse(Nc,1);TempVec(Root)=1;
    Mb=[Mb TempVec];
    clear TempVec
end
while sum(CliqCount)<Nc;
    NxtRoot=[];
    for n=1:length(Root);
        DirectChild=find(Pred==Root(n));
        TempVec=sparse(Nc,1);TempVec([Root(n);DirectChild(:)])=1;
        Mb=[Mb TempVec];
        CliqCount(DirectChild)=1;
        for m=1:length(DirectChild);
            if Sb(DirectChild(m))~=0;
                NxtRoot=[NxtRoot;DirectChild(m)];
            end
        end
        clear DirectChild TempVec
    end
    Root=unique(NxtRoot);
    k=k+1;
end

Nb=size(Mb,2);
H=sparse(Nb,Nb);



% if sum(IdentifyJoint==0)==0;
    for n=1:Nb;
        Indx=find(Mb(:,n)==1);
        JointSum=sum(Mb(Indx,:));
        Neigh=find(JointSum>=1);
        H(n,Neigh)=sparse(1);
    end
% else
%     H=[];
% end
H=H+H';H=H-diag(diag(H));