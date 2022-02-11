function Z=HierarchyConstruct4(Rpm,Dpm,Tc,Adjv,Mv);
%
% Constructs intra- and inter-cluster hierarchy by utilizing Bubble
% hierarchy structure of a maximal planar graph, namely Planar Maximally Filtered Graph (PMFG). 
%
% Input

% Rpm = NxN Weighted adjacency matrix of PMFG.
% Dpm = NxN shortest path length matrix of PMFG. 
% Tc = Nx1 cluster membership vector from DBHT clustering. Tc(n)=z_i
% indicate cluster of nth vertex. 
% Adjv = Bubble cluster membership matrix from BubbleCluster8. 
% Mv = Bubble membership of vertices from BubbleCluster8. 
%
% Output
%
% Z = (N-1)x3 linkage matrix, in the same format as the output from matlab
% function 'linkage'. To plot the respective dendrogram, use dendrogram(Z).
% Use 'help linkage' for the details 

N=size(Dpm,1);
kvec=unique(Tc);
LabelVec1=[1:N];LinkageDist=[0];
E=sparse(1:N,Tc,ones(N,1),N,max(Tc));
Z=[];
% Intra-cluster hierarchy construction
for n=1:length(kvec);
    Mc=bsxfun(@times,E(:,kvec(n)),Mv);%Get the list of bubbles which coincide with nth cluster
    Mvv=BubbleMember(Dpm,Rpm,Mv,Mc);%Assign each vertex in the nth cluster to a specific bubble.
    Bub=find(sum(Mvv)>0);%Get the list of bubbles which contain the vertices of nth cluster 
    nc=sum(Tc==kvec(n))-1;
    %Apply the linkage within the bubbles.
    for m=1:length(Bub);
        V=find(Mvv(:,Bub(m))~=0);%Retrieve the list of vertices assigned to mth bubble.
        if length(V)>1;
            dpm=Dpm(V,V);%Retrieve the distance matrix for the vertices in V
            LabelVec=LabelVec1(V);%Initiate the label vector which labels for the clusters.
            LabelVec2=LabelVec1;
            for v=1:(length(V)-1);
                [PairLink,dvu]=LinkageFunction(dpm,LabelVec);%Look for the pair of clusters which produces the best linkage
                LabelVec((LabelVec==PairLink(1))|(LabelVec==PairLink(2)))=max(LabelVec1)+1;%Merge the cluster pair by updating the label vector with a same label.
                LabelVec2(V)=LabelVec;
                Z=DendroConstruct(Z,LabelVec1,LabelVec2,1/nc);
                nc=nc-1;
                LabelVec1=LabelVec2;
                clear PairLink dvu Vect
            end
            clear LabelVec dpm rpm LabelVec2
        end
        clear V 
    end

    V=find(E(:,kvec(n))~=0);
    dpm=Dpm(V,V);
    
    %Perform linkage merging between the bubbles
    LabelVec=LabelVec1(V);%Initiate the label vector which labels for the clusters.
    LabelVec2=LabelVec1;
    for b=1:(length(Bub)-1);
        [PairLink,dvu]=LinkageFunction(dpm,LabelVec);
        %[PairLink,dvu]=LinkageFunction(rpm,LabelVec);
        LabelVec((LabelVec==PairLink(1))|(LabelVec==PairLink(2)))=max(max(LabelVec1))+1;%Merge the cluster pair by updating the label vector with a same label.
        LabelVec2(V)=LabelVec;
        Z=DendroConstruct(Z,LabelVec1,LabelVec2,1/nc);
        nc=nc-1;
        LabelVec1=LabelVec2;
        clear PairLink dvu Vect
    end
    
    clear LabelVec V dpm rpm LabelVec2
end


%Inter-cluster hierarchy construction

LabelVec2=LabelVec1;
dcl=ones(1,length(LabelVec1));
for n=1:(length(kvec)-1);
    [PairLink,dvu]=LinkageFunction(Dpm,LabelVec1);
    %[PairLink,dvu]=LinkageFunction(Rpm,LabelVec);
    LabelVec2((LabelVec1==PairLink(1))|(LabelVec1==PairLink(2)))=max(LabelVec1)+1;%Merge the cluster pair by updating the label vector with a same label.
    dvu=unique(dcl(LabelVec1==PairLink(1)))+unique(dcl(LabelVec1==PairLink(2)));
    dcl((LabelVec1==PairLink(1))|(LabelVec1==PairLink(2)))=dvu;
    Z=DendroConstruct(Z,LabelVec1,LabelVec2,dvu);
    LabelVec1=LabelVec2;
    clear PairLink dvu
end

clear LabelVec1 

if length(unique(LabelVec2))>1;
    disp('Something Wrong in Merging. Check the codes.');
    return
end

%%
function [PairLink,dvu]=LinkageFunction(d,labelvec);

lvec=unique(labelvec);

Links=[];

for r=1:(length(lvec)-1);
    vecr=(labelvec==lvec(r));
    for c=(r+1):length(lvec);
        vecc=(labelvec==lvec(c));
        dd=d((vecr|vecc),(vecr|vecc));
        Links=[Links;[lvec(r) lvec(c) max(dd(dd~=0))]];
        clear vecc
    end
end

[dvu imn]=min(Links(:,3));PairLink=Links(imn,1:2);


%%
function Mvv=BubbleMember(Dpm,Rpm,Mv,Mc);

Mvv=sparse(size(Mv,1),size(Mv,2));
vu=find(sum(Mc')>1);
v=find(sum(Mc')==1);

Mvv(v,:)=Mc(v,:);
% 
% for n=1:length(vu);
%     bub=find(Mc(vu(n),:)~=0);
%     vec=Dpm(:,vu(n));
%     vec=vec'*(Mv(:,bub)*diag(1./sum(Mv(:,bub))));
%     [mn imn]=min(vec);imn=imn(1);
%     Mvv(vu(n),bub(imn))=1;
%     clear vec bub vec mn imn
% end

for n=1:length(vu);
    bub=find(Mc(vu(n),:)~=0);
    vu_bub=sum(bsxfun(@times,Rpm(:,vu(n)),Mv(:,bub)))';
    all_bub=diag(Mv(:,bub)'*Rpm*Mv(:,bub))/2;
    frac=vu_bub./all_bub;
    [mx imx]=max(frac);
    Mvv(vu(n),bub(imx(1)))=1;
    clear v_bub all_bub frac bub vec mx imx
end
    

%%
function Z=DendroConstruct(Zi,LabelVec1,LabelVec2,LinkageDist);


indx=(bsxfun(@eq,LabelVec1',LabelVec2')~=1);
if length(unique(LabelVec1(indx)))~=2;
    disp('Check the codes');
    return
end
Z=[Zi;[sort(unique(LabelVec1(indx))) LinkageDist]];

%%
function [PairLink,dvu]=ClusterLinkageFunction(Dpm,LabelVec,Tc,Adjv);

kvec=unique(LabelVec);
Links=[];

for r=1:(length(kvec)-1);
    [v,ignore]=find(Adjv(:,unique(Tc(LabelVec==kvec(r))))~=0);clear ignore
    vecr=unique(v);clear v
    for c=(r+1):length(kvec);
        [v,ignore]=find(Adjv(:,unique(Tc(LabelVec==kvec(c))))~=0);clear ignore
        vecc=unique(v);clear v
        d=Dpm(union(vecr,vecc),union(vecr,vecc));
        %Links=[Links;kvec(r) kvec(c) mean(d(d~=0))];
        Links=[Links;kvec(r) kvec(c) max(d(d~=0))];
    end
end

[dvu imn]=min(Links(:,3));PairLink=Links(imn,1:2);