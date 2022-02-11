function [K3,E,clique]=clique3(A);

% Computes the list of 3-cliques. 
% 
% Input
%
% A = NxN sparse adjacency matrix
%
% Output
%
% clique = Nc x 3 matrix. Each row vector contains the list of vertices for
% a 3-clique. 

A=A-diag(diag(A));
A=(A~=0);

A2=A^2;

P=(A2~=0).*(A~=0);

P=sparse(triu(P));

[r,c]=find(P~=0);
K3=cell(length(r),1);
for n=1:length(r);
    i=r(n);j=c(n);
    a=A(i,:).*A(j,:);
    indx=find(a~=0);
    K3{n}=indx;
    N3(n)=length(indx);
end 

E=[r c];
clique=[0 0 0];

for n=1:length(r);
    temp=K3{n};
    for m=1:length(temp);
        candidate=sort([E(n,:) temp(m)],'ascend');
        a=(clique(:,1)==candidate(1));b=(clique(:,2)==candidate(2));c=(clique(:,3)==candidate(3));
        check=(a.*b).*c;
        check=sum(check);
        if check==0;
            clique=[clique;candidate];
        end
        clear candidate check a b c
    end
end
            
[ignore isort]=sort(clique(:,1),'ascend');

clique=clique(isort,:);
clique=clique(2:size(clique,1),:);
