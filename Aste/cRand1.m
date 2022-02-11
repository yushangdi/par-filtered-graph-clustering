function [AR]=cRand1(true_mem,mem)

if nargin==1
    T=true_mem; %contingency table pre-supplied
elseif nargin==2
    %build the contingency table from membership arrays
    R=max(true_mem);
    C=max(mem);
    n=length(mem);N=n;

    %identify the missing labels
    list_t=ismember(1:R,true_mem);
    list_m=ismember(1:C,mem);
    T=Contingency(true_mem,mem);
    T=T(list_t,list_m);
end

%-----------------------calculate Rand index and others----------
n=sum(sum(T));N=n;
C=T;
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(C.^2));      %sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

%---------------------auxiliary functions---------------------
function Cont=Contingency(Mem1,Mem2)

if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end