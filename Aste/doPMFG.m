function PMFG = doPMFG(W)
%
%  Calculaed the PMFG graph from a matrix of weights W. 
%  W should be sparse, square and symmetric.
%  uses "matlab_bgl" package from 
%  http://www.stanford.edu/%7Edgleich/programs/matlab_bgl/index.htmll
%

N = size(W,1);
if N == 1
    PMFG = sparse(1);
    return
end
[i,j,w] = find(sparse(W));
kk = find(i < j);
ijw= [i(kk),j(kk),w(kk)];
ijw = -sortrows(-ijw,3); %make a sorted list of edges
PMFG = sparse(N,N);
for ii =1:min(6,size(ijw,1)) % the first 6 edges from the list are certanly allowded
    PMFG(ijw(ii,1),ijw(ii,2)) = ijw(ii,3);
    PMFG(ijw(ii,2),ijw(ii,1)) = ijw(ii,3);
end
E = 6; % number of edges in PMFG at this stage
PMFG1 = PMFG;
while( E < 3*(N-2) ) % continue while all edges for a maximal planar graph are inserted
    ii = ii+1;
    PMFG1(ijw(ii,1),ijw(ii,2))=ijw(ii,3); % try to insert the next edge from the sorted list
    PMFG1(ijw(ii,2),ijw(ii,1))=ijw(ii,3); % insert its reciprocal
    if boyer_myrvold_planarity_test(PMFG1~=0) % is the resulting graph planar?
        PMFG = PMFG1; % Yes: insert the edge in PMFG
        E = E+1;
    else
        PMFG1 = PMFG; % No: discard the edge
    end
    if floor(ii/1000)==ii/1000; 
       % save('PMFG.mat','PMFG','ii')
       fprintf('%d    :   %2.2f per-cent done\n',ii,E/(3*(N-2))*100);
        if ii > (N*(N-1)/2)
            fprintf('someting wrong! PMFG not found \n')
            return
        end
    end
    %[E,ijw(ii,1),ijw(ii,2),ijw(ii,3)]
end

% if sum(sum(PMFG<0))>0;
%     indx0=(PMFG<0);
%     vec=PMFG(indx0);
%     rmn0=min(vec);
%     vec=PMFG(PMFG>rmn0);
%     rmn1=min(vec);
%     rmx=max(PMFG(PMFG~=0));
%     PMFG(PMFG~=0)=[PMFG(PMFG~=0)-rmn0+(rmn1-rmn0)/2]/(rmx-rmn0);
%     clear rmn0 rmn1 rmx vec
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% to plot it %%%%%%%
% xy = chrobak_payne_straight_line_drawing(PMFG);
% figure
% gplotwl(PMFG,xy,Labels,'-b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
