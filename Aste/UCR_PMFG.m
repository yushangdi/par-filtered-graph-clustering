function UCR_PMFG(name, inputdirname, outputdirname)

    % clear
    % close all

    %%%%Prepare inputs for DBHT %%%%%%%%%
    % addpath /home/ubuntu/matlab_bgl/matlab_bgl 
    addpath '/Users/sy/Desktop/MIT/MatlabBGL/matlab_bgl'
%     inputdirname = "./";
    outputdirname = strcat(outputdirname, "/");
    load(strcat(inputdirname, name, ".mat"))
    name = strcat(name, "-pmfg")
    D = sqrt(2*(1-R));

    fprintf('Computing PMFG....\n');
    mTemp = tic;
    [T8,Rpm,Adjv,Dpm,Mv]=DBHT(D,R);% DBHT  with PMFG clustering
    tic;
    Z=HierarchyConstruct4(Rpm,Dpm,T8,Adjv,Mv);clear Adjv Dpm Mv % DBHT hierarchy
    t = toc;
    fprintf('hierarchy time: %.5f  \n',t)
    t = toc(mTemp);
    fprintf('time: %.5f  \n',t)

    fprintf('Found %d clusters \n',length(unique(T8)))
    save(strcat(outputdirname, name, "-Z"),'Z','-ascii');
%     exit;
end