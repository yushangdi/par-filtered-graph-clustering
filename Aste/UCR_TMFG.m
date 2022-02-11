function UCR_TMFG(name, inputdirname, outputdirname)

    % clear
    % close all
    
    %%%%Prepare inputs for DBHT %%%%%%%%%
%     addpath /home/ubuntu/matlab_bgl/matlab_bgl 
    addpath '/Users/sy/Desktop/MIT/MatlabBGL/matlab_bgl'
    outputdirname = strcat(outputdirname, "/");
    load(strcat(inputdirname, name, ".mat"))
    name = strcat(name,"", "-tmfg")
    D = sqrt(2*(1-R));
    
    fprintf('Computing TMFG....\n');
    mTemp = tic;
    [T8,Rpm,Adjv,Dpm,Mv,Z]=DBHTs(D,R);% DBHT with TMFG clustering
    t = toc(mTemp);
    fprintf('time: %.5f  \n',t)

    fprintf('Found %d clusters \n',length(unique(T8)))
    save(strcat(outputdirname, name, "-Z"),'Z','-ascii');
%     exit;
end