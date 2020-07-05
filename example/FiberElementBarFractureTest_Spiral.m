clear;
clc;

% Set default to desired settings
set(0,'DefaultAxesFontName','Century Schoolbook');
set(0,'DefaultTextFontName','Century Schoolbook');
set(0,'DefaultAxesFontWeight', 'bold', ...
    'DefaultAxesFontSize', 10, ...
    'DefaultTextFontSize',10, ...
    'DefaultAxesFontAngle', 'normal', ...
    'DefaultTextFontWeight', 'normal');
cardinal = [140,21,21]/256;
beige = [157,149,115]/256;
dard_gray = beige;

%%
DVInputFileName = 'FractureTestPropertiesSpiral.xlsx';

%% Define design variables
ex = importdata(DVInputFileName);
test.ID = ex.data(:,1);
test.TSN = ex.data(:,3);
test.L = ex.data(:,4);
test.D = ex.data(:,5);
test.c = ex.data(:,6);
test.fc = ex.data(:,7);
test.ft = ex.data(:,8);
test.Ec = ex.data(:,9);
test.Et = ex.data(:,10);
test.nsl = ex.data(:,11);
test.dbl = ex.data(:,12);
test.Asl = ex.data(:,13);
test.fyl = ex.data(:,14);
test.ful = ex.data(:,15);
test.TY = test.ful./test.fyl;
test.esy = ex.data(:,16);
test.esh = ex.data(:,17);
test.esu = ex.data(:,18);
test.Esl = ex.data(:,19);
test.Esh = ex.data(:,20);
test.dbt = ex.data(:,21);
test.Ashi = ex.data(:,22);
test.fyt = ex.data(:,23);
test.Est = ex.data(:,24);
test.s = ex.data(:,25);
test.P = ex.data(:,26);
test.bs = ex.data(:,27);
test.R = ex.data(:,28);
test.barlayout = ex.data(:,29);
test.Sy = ex.data(:,30);
test.CR = ex.data(:,31);
test.ER = ex.data(:,32);
test.phsi = ex.data(:,33);
test.syaci = ex.data(:,35);
test.suaci = ex.data(:,36);

%% Compute other modeling parameters
% compute the equivalent lateral buckleproof capacity
test.kt = 8*pi*test.Est.*test.Ashi./3./test.nsl./(test.D-2*test.c-test.dbt);
% the required spring stiffness for buckle modes
test.keff = 0.04*test.Esl.*(pi*test.dbl.^4/64);
test.kelastic = pi^4*test.keff./(test.s.^3);
test.kratio = test.kt./test.kelastic;
test.kmodefactor = [0.5 0.1173 0.0377 0.0154 0.0075 0.0040 0.0024 0.0015 0.0010 0.0007];
test.kequ = test.kelastic*test.kmodefactor;
for i = 1:1:length(test.ID)
    test.modenum(i,1) = max(1,min(find(test.kmodefactor<test.kratio(i,1))));
end
test.betabuckle = min([1*ones(numel(test.ID),1), ...
    max([0.75*ones(numel(test.ID),1), ...
    0.9571+0.1209*test.fyl/60+0.1276*test.TY-0.0698* ...
    test.s./test.dbl],[],2)],[],2);
% buckle strain (per Dhakal and Maekawa equation)
test.ebuckle = test.fyl./test.Esl.*(55-2.3*sqrt(test.fyl/14.5).* ...
    test.modenum.*test.s./test.dbl);
test.sdb = test.s./test.dbl;
test.sdb_eff = test.phsi.*test.sdb;
% deformation properties
Hdb = 0.07;
RminH = 1.46;
% MP curve hardening
test.MPR1 = 3.0*test.fyl./test.fyl; % fix at 3.0
test.MPR2 = exp(3.4439+0.1790*log(test.fyl/60));
test.MPR3 = exp(-1.0389-2.1958*log(test.fyl/60)-1.9579*log(test.TY)- ...
    1.7812*log(test.esu)-0.2175*log(test.dbl));
% fracture model
test.c_mono = exp(-3.7513-1.7857*log(test.esu)+0.2892*log(test.dbl));
test.c_cycl = exp(4.64+1.36*log(test.fyl/60)+ ...
    1.74*log(test.esu)+0.89*log(test.dbl)-0.07*log(RminH/1.46));
test.c_symm = 1.02*test.fyl./test.fyl; % fix at 1.02
test.k1 = exp(2.12-0.66*log(test.dbl));
test.k2 = exp(1.8900-1.3793*log(test.TY)-0.5085*log(test.dbl));
test.b1 = exp(-2.43-1.66*log(test.TY)- ...
    0.42*log(RminH/1.46));
test.b2 = exp(-2.83-0.43*log(test.esu)- ...
    0.88*log(test.sdb_eff)- ...
    0.12*log(RminH/1.46));
% bond slip
test.ub = 12*sqrt(-test.fc*1000)/1000;
test.ld = test.fyl.*test.dbl/4./test.ub;

% define load tag
% 0: loading protocol in tests; 1: FEMA 461; 2: modified FEMA461; 
% 3: ACI 374; 4: modified ACI 374
loadtag = 0;

numIP = 6;
% define specimen tag
tag1 = 1;
tag2 = 26; % <= test.numSec
for currentID = tag1:1:tag2
    %% Create excutable file of design variables
    f = fopen('DesignVariableSpiral.tcl','wt');
    fprintf(f,'set L %6.4f; \n',test.L(currentID,1));
    fprintf(f,'set D %6.4f; \n',test.D(currentID,1));
    fprintf(f,'set c %6.4f; \n',test.c(currentID,1));
    fprintf(f,'set fc %6.4f; \n',test.fc(currentID,1));
    fprintf(f,'set ft %6.4f; \n',test.ft(currentID,1));
    fprintf(f,'set Ec %6.4f; \n',test.Ec(currentID,1));
    fprintf(f,'set Et %6.4f; \n',test.Et(currentID,1));
    fprintf(f,'set nsl %d; \n',test.nsl(currentID,1));
    fprintf(f,'set dbl %6.4f; \n',test.dbl(currentID,1));
    fprintf(f,'set Asl %6.4f; \n',test.Asl(currentID,1));
    fprintf(f,'set fyl %6.4f; \n',test.fyl(currentID,1));
    fprintf(f,'set ful %6.4f; \n',test.ful(currentID,1));
    fprintf(f,'set esy %6.4f; \n',test.esy(currentID,1));
    fprintf(f,'set esh %6.4f; \n',test.esh(currentID,1));
    fprintf(f,'set esu %6.4f; \n',test.esu(currentID,1));
    fprintf(f,'set Esl %6.4f; \n',test.Esl(currentID,1));
    fprintf(f,'set Esh %6.4f; \n',test.Esh(currentID,1));
    fprintf(f,'set dbt %6.4f; \n',test.dbt(currentID,1));
    fprintf(f,'set Ast %6.4f; \n',test.Ashi(currentID,1));
    fprintf(f,'set fyt %6.4f; \n',test.fyt(currentID,1));
    fprintf(f,'set Est %6.4f; \n',test.Est(currentID,1));
    fprintf(f,'set s %6.4f; \n',test.s(currentID,1));
    fprintf(f,'set P %6.4f; \n',-abs(test.P(currentID,1)));
    fprintf(f,'set betabuckle %6.4f; \n',test.betabuckle(currentID,1));
    fprintf(f,'set modenum %6.4f; \n',test.modenum(currentID,1));
    fprintf(f,'set MPR1 %6.4f; \n',test.MPR1(currentID,1));
    fprintf(f,'set MPR2 %6.4f; \n',test.MPR2(currentID,1));
    fprintf(f,'set MPR3 %6.4f; \n',test.MPR3(currentID,1));
    fprintf(f,'set c_mono %6.4f; \n',test.c_mono(currentID,1));
    fprintf(f,'set c_cycl %6.4f; \n',test.c_cycl(currentID,1));
    fprintf(f,'set c_symm %6.4f; \n',test.c_symm(currentID,1));
    fprintf(f,'set k1 %6.4f; \n',test.k1(currentID,1));
    fprintf(f,'set k2 %6.4f; \n',test.k2(currentID,1));
    fprintf(f,'set b1 %6.4f; \n',test.b1(currentID,1));
    fprintf(f,'set b2 %6.4f; \n',test.b2(currentID,1));
    fprintf(f,'set bs %6.4f; \n',test.bs(currentID,1));
    fprintf(f,'set BR %6.4f; \n',test.R(currentID,1));
    fprintf(f,'set barlayout %d; \n',test.barlayout(currentID,1));
    fprintf(f,'set Sy %6.4f; \n',test.syaci(currentID,1));
    fprintf(f,'set Su %6.4f; \n',test.suaci(currentID,1));
    fprintf(f,'set numIntgrPts %d; \n',numIP);
    fprintf(f,'set DFTag 0; \n'); % no fracture
    if currentID == 11
        fprintf(f,'set LWTag 1; \n'); % lightweight concrete
    else
        fprintf(f,'set LWTag 0; \n');
    end
    fprintf(f,'set CR %6.4f; \n',test.CR(currentID,1));
    fprintf(f,'set ER %6.4f; \n',test.ER(currentID,1));
    fclose(f);
    
    %% Create loading file
    ex = importdata(fullfile('./TestData', ...
        [num2str(test.TSN(currentID,1)) '.txt']));
    disp = ex(:,1);
    % maximum length per step (0.1% of column height)
    delta_disp = 0.0005*test.L(currentID,1);
    disp_m = disp(1,1);
    if currentID > 18 && currentID <25
        v = zeros(length(disp),1);
    else
        v = ex(:,2);
    end
    v_m = v(1,1);
    for ldtag = 2:1:length(disp)
        cur_delta = abs(disp(ldtag)-disp(ldtag-1));
        if cur_delta < delta_disp
            disp_m = [disp_m;disp(ldtag)];
            v_m = [v_m;v(ldtag)];
        else
            if (disp(ldtag-1)<disp(ldtag))
                temp_disp = [disp(ldtag-1)+delta_disp:delta_disp:disp(ldtag)]';
            else
                temp_disp = [disp(ldtag-1)-delta_disp:-delta_disp:disp(ldtag)]';
            end
            if temp_disp(end) ~= disp(ldtag)
                temp_disp = [temp_disp;disp(ldtag)];
            end
            temp_v = interp1([disp(ldtag-1);disp(ldtag)], ...
                [v(ldtag-1);v(ldtag)],temp_disp);
            disp_m = [disp_m;temp_disp];
            v_m = [v_m;temp_v];
            clear temp_disp temp_v;
        end
    end
    % Remove unmoving points
    xtag = 1;
    while xtag < length(disp_m)
        if disp_m(xtag) == disp_m(xtag+1)
            disp_m = [disp_m(1:xtag);disp_m(xtag+2:end)];
            v_m = [v_m(1:xtag);v_m(xtag+2:end)];
        else
            xtag = xtag+1;
        end
    end
    % P-Delta re-correction
    if currentID > 8 && currentID < 24
        v_m = v_m-test.P(currentID,1)*disp_m/test.L(currentID,1);
    end
    test.disp{currentID,1} = disp_m;
    test.shear{currentID,1} = v_m;
    f = fopen('LoadingParameterSpiral.tcl','wt');
    fprintf(f,'set numIncr %d; \n',800);
    fprintf(f,'set numEle %d; \n',1);
    fprintf(f,'set Dincr %8.6f; \n',0.001);
    fprintf(f,'set LoadHistory { \n');
    for ldstep = 1:1:length(disp_m)
        fprintf(f,'%8.4f \n',disp_m(ldstep));
    end
    fprintf(f,'}; \n');
    fclose(f);
    
    %%
    ! OpenSees.exe ColumnCyclicTestSpiral.tcl
    
    %% Collect data
    ex = importdata(fullfile('./CyclicOutputSpiral/','disp.out'));
    simu.disp{currentID,1} = ex(:,2);
    ex = importdata(fullfile('./CyclicOutputSpiral/','force.out'));
    simu.shear{currentID,1} = -ex(:,2);
    % concrete strain
    ex = importdata(fullfile('./CyclicOutputSpiral/','ConcrCoreBot.out'));
    simu.stress.concrcorebot{currentID,1} = ex(:,2);
    simu.strain.concrcorebot{currentID,1} = ex(:,3);
    ex = importdata(fullfile('./CyclicOutputSpiral/','ConcrCoreTop.out'));
    simu.stress.concrcoretop{currentID,1} = ex(:,2);
    simu.strain.concrcoretop{currentID,1} = ex(:,3);
    ex = importdata(fullfile('./CyclicOutputSpiral/','ConcrCoverBot.out'));
    simu.stress.concrcoverbot{currentID,1} = ex(:,2);
    simu.strain.concrcoverbot{currentID,1} = ex(:,3);
    ex = importdata(fullfile('./CyclicOutputSpiral/','ConcrCoverTop.out'));
    simu.stress.concrcovertop{currentID,1} = ex(:,2);
    simu.strain.concrcovertop{currentID,1} = ex(:,3);
    % steel strain
    ex = importdata(fullfile('./CyclicOutputSpiral/','SS0.out'));
    simu.stress.steelbot{currentID,1} = ex(:,2);
    simu.strain.steelbot{currentID,1} = ex(:,3);    
    ex = importdata(fullfile('./CyclicOutputSpiral/', ...
        ['SS' num2str(round(test.nsl(currentID,1)/2)) '.out']));
    simu.stress.steeltop{currentID,1} = ex(:,2);
    simu.strain.steeltop{currentID,1} = ex(:,3);
    % bar-slip and shear drift contribution
    ex = importdata('./CyclicOutputSpiral/BarSlipForce.out');
    simu.BarSlipForce{currentID,1} = -ex(:,2);
    ex = importdata('./CyclicOutputSpiral/BarSlipDisp.out');
    simu.BarSlipDrift{currentID,1} = -ex(:,3);
    simu.ShearDrift{currentID,1} = -ex(:,4);
    % bar-slip concrete
    ex = importdata('./CyclicOutputSpiral/BarSlipConcrCoreBot.out');
    simu.stress.BarSlipconcrcorebot{currentID,1} = ex(:,2);
    simu.strain.BarSlipconcrcorebot{currentID,1} = ex(:,3);
    ex = importdata('./CyclicOutputSpiral/BarSlipConcrCoreTop.out');
    simu.stress.BarSlipconcrcoretop{currentID,1} = ex(:,2);
    simu.strain.BarSlipconcrcoretop{currentID,1} = ex(:,3);
    
    % compute fracture index
    es_hist = [simu.strain.steelbot{currentID,1}, ...
        simu.strain.steeltop{currentID,1}];
    ss_hist = [simu.stress.steelbot{currentID,1}, ...
        simu.stress.steeltop{currentID,1}];
    num_rlz = 200;
    [FI,rlz] = DuctileFractureModel(es_hist,ss_hist, ...
        test.c_mono(currentID,1),test.c_symm(currentID,1), ...
        test.c_cycl(currentID,1),test.Esl(currentID,1), ...
        test.esu(currentID,1),test.k1(currentID,1),test.k2(currentID,1), ...
        test.dbl(currentID,1),test.b1(currentID,1),test.b2(currentID,1),num_rlz,0);
    % fracture probability (empirical)
    for i = 1:1:2
        ftag = find(max(FI{i,1})>=1);
        for j = 1:1:length(ftag)
            ftag(j) = min(find(FI{i,1}(:,ftag(j))>=1.0));
        end
        [frac_a,frac_b] =ecdf(ftag);
        frac_step{i,1} = [frac_a*length(frac_a)/num_rlz,frac_b];
    end
    FI_1 = FI{1,1}(:,1);
    for mtag = 1:1:length(FI_1)
        FI_1(mtag,1) = max(FI_1(1:mtag));
    end
    pf_1 = logncdf(FI_1,0,1.2*sqrt(FI_1));
    FI_2 = FI{2,1}(:,1);
    for mtag = 1:1:length(FI_2)
        FI_2(mtag,1) = max(FI_2(1:mtag));
    end
    pf_2 = logncdf(FI_2,0,1.2*sqrt(FI_2));

    % plot
    figure;
    plot(test.disp{currentID,1}/test.L(currentID,1), ...
        test.shear{currentID,1},'-','color','k','linewidth',1);
    hold on;
    plot(simu.disp{currentID,1}/test.L(currentID,1), ...
        simu.shear{currentID,1},'--','color',beige,'linewidth',1);
    xlabel('Lateral drift');
    ylabel('Base shear (kip)');
    grid on;
    legend('Test','Simulation','location','southeast');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_hyst.fig']);
    close(gcf);
    % Bar-slip hysteresis
    figure;
    plot(simu.BarSlipDrift{currentID,1}, ...
        simu.shear{currentID,1},'-k','linewidth',1);
    xlabel('Bar-slip drift');
    ylabel('Base shear (kip)');
    grid on;
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_hystbarslip.fig']);
    close(gcf);
    % Shear hysteresis
    figure;
    plot(simu.ShearDrift{currentID,1}, ...
        simu.shear{currentID,1},'-k','linewidth',1);
    xlabel('Shear drift');
    ylabel('Base shear (kip)');
    grid on;
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_hystshear.fig']);
    close(gcf);
    % Drift vs. strain
    figure;
    plot(simu.disp{currentID,1}/test.L(currentID,1), ...
        simu.strain.steelbot{currentID,1},'-k','linewidth',1);
    hold on;
    plot(simu.disp{currentID,1}/test.L(currentID,1), ...
        simu.strain.steeltop{currentID,1},'-','color',beige,'linewidth',1);
    grid on;
    xlabel('Lateral drift');
    ylabel('Steel strain');
    legend('Negative-side rebar','Positive-side rebar','location','north');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_driftstrain.fig']);
    close(gcf);
    % Steel strain-stress
    figure;
    plot(simu.strain.steelbot{currentID,1}, ...
        simu.stress.steelbot{currentID,1},'-k','linewidth',1);
    hold on;
    plot(simu.strain.steeltop{currentID,1}, ...
        simu.stress.steeltop{currentID,1},'-','color',beige,'linewidth',1);
    grid on;
    xlabel('Steel strain');
    ylabel('Steel stress (ksi)');
    legend('Negative-side rebar','Positive-side rebar','location','north');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_steelstresstrain.fig']);
    close(gcf);
    % Fracture probability
    figure;
    subplot(3,1,1);
    plot(disp_m/test.L(currentID,1),'-k');
    hold on;
    plot(simu.disp{currentID,1}/test.L(currentID,1),'--', ...
        'color',beige);
    grid on;
    xlabel('Load step');
    ylabel('Lateral drift');
    legend('Test','Simulation','location','northwest');
    xlim([0 length(disp_m)]);
    subplot(3,1,2);
    plot(simu.strain.steelbot{currentID,1},'-k');
    hold on;
    plot(simu.strain.steeltop{currentID,1},'-','color',beige);
    grid on;
    xlabel('Load step');
    ylabel('Steel strain');
    legend('Negative-side rebar','Positive-side rebar','location','northwest');
    xlim([0 length(disp_m)]);
    subplot(3,1,3);
    stairs([frac_step{1,1}(:,2);length(disp_m)], ...
        [frac_step{1,1}(:,1);frac_step{1,1}(end,1)],'-k');
    hold on;
    stairs([frac_step{2,1}(:,2);length(disp_m)], ...
        [frac_step{2,1}(:,1);frac_step{2,1}(end,1)],'-','color',beige);
%     plot(pf_1,'--k');
%     plot(pf_2,'--','color',beige);
    grid on;
    xlabel('Load step');
    ylabel('Fracture probability');
    legend('Negative-side rebar (empirical)', ...
        'Positive-side rebar (empirical)','location','northwest');
    xlim([0 length(disp_m)]);
    ylim([0 1]);
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_fractureprob.fig']);
    close(gcf);
    % Fracture index history
    for i = 1:1:2
        for j = 1:1:num_rlz
            if max(FI{i,1}(:,j)) >= 1
                FI{i,1}(min(find(FI{i,1}(:,j)>=1)):end,j) = 1;
            end
        end
    end
    figure;
    subplot(3,1,1);
    plot(disp_m/test.L(currentID,1),'-k');
    hold on;
    plot(simu.disp{currentID,1}/test.L(currentID,1),'--', ...
        'color',beige);
    grid on;
    xlabel('Load step');
    ylabel('Lateral drift');
    legend('Test','Simulation','location','northwest');
    xlim([0 length(disp_m)]);
    subplot(3,1,2);
    plot(FI{1,1},'-','color',beige);
    hold on;
    plot(median(FI{1,1},2),'-k','linewidth',1);
    grid on;
    xlim([0 length(disp_m)]);
    xlabel('Load step');
    ylabel('Fracture index (Negative-side rebar)');
    legend('Monte-Carlo sample','Median','location','northwest');
    subplot(3,1,3);
    plot(FI{2,1},'-','color',beige);
    hold on;
    plot(median(FI{2,1},2),'-k','linewidth',1);
    grid on;
    xlim([0 length(disp_m)]);
    xlabel('Load step');
    ylabel('Fracture index (Positive-side rebar)');
    legend('Monte-Carlo sample','Median','location','northwest');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_fihist.fig']);
    close(gcf);
    % drift decomposition
    figure;
    area(simu.disp{currentID,1}/test.L(currentID,1), ...
        'FaceColor',beige,'linewidth',0.5)
    hold on;
    area(simu.disp{currentID,1}/test.L(currentID,1)- ...
        simu.ShearDrift{currentID,1},'FaceColor',[0.9,0.9,0.9],'linewidth',0.5)
    area(simu.disp{currentID,1}/test.L(currentID,1)- ...
        simu.ShearDrift{currentID,1}- ...
        simu.BarSlipDrift{currentID,1},'FaceColor',[0.7,0.7,0.7],'linewidth',0.5)
    grid on;
    xlabel('Load step');
    ylabel('Lateral drift');
    legend('Shear','Bar-slip','Flexure');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_driftdecomp.fig']);
    close(gcf);
    
    
    %% Simulate fracture (median model)
    f = fopen('DesignVariableSpiral.tcl','wt');
    fprintf(f,'set L %6.4f; \n',test.L(currentID,1));
    fprintf(f,'set D %6.4f; \n',test.D(currentID,1));
    fprintf(f,'set c %6.4f; \n',test.c(currentID,1));
    fprintf(f,'set fc %6.4f; \n',test.fc(currentID,1));
    fprintf(f,'set ft %6.4f; \n',test.ft(currentID,1));
    fprintf(f,'set Ec %6.4f; \n',test.Ec(currentID,1));
    fprintf(f,'set Et %6.4f; \n',test.Et(currentID,1));
    fprintf(f,'set nsl %d; \n',test.nsl(currentID,1));
    fprintf(f,'set dbl %6.4f; \n',test.dbl(currentID,1));
    fprintf(f,'set Asl %6.4f; \n',test.Asl(currentID,1));
    fprintf(f,'set fyl %6.4f; \n',test.fyl(currentID,1));
    fprintf(f,'set ful %6.4f; \n',test.ful(currentID,1));
    fprintf(f,'set esy %6.4f; \n',test.esy(currentID,1));
    fprintf(f,'set esh %6.4f; \n',test.esh(currentID,1));
    fprintf(f,'set esu %6.4f; \n',test.esu(currentID,1));
    fprintf(f,'set Esl %6.4f; \n',test.Esl(currentID,1));
    fprintf(f,'set Esh %6.4f; \n',test.Esh(currentID,1));
    fprintf(f,'set dbt %6.4f; \n',test.dbt(currentID,1));
    fprintf(f,'set Ast %6.4f; \n',test.Ashi(currentID,1));
    fprintf(f,'set fyt %6.4f; \n',test.fyt(currentID,1));
    fprintf(f,'set Est %6.4f; \n',test.Est(currentID,1));
    fprintf(f,'set s %6.4f; \n',test.s(currentID,1));
    fprintf(f,'set P %6.4f; \n',-abs(test.P(currentID,1)));
    fprintf(f,'set betabuckle %6.4f; \n',test.betabuckle(currentID,1));
    fprintf(f,'set modenum %6.4f; \n',test.modenum(currentID,1));
    fprintf(f,'set MPR1 %6.4f; \n',test.MPR1(currentID,1));
    fprintf(f,'set MPR2 %6.4f; \n',test.MPR2(currentID,1));
    fprintf(f,'set MPR3 %6.4f; \n',test.MPR3(currentID,1));
    fprintf(f,'set c_mono %6.4f; \n',test.c_mono(currentID,1));
    fprintf(f,'set c_cycl %6.4f; \n',test.c_cycl(currentID,1));
    fprintf(f,'set c_symm %6.4f; \n',test.c_symm(currentID,1));
    fprintf(f,'set k1 %6.4f; \n',test.k1(currentID,1));
    fprintf(f,'set k2 %6.4f; \n',test.k2(currentID,1));
    fprintf(f,'set b1 %6.4f; \n',test.b1(currentID,1));
    fprintf(f,'set b2 %6.4f; \n',test.b2(currentID,1));
    fprintf(f,'set bs %6.4f; \n',test.bs(currentID,1));
    fprintf(f,'set BR %6.4f; \n',test.R(currentID,1));
    fprintf(f,'set barlayout %d; \n',test.barlayout(currentID,1));
    fprintf(f,'set Sy %6.4f; \n',test.syaci(currentID,1));
    fprintf(f,'set Su %6.4f; \n',test.suaci(currentID,1));
    fprintf(f,'set numIntgrPts %d; \n',numIP);
    fprintf(f,'set DFTag 1; \n'); % no fracture
    if currentID == 11
        fprintf(f,'set LWTag 1; \n'); % lightweight concrete
    else
        fprintf(f,'set LWTag 0; \n');
    end
    fprintf(f,'set CR %6.4f; \n',test.CR(currentID,1));
    fprintf(f,'set ER %6.4f; \n',test.ER(currentID,1));
    fclose(f);
    
    %%
    ! OpenSees.exe ColumnCyclicTestSpiral.tcl
    
    %% Collect data
    ex = importdata(fullfile('./CyclicOutputSpiral/','disp.out'));
    simu.disp{currentID,1} = ex(:,2);
    ex = importdata(fullfile('./CyclicOutputSpiral/','force.out'));
    simu.shear{currentID,1} = -ex(:,2);
    % bar-slip and shear drift contribution
    ex = importdata('./CyclicOutputSpiral/BarSlipForce.out');
    simu.BarSlipForce{currentID,1} = -ex(:,2);
    ex = importdata('./CyclicOutputSpiral/BarSlipDisp.out');
    simu.BarSlipDrift{currentID,1} = -ex(:,3);
    simu.ShearDrift{currentID,1} = -ex(:,4);
    % steel strain
    ex = importdata(fullfile('./CyclicOutputSpiral/','SS0.out'));
    simu.stress.steelbot{currentID,1} = ex(:,2);
    simu.strain.steelbot{currentID,1} = ex(:,3);    
    ex = importdata(fullfile('./CyclicOutputSpiral/', ...
        ['SS' num2str(round(test.nsl(currentID,1)/2)) '.out']));
    simu.stress.steeltop{currentID,1} = ex(:,2);
    simu.strain.steeltop{currentID,1} = ex(:,3);
    % FI
    for i = 1:1:test.nsl(currentID,1)
        temp = importdata(['./CyclicOutputSpiral/', 'FI',num2str(i-1),'.out']);
        hist_fi(:,i) = temp(:,2);
        if max(hist_fi(:,i))>=1.0
            simuftag(i) = min(find(hist_fi(:,i)>=1.0));
        else
            simuftag(i) = 0;
        end
        temp = importdata(['./CyclicOutputSpiral/', 'SS',num2str(i-1),'.out']);
        hist_strain(:,i) = temp(:,3);
        hist_stress(:,i) = temp(:,2);
    end
    
    % plot
    figure;
    plot(test.disp{currentID,1}/test.L(currentID,1), ...
        test.shear{currentID,1},'-','color','k','linewidth',0.5);
    hold on;
    plot(simu.disp{currentID,1}/test.L(currentID,1), ...
        simu.shear{currentID,1},'--','color',beige,'linewidth',0.5);
    plot(simu.disp{currentID,1}(simuftag(simuftag>0))/test.L(currentID,1), ...
        simu.shear{currentID,1}(simuftag(simuftag>0)),'o', ...
        'MarkerFaceColor',cardinal);
    xlabel('Lateral drift');
    ylabel('Base shear (kip)');
    grid on;
    legend('Test','Simulation','location','southeast');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_hyst_wf.fig']);
    close(gcf);
    % shear history
    figure;
    plot(v_m,'-k','linewidth',1);
    hold on;
    plot(simu.shear{currentID,1},'--','color',beige,'linewidth',1);
    plot(simuftag(simuftag>0),simu.shear{currentID,1}(simuftag(simuftag>0)),'o', ...
        'MarkerFaceColor',cardinal)
    grid on;
    xlabel('Load step');
    ylabel('Base shear (kip)');
    legend('Test','Simulation','location','best');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_shearhist_wf.fig']);
    close(gcf);
    % Fracture index
    figure;
    plot(hist_fi,'-k');
    grid on;
    xlabel('Load step');
    ylabel('Fracture index history');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_fihist_wf.fig']);
    close(gcf);
    % drift decomposition
    figure;
    area(simu.disp{currentID,1}/test.L(currentID,1), ...
        'FaceColor',beige,'linewidth',0.5)
    hold on;
    area(simu.disp{currentID,1}/test.L(currentID,1)- ...
        simu.ShearDrift{currentID,1},'FaceColor',[0.9,0.9,0.9],'linewidth',0.5)
    area(simu.disp{currentID,1}/test.L(currentID,1)- ...
        simu.ShearDrift{currentID,1}- ...
        simu.BarSlipDrift{currentID,1},'FaceColor',[0.7,0.7,0.7],'linewidth',0.5)
    grid on;
    xlabel('Load step');
    ylabel('Lateral drift');
    legend('Shear','Bar-slip','Flexure');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_driftdecomp_wf.fig']);
    close(gcf);
    % Steel strain-stress
    figure;
    plot(simu.strain.steelbot{currentID,1}, ...
        simu.stress.steelbot{currentID,1},'-k','linewidth',1);
    hold on;
    plot(simu.strain.steeltop{currentID,1}, ...
        simu.stress.steeltop{currentID,1},'-','color',beige,'linewidth',1);
    grid on;
    xlabel('Steel strain');
    ylabel('Steel stress (ksi)');
    legend('Negative-side rebar','Positive-side rebar','location','north');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_steelstrainstress_wf.fig']);
    close(gcf);
    
    clear hist_fi ftag simuftag hist_strain hist_stress;
    clc;
end