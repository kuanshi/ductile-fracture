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
DVInputFileName = 'FractureTestPropertiesRectangular.xlsx';

%% Define design variables
ex = importdata(DVInputFileName);
test.ID = ex.data(:,1);
test.TSN = ex.data(:,4);
test.UnitTag = ex.data(:,5);
test.b = ex.data(:,6);
test.h = ex.data(:,7);
test.c = ex.data(:,8);
test.fc = ex.data(:,9);
test.fyl = ex.data(:,10);
test.TY = ex.data(:,11);
test.Es = ex.data(:,12);
test.Epyr = ex.data(:,13);
test.esh = ex.data(:,14);
test.eult = ex.data(:,15);
test.fyt = ex.data(:,16);
test.nlt = ex.data(:,17);
test.nlb = ex.data(:,18);
test.nlm = ex.data(:,19);
test.nt = ex.data(:,20);
test.s = ex.data(:,21);
test.Ashi = ex.data(:,22);
test.Asli = ex.data(:,23);
test.anaType = ex.data(:,24);
test.axialLoad = ex.data(:,25);
test.maxK = ex.data(:,26);
test.maxD = ex.data(:,27);
test.dirc = ex.data(:,28);
test.Lcol = ex.data(:,29);
test.numEle = ex.data(:,30);
test.numIntgrPts = ex.data(:,31);
test.numIncr = ex.data(:,32);
test.eleType = ex.data(:,33);
test.db = ex.data(:,34);
test.bs = ex.data(:,35);
test.R = ex.data(:,36);
test.phsi = ex.data(:,37);
test.susy = ex.data(:,38);
test.syaci = ex.data(:,39);
test.suaci = ex.data(:,41);
test.numSec = numel(test.ID);

%% Convert metric units to imperical units
for i = 1:1:test.numSec
    if test.UnitTag(i,1) == 0
        mmin = 25.4;
        mpaksi = 6.895;
        knkip = 4.448;
        test.b(i,1) = test.b(i,1)/mmin;
        test.h(i,1) = test.h(i,1)/mmin;
        test.c(i,1) = test.c(i,1)/mmin;
        test.fc(i,1) = test.fc(i,1)/mpaksi;
        test.fyl(i,1) = test.fyl(i,1)/mpaksi;
        test.Es(i,1) = test.Es(i,1)/mpaksi;
        test.fyt(i,1) = test.fyt(i,1)/mpaksi;
        test.s(i,1) = test.s(i,1)/mmin;
        test.Ashi(i,1) = test.Ashi(i,1)/mmin/mmin;
        test.Asli(i,1) = test.Asli(i,1)/mmin/mmin;
        test.axialLoad(i,1) = test.axialLoad(i,1)/knkip;
        test.Lcol(i,1) = test.Lcol(i,1)/mmin;
        test.db(i,1) = test.db(i,1)/mmin;
    end
end

%% Compute other modeling parameters
% elastic concrete modulus
test.Ec = 57*sqrt(-test.fc*1000);
% slip related terms
test.ub = 12*sqrt(-test.fc*1000)/1000;
test.ubp = 6*sqrt(-test.fc*1000)/1000;
test.le = test.h-2*test.c-sqrt(4*test.Ashi/pi);
% compute the equivalent lateral buckleproof capacity
test.kt = test.Es.*test.Ashi./test.le.*(test.nt./test.nlt);
test.kt(3:4,1) = test.kt(3:4,1)/sqrt(2);
% the required spring stiffness for buckle modes
test.keff = 0.04*test.Es.*(pi*test.db.^4/64);
test.kelastic = pi^4*test.keff./(test.s.^3);
test.kratio = test.kt./test.kelastic;
% test.kmodefactor = [0.75 0.1649 0.0976 0.0448 0.0084 0.0063 0.0037 0.0031 0.0013 0.0009];
test.kmodefactor = [0.5 0.1173 0.0377 0.0154 0.0075 0.0040 0.0024 0.0015 0.0010 0.0007];
test.kequ = test.kelastic*test.kmodefactor;
for i = 1:1:length(test.ID)
    test.modenum(i,1) = max(1,min(find(test.kmodefactor<test.kratio(i,1))));
end
test.betabuckle = min([1*ones(numel(test.ID),1), ...
    max([0.75*ones(numel(test.ID),1), ...
    0.9571+0.1209*test.fyl/60+0.1276*test.TY-0.0698*test.s./test.db],[],2)],[],2);
% buckle strain (per Dhakal and Maekawa equation)
test.ebuckle = test.fyl./test.Es.*(55-2.3*sqrt(test.fyl/14.5).*test.modenum.*test.s./test.db);
% effective s/db
test.sdb = test.s./test.db;
test.sdb_eff = test.phsi.*test.sdb;
% deformation properties
Hdb = 0.07;
RminH = 1.46;
% MP curve hardening
test.MPR1 = 3.0*test.fyl./test.fyl; % fix at 3.0
test.MPR2 = exp(3.4439+0.1790*log(test.fyl/60));
test.MPR3 = exp(-1.0389-2.1958*log(test.fyl/60)-1.9579*log(test.TY)- ...
    1.7812*log(test.eult)-0.2175*log(test.db));
% fracture model
test.c_mono = exp(-3.7513-1.7857*log(test.eult)+0.2892*log(test.db));
% test.c_cycl = exp(0.9342+0.0732*log(test.fyl)+0.9444*log(test.db))+0.4984*log(Hdb/0.07);
% test.c_cycl = exp(4.24+1.28*log(test.fyl/60)+ ...
%     1.56*log(test.eult)+0.89*log(test.db));
test.c_cycl = exp(4.64+1.36*log(test.fyl/60)+ ...
    1.74*log(test.eult)+0.89*log(test.db)-0.07*log(RminH/1.46));
test.c_symm = 1.02*test.fyl./test.fyl; % fix at 1.02
test.k1 = exp(2.12-0.66*log(test.db));
test.k2 = exp(1.8900-1.3793*log(test.TY)-0.5085*log(test.db));
test.b1 = exp(-2.43-1.66*log(test.TY)- ...
    0.42*log(RminH/1.46));
test.b2 = exp(-2.83-0.43*log(test.eult)- ...
    0.88*log(test.sdb_eff)- ...
    0.12*log(RminH/1.46));

% define load tag
% 0: loading protocol in tests; 1: FEMA 461; 2: modified FEMA461; 
% 3: ACI 374; 4: modified ACI 374
loadtag = 0;

% define specimen tag
tag1 = 3;
tag2 = 24; % <= test.numSec
for currentID = tag1:1:tag2
    %% Create excutable file of design variables
    f = fopen('DesignVariableRectangular.tcl','wt');
    fprintf(f,'set b %6.4f; \n',test.b(currentID,1));
    fprintf(f,'set h %6.4f; \n',test.h(currentID,1));
    fprintf(f,'set Acol [expr $b*$h]; \n');
    fprintf(f,'set c %6.4f; \n',test.c(currentID,1));
    fprintf(f,'set fc %6.4f; \n',test.fc(currentID,1));
    fprintf(f,'set fyl %6.4f; \n',test.fyl(currentID,1));
    fprintf(f,'set TY %6.4f; \n',test.TY(currentID,1));
    fprintf(f,'set ful [expr $TY*$fyl]; \n');
    fprintf(f,'set Es %6.4f; \n',test.Es(currentID,1));
    fprintf(f,'set esh %6.4f; \n',test.esh(currentID,1));
    fprintf(f,'set eult %6.4f; \n',test.eult(currentID,1));
    fprintf(f,'set steel_b [expr $fyl*($TY-1.0)/$Es/($eult-$fyl/$Es)]; \n');
    fprintf(f,'set Epyr %6.4f; \n',test.Epyr(currentID,1));
    fprintf(f,'set Esh [expr $Epyr*$Es]; \n');
    fprintf(f,'set fyt %6.4f; \n',test.fyt(currentID,1));
    fprintf(f,'set nlt %d; \n',test.nlt(currentID,1));
    fprintf(f,'set nlb %d; \n',test.nlb(currentID,1));
    fprintf(f,'set nlm %d; \n',test.nlm(currentID,1));
    if (currentID == 24) || ((currentID > 16) && (currentID < 20))
        fprintf(f,'set nl [expr $nlt+$nlb]; \n');
    else
        fprintf(f,'set nl [expr $nlt+$nlb+2*$nlm]; \n');
    end
    fprintf(f,'set nt %d; \n',test.nt(currentID,1));
    fprintf(f,'set s %6.4f; \n',test.s(currentID,1));
    fprintf(f,'set Ashi %6.4f; \n',test.Ashi(currentID,1));
    fprintf(f,'set Asli %6.4f; \n',test.Asli(currentID,1));
    fprintf(f,'set Ash [expr $nt*$Ashi]; \n');
    fprintf(f,'set rou [expr $Ash/$b/$s]; \n');
    fprintf(f,'set SecTag %d; \n',test.ID(currentID,1));
    fprintf(f,'set db %d; \n',test.db(currentID,1));
    fprintf(f,'set ub [expr 12*sqrt(-$fc*1000)/1000]; \n');
    fprintf(f,'set theta_sy [expr $fyl/$Es*$fyl*$db/8.0/$ub/($h-2*$c)]; \n');
    fprintf(f,'set theta_su [expr $db/8.0/$ub/($h-2*$c)*($fyl/$Es*$fyl+2*($fyl/$Es+$eult)*($ful-$fyl))]; \n');
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
    fprintf(f,'set SuSy %6.4f; \n',test.susy(currentID,1));
    fprintf(f,'set Sy %6.4f; \n',test.syaci(currentID,1));
    fprintf(f,'set Su %6.4f; \n',test.suaci(currentID,1));
    fprintf(f,'set DFTag 0; \n'); % no fracture
    fclose(f);
    
    %% Create loading file
    switch loadtag
        case 0
            ex = importdata(fullfile('./TestData', ...
                [num2str(test.TSN(currentID,1)) '.txt']));
            if test.UnitTag(currentID,1) == 0
                disp = ex(:,1)/mmin;
            else
                disp = ex(:,1);
            end
            test.disp{currentID,1} = ex(:,1);
            test.shear{currentID,1} = ex(:,2);
        case 1
            disp = [];
            ldg_step = [0; 0.0006; -0.0006; 0.0006; -0.0006; 0.0006; -0.0006;
                0.0061; -0.0061; 0.0061; -0.0061; 0.0061; -0.0061;
                0.0086; -0.0086; 0.0086; -0.0086; 0.0086; -0.0086;
                0.012; -0.012; 0.012; -0.012; 0.012; -0.012;
                0.017; -0.017; 0.017; -0.017; 0.017; -0.017;
                0.024; -0.024; 0.024; -0.024; 0.024; -0.024;
                0.033; -0.033; 0.033; -0.033; 0.033; -0.033;
                0.046; -0.046; 0.046; -0.046; 0.046; -0.046;
                0.065; -0.065; 0.065; -0.065; 0.065; -0.065];
            for substep = 1:1:length(ldg_step)-1
                disp = [disp; test.Lcol(currentID,1)* ...
                    linspace(ldg_step(substep),ldg_step(substep+1),20)'];
            end
            test.disp{currentID,1} = disp;
        case 2
            disp = [];
            ldg_step = [0; 0.04; -0.04; 0.04; -0.04];
            for substep = 1:1:length(ldg_step)-1
                disp = [disp; test.Lcol(currentID,1)* ...
                    linspace(ldg_step(substep),ldg_step(substep+1),100)'];
            end
            test.disp{currentID,1} = disp;
    end
    if currentID > 0
        curpeaktag = 1;
        curzerotag = 1;
        curpeakamp = 0.0001*test.Lcol(currentID);
        ldg_step = [0];
        ldg_steptag = [1];
        disp_m = [];
        v_m = [];
        for temptag = 1:1:length(disp)-1
            if disp(temptag)*disp(temptag+1)<=0 && ...
                    max(abs(disp(curzerotag:temptag)))>curpeakamp
                [~,curpeaktag] = max(abs(disp(curzerotag:temptag)));
                ldg_step = [ldg_step;disp(curzerotag+curpeaktag-1)];
                ldg_steptag = [ldg_steptag; curzerotag+curpeaktag-1];
                curpeakamp = abs(disp(curzerotag+curpeaktag-1))*0.5;
                curzerotag = temptag;
            end
        end
        [~,curpeaktag] = max(abs(disp(curzerotag:end)));
        ldg_step = [ldg_step; disp(curzerotag+curpeaktag-1)];
        ldg_steptag = [ldg_steptag; curzerotag+curpeaktag-1];
        tag = 1;
        while tag < length(ldg_step)
            if ldg_step(tag) == ldg_step(tag+1)
                ldg_step = [ldg_step(1:tag);ldg_step(tag+2:end)];
                ldg_steptag = [ldg_steptag(1:tag);ldg_steptag(tag+2:end)];
            else
                tag = tag+1;
            end
        end
        xtag = 1;
        while xtag < length(ldg_steptag)
            if ldg_steptag(xtag) == ldg_steptag(xtag+1)
                ldg_steptag = [ldg_steptag(1:xtag);ldg_steptag(xtag+2:end)];
                ldg_step = [ldg_step(1:xtag);ldg_step(xtag+2:end)];
            else
                xtag = xtag+1;
            end
        end
        for substep = 1:1:length(ldg_step)-1
            tempvector = linspace(ldg_step(substep),ldg_step(substep+1),500)';
            disp_m = [disp_m; tempvector(1:end-1)];
            v_m = [v_m; interp1(disp([ldg_steptag(substep),ldg_steptag(substep+1)]), ...
                test.shear{currentID,1}([ldg_steptag(substep),ldg_steptag(substep+1)]), ...
                tempvector(1:end-1))];
        end
    else
        disp_m = disp;
    end
    
%     ex = importdata(fullfile('./TestData', ...
%         [num2str(test.TSN(currentID,1)) '.txt']));
%     disp = ex(:,1);
%     % maximum length per step (0.1% of column height)
%     delta_disp = 0.00001*test.Lcol(currentID,1);
%     disp_m = disp(1,1);
%     v = ex(:,2);
%     v_m = v(1,1);
%     for ldtag = 2:1:length(disp)
%         cur_delta = abs(disp(ldtag)-disp(ldtag-1));
%         if cur_delta < delta_disp
%             disp_m = [disp_m;disp(ldtag)];
%             v_m = [v_m;v(ldtag)];
%         else
%             if (disp(ldtag-1)<disp(ldtag))
%                 temp_disp = [disp(ldtag-1)+delta_disp:delta_disp:disp(ldtag)]';
%             else
%                 temp_disp = [disp(ldtag-1)-delta_disp:-delta_disp:disp(ldtag)]';
%             end
%             if temp_disp(end) ~= disp(ldtag)
%                 temp_disp = [temp_disp;disp(ldtag)];
%             end
%             temp_v = interp1([disp(ldtag-1);disp(ldtag)], ...
%                 [v(ldtag-1);v(ldtag)],temp_disp);
%             disp_m = [disp_m;temp_disp];
%             v_m = [v_m;temp_v];
%             clear temp_disp temp_v;
%         end
%     end
%     test.disp{currentID,1} = disp_m;
%     test.shear{currentID,1} = v_m;
%     
    if currentID < 17 && currentID > 12
        test.shear{currentID,1} = test.shear{currentID,1}+ ...
            0.5*test.axialLoad(currentID,1)*test.disp{currentID,1}/test.Lcol(currentID,1);
    end
    
    f = fopen('LoadingParameterRectangular.tcl','wt');
    fprintf(f,'set Lcol %6.4f; \n',test.Lcol(currentID,1));
    fprintf(f,'set numIncr %d; \n',test.numIncr(currentID,1));
    fprintf(f,'set numEle %d; \n',test.numEle(currentID,1));
    fprintf(f,'set numIntgrPts %d; \n',6);
    fprintf(f,'set P %8.6f; \n',test.axialLoad(currentID,1));
    fprintf(f,'set Dincr %8.6f; \n',0.001);
    fprintf(f,'set LoadHistory { \n');
    for ldstep = 1:1:length(disp_m)
        fprintf(f,'%8.4f \n',disp_m(ldstep));
    end
    fprintf(f,'}; \n');
    fclose(f);
    
    %%
    ! OpenSees.exe ColumnCyclicTestRectangular.tcl
    
    %% Collect data
    ex = importdata(fullfile('./CyclicOutputRectangular/','disp.out'));
    if test.UnitTag(currentID,1) == 0
        simu.disp{currentID,1} = ex(:,2)*mmin;
    else
        simu.disp{currentID,1} = ex(:,2);
    end
    ex = importdata(fullfile('./CyclicOutputRectangular/','force.out'));
    if test.UnitTag(currentID,1) == 0
        simu.shear{currentID,1} = -ex(:,2)*knkip;
    else
        simu.shear{currentID,1} = -ex(:,2);
    end
    % concrete strain
    ex = importdata(fullfile('./CyclicOutputRectangular/','ConcrCoreBot.out'));
    simu.stress.concrcorebot{currentID,1} = ex(:,2);
    simu.strain.concrcorebot{currentID,1} = ex(:,3);
    ex = importdata(fullfile('./CyclicOutputRectangular/','ConcrCoreTop.out'));
    simu.stress.concrcoretop{currentID,1} = ex(:,2);
    simu.strain.concrcoretop{currentID,1} = ex(:,3);
    ex = importdata(fullfile('./CyclicOutputRectangular/','ConcrCoverBot.out'));
    simu.stress.concrcoverbot{currentID,1} = ex(:,2);
    simu.strain.concrcoverbot{currentID,1} = ex(:,3);
    ex = importdata(fullfile('./CyclicOutputRectangular/','ConcrCoverTop.out'));
    simu.stress.concrcovertop{currentID,1} = ex(:,2);
    simu.strain.concrcovertop{currentID,1} = ex(:,3);
    % steel strain
    ex = importdata(fullfile('./CyclicOutputRectangular/','SteelBot.out'));
    simu.stress.steelbot{currentID,1} = ex(:,2);
    simu.strain.steelbot{currentID,1} = ex(:,3);    
    ex = importdata(fullfile('./CyclicOutputRectangular/','SteelTop.out'));
    simu.stress.steeltop{currentID,1} = ex(:,2);
    simu.strain.steeltop{currentID,1} = ex(:,3);
    ex = importdata(fullfile('./CyclicOutputRectangular/','SteelBot2.out'));
    simu.stress.steelbot2{currentID,1} = ex(:,2);
    simu.strain.steelbot2{currentID,1} = ex(:,3);    
    ex = importdata(fullfile('./CyclicOutputRectangular/','SteelTop2.out'));
    simu.stress.steeltop2{currentID,1} = ex(:,2);
    simu.strain.steeltop2{currentID,1} = ex(:,3);
    % bar-slip and shear drift contribution
    ex = importdata('./CyclicOutputRectangular/BarSlipForce.out');
    simu.BarSlipForce{currentID,1} = -ex(:,2);
    ex = importdata('./CyclicOutputRectangular/BarSlipDisp.out');
    simu.BarSlipDrift{currentID,1} = -ex(:,3);
    simu.ShearDrift{currentID,1} = -ex(:,4);
    % compute fracture index
    es_hist = [simu.strain.steelbot{currentID,1}, ...
        simu.strain.steeltop{currentID,1}];
    ss_hist = [simu.stress.steelbot{currentID,1}, ...
        simu.stress.steeltop{currentID,1}];
    num_rlz = 200;
    [FI,rlz] = DuctileFractureModel(es_hist,ss_hist, ...
        test.c_mono(currentID,1),test.c_symm(currentID,1), ...
        test.c_cycl(currentID,1),test.Es(currentID,1), ...
        test.eult(currentID,1),test.k1(currentID,1),test.k2(currentID,1), ...
        test.db(currentID,1),test.b1(currentID,1),test.b2(currentID,1),num_rlz,0);   
    % fracture probability
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
    pf_1 = logncdf(FI_1,0,1.6);
    FI_2 = FI{2,1}(:,1);
    for mtag = 1:1:length(FI_2)
        FI_2(mtag,1) = max(FI_2(1:mtag));
    end
    pf_2 = logncdf(FI_2,0,1.6);
    
    % plot
    figure;
    plot(test.disp{currentID,1}/test.Lcol(currentID,1), ...
        test.shear{currentID,1},'-','color','k','linewidth',1);
    hold on;
    plot(simu.disp{currentID,1}/test.Lcol(currentID,1), ...
        simu.shear{currentID,1},'--','color',beige,'linewidth',1);
    xlabel('Lateral drift');
    ylabel('Base shear (kip)');
    ylim([-100 100]);
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
    ylim([-100 100]);
    grid on;
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_hystbarslip.fig']);
    close(gcf);
    % Shear hysteresis
    figure;
    plot(simu.ShearDrift{currentID,1}, ...
        simu.shear{currentID,1},'-k','linewidth',1);
    xlabel('Shear drift');
    ylabel('Base shear (kip)');
    ylim([-100 100]);
    grid on;
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_hystshear.fig']);
    close(gcf);
    % Drift vs. strain
    figure;
    plot(simu.disp{currentID,1}/test.Lcol(currentID,1), ...
        simu.strain.steelbot{currentID,1},'-k','linewidth',1);
    hold on;
    plot(simu.disp{currentID,1}/test.Lcol(currentID,1), ...
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
    plot(disp_m/test.Lcol(currentID,1),'-k');
    hold on;
    plot(simu.disp{currentID,1}/test.Lcol(currentID,1),'--', ...
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
    plot(disp_m/test.Lcol(currentID,1),'-k');
    hold on;
    plot(simu.disp{currentID,1}/test.Lcol(currentID,1),'--', ...
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
    area(simu.disp{currentID,1}/test.Lcol(currentID,1), ...
        'FaceColor',beige,'linewidth',0.5)
    hold on;
    area(simu.disp{currentID,1}/test.Lcol(currentID,1)- ...
        simu.ShearDrift{currentID,1},'FaceColor',[0.9,0.9,0.9],'linewidth',0.5)
    area(simu.disp{currentID,1}/test.Lcol(currentID,1)- ...
        simu.ShearDrift{currentID,1}- ...
        simu.BarSlipDrift{currentID,1},'FaceColor',[0.7,0.7,0.7],'linewidth',0.5)
    grid on;
    xlabel('Load step');
    ylabel('Lateral drift');
    legend('Shear','Bar-slip','Flexure');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_driftdecomp.fig']);
    close(gcf);
    
    %% Simulate fracture (median model)
    f = fopen('DesignVariableRectangular.tcl','wt');
    fprintf(f,'set b %6.4f; \n',test.b(currentID,1));
    fprintf(f,'set h %6.4f; \n',test.h(currentID,1));
    fprintf(f,'set Acol [expr $b*$h]; \n');
    fprintf(f,'set c %6.4f; \n',test.c(currentID,1));
    fprintf(f,'set fc %6.4f; \n',test.fc(currentID,1));
    fprintf(f,'set fyl %6.4f; \n',test.fyl(currentID,1));
    fprintf(f,'set TY %6.4f; \n',test.TY(currentID,1));
    fprintf(f,'set ful [expr $TY*$fyl]; \n');
    fprintf(f,'set Es %6.4f; \n',test.Es(currentID,1));
    fprintf(f,'set esh %6.4f; \n',test.esh(currentID,1));
    fprintf(f,'set eult %6.4f; \n',test.eult(currentID,1));
    fprintf(f,'set steel_b [expr $fyl*($TY-1.0)/$Es/($eult-$fyl/$Es)]; \n');
    fprintf(f,'set Epyr %6.4f; \n',test.Epyr(currentID,1));
    fprintf(f,'set Esh [expr $Epyr*$Es]; \n');
    fprintf(f,'set fyt %6.4f; \n',test.fyt(currentID,1));
    fprintf(f,'set nlt %d; \n',test.nlt(currentID,1));
    fprintf(f,'set nlb %d; \n',test.nlb(currentID,1));
    fprintf(f,'set nlm %d; \n',test.nlm(currentID,1));
    if (currentID == 24) || ((currentID > 16) && (currentID < 20))
        fprintf(f,'set nl [expr $nlt+$nlb]; \n');
    else
        fprintf(f,'set nl [expr $nlt+$nlb+2*$nlm]; \n');
    end
    fprintf(f,'set nt %d; \n',test.nt(currentID,1));
    fprintf(f,'set s %6.4f; \n',test.s(currentID,1));
    fprintf(f,'set Ashi %6.4f; \n',test.Ashi(currentID,1));
    fprintf(f,'set Asli %6.4f; \n',test.Asli(currentID,1));
    fprintf(f,'set Ash [expr $nt*$Ashi]; \n');
    fprintf(f,'set rou [expr $Ash/$b/$s]; \n');
    fprintf(f,'set SecTag %d; \n',test.ID(currentID,1));
    fprintf(f,'set db %d; \n',test.db(currentID,1));
    fprintf(f,'set ub [expr 12*sqrt(-$fc*1000)/1000]; \n');
    fprintf(f,'set theta_sy [expr $fyl/$Es*$fyl*$db/8.0/$ub/($h-2*$c)]; \n');
    fprintf(f,'set theta_su [expr $db/8.0/$ub/($h-2*$c)*($fyl/$Es*$fyl+2*($fyl/$Es+$eult)*($ful-$fyl))]; \n');
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
    fprintf(f,'set Sy %6.4f; \n',test.syaci(currentID,1));
    fprintf(f,'set Su %6.4f; \n',test.suaci(currentID,1));
    fprintf(f,'set SuSy %6.4f; \n',test.susy(currentID,1));
    fprintf(f,'set DFTag 1; \n');
    fclose(f);
    
    %%
    ! OpenSees.exe ColumnCyclicTestRectangular.tcl
    
    %% Collect data
    ex = importdata(fullfile('./CyclicOutputRectangular/','disp.out'));
    if test.UnitTag(currentID,1) == 0
        simu.disp{currentID,1} = ex(:,2)*mmin;
    else
        simu.disp{currentID,1} = ex(:,2);
    end
    ex = importdata(fullfile('./CyclicOutputRectangular/','force.out'));
    if test.UnitTag(currentID,1) == 0
        simu.shear{currentID,1} = -ex(:,2)*knkip;
    else
        simu.shear{currentID,1} = -ex(:,2);
    end
    % bar-slip and shear drift contribution
    ex = importdata('./CyclicOutputRectangular/BarSlipForce.out');
    simu.BarSlipForce{currentID,1} = -ex(:,2);
    ex = importdata('./CyclicOutputRectangular/BarSlipDisp.out');
    simu.BarSlipDrift{currentID,1} = -ex(:,3);
    simu.ShearDrift{currentID,1} = -ex(:,4);
    % FI
    FI1 = importdata('./CyclicOutputRectangular/FISteelBot.out');
    FI2 = importdata('./CyclicOutputRectangular/FISteelTop.out');
    FItag1 = min(find(FI1(:,2)>=1));
    FItag2 = min(find(FI2(:,2)>=1));
    simu.FI1{currentID,1} = FI1;
    simu.FI2{currentID,1} = FI2;
    ex = importdata('./CyclicOutputRectangular/SteelBot.out');
    simu.strainf.steelbot{currentID,1} = ex(:,3);
    simu.stressf.steelbot{currentID,1} = ex(:,2);
    ex = importdata('./CyclicOutputRectangular/SteelTop.out');
    simu.strainf.steeltop{currentID,1} = ex(:,3);
    simu.stressf.steeltop{currentID,1} = ex(:,2);
    
    % plot
    figure;
    plot(test.disp{currentID,1}/test.Lcol(currentID,1), ...
        test.shear{currentID,1},'-','color','k','linewidth',0.5);
    hold on;
    plot(simu.disp{currentID,1}/test.Lcol(currentID,1), ...
        simu.shear{currentID,1},'--','color',beige,'linewidth',0.5);
    if ~isempty(FItag1)
        plot(simu.disp{currentID,1}(FItag1)/test.Lcol(currentID,1), ...
            simu.shear{currentID,1}(FItag1,1),'o', ...
            'MarkerFaceColor',cardinal);
    end
    if ~isempty(FItag2)
        plot(simu.disp{currentID,1}(FItag2)/test.Lcol(currentID,1), ...
            simu.shear{currentID,1}(FItag2,1),'o', ...
            'MarkerFaceColor',cardinal);
    end
    xlabel('Lateral drift');
    ylabel('Base shear (kip)');
    ylim([-100 100]);
    grid on;
    legend('Test','Simulation','location','southeast');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_hyst_wf.fig']);
    close(gcf);
    % Shear history
    figure;
    plot(v_m,'-k','linewidth',1);
    hold on;
    plot(simu.shear{currentID,1},'--','color',beige,'linewidth',1);
    grid on;
    if ~isempty(FItag1)
        plot(FItag1,simu.shear{currentID,1}(FItag1,1),'o', ...
            'MarkerFaceColor',cardinal);
    end
    if ~isempty(FItag2)
        plot(FItag2,simu.shear{currentID,1}(FItag2,1),'o', ...
            'MarkerFaceColor',cardinal);
    end  
    xlabel('Load step');
    ylabel('Base shear (kip)');
    ylim([-100 100]);
    legend('Test','Simulation','location','best');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_shearhist_wf.fig']);
    close(gcf);
    % Fracture index history
    figure;
    plot(FI1(:,2),'-k','linewidth',1);
    hold on;
    plot(FI2(:,2),'-','color',beige,'linewidth',1);
    if ~isempty(FItag1)
        plot(FItag1,FI1(FItag1,2),'o','MarkerFaceColor',cardinal);
    end
    if ~isempty(FItag2)
        plot(FItag2,FI2(FItag2,2),'o','MarkerFaceColor',cardinal);
    end  
    grid on;
    xlabel('Load step');
    ylabel('Fracture index');
    ylim([0 1.2]);
    legend('Negative-side rebar','Positive-side rebar','location','northwest');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_fihist_wf.fig']);
    close(gcf);
    % drift decomposition
    figure;
    area(simu.disp{currentID,1}/test.Lcol(currentID,1), ...
        'FaceColor',beige,'linewidth',0.5);
    hold on;
    area(simu.disp{currentID,1}/test.Lcol(currentID,1)- ...
        simu.ShearDrift{currentID,1},'FaceColor',[0.9,0.9,0.9],'linewidth',0.5)
    area(simu.disp{currentID,1}/test.Lcol(currentID,1)- ...
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
    plot(simu.strainf.steelbot{currentID,1}, ...
        simu.stressf.steelbot{currentID,1},'-k','linewidth',1);
    hold on;
    plot(simu.strainf.steeltop{currentID,1}, ...
        simu.stressf.steeltop{currentID,1},'-','color',beige,'linewidth',1);
    grid on;
    xlabel('Steel strain');
    ylabel('Steel stress (ksi)');
    legend('Negative-side rebar','Positive-side rebar','location','north');
    saveas(gcf,['./FigOut/V' num2str(test.TSN(currentID,1)) '_steelstrainstress_wf.fig']);
    close(gcf);
    clc;
end