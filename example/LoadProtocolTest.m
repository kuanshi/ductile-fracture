clear;
clc;

test.eult = 0.09;
test.db = 1.0;
test.fyl = 66.0;
test.ful = 92.4;
test.TY = test.ful./test.fyl;
RminH = 1.46;
test.sdb = 1;
test.Es = 29000;

ModelCoefficient
test.c_mono = exp(mc(2,:)+mc(3,:)*log(test.eult)+mc(4,:)*log(test.db));
test.c_cycl = exp(mc(5,:)+mc(6,:)*log(test.fyl/60)+ ...
    mc(7,:)*log(test.eult)+mc(8,:)*log(test.db));
test.c_symm = mc(1,:)*test.fyl./test.fyl;
test.k1 = exp(2.21-0.32*log(test.TY)-0.66*log(test.db));
test.k2 = exp(1.29+0.64*log(test.fyl/60)-0.46*log(test.db));
test.b1 = exp(-2.53-1.90*log(test.TY)- ...
    1.36*log(test.db)-0.43*log(RminH/1.46));
test.b2 = exp(-3.29-0.49*log(test.eult)- ...
    0.70*log(test.sdb)-0.14*log(RminH/1.46));

%% Monotonic loads
monodata = importdata('Gr60MaterialMono.mat');
es_hist = monodata.es';
ss_hist = monodata.fs';
for curtag = 1:1:length(test.c_cycl)
    [FI,~] = DuctileFractureModel_newBAM(es_hist,ss_hist, ...
        test.c_mono(curtag),test.c_symm(curtag), ...
        test.c_cycl(curtag),test.Es,test.eult,test.k1,test.k2, ...
        test.db,test.b1,test.b2,1,1);
    FI_mono(:,curtag) = FI{1,1};
end
clear FI;

figure;
plot(FI_mono);
grid on;

%% Cyclic loads
cycldata = importdata('Gr60MaterialCycl.mat');
es_hist = cycldata.es*0.8;
ss_hist = cycldata.fs;
for curtag = 1:1:length(test.c_cycl)
    [FI,~] = DuctileFractureModel_newBAM(es_hist,ss_hist, ...
        test.c_mono(curtag),test.c_symm(curtag), ...
        test.c_cycl(curtag),test.Es,test.eult,test.k1,test.k2, ...
        test.db,test.b1,test.b2,1,1);
    FI_cycl(:,curtag) = FI{1,1};
end
clear FI;

FIm_cycl = max(FI_cycl);
figure;
plot(FI_cycl);
grid on;

%% Random loads
randdata = importdata('Gr60MaterialRand.mat');
es_hist = randdata.es;
ss_hist = randdata.fs;
for curtag = 1:1:length(test.c_cycl)
    [FI,~] = DuctileFractureModel_newBAM(es_hist,ss_hist, ...
        test.c_mono(curtag),test.c_symm(curtag), ...
        test.c_cycl(curtag),test.Es,test.eult,test.k1,test.k2, ...
        test.db,test.b1,test.b2,1,1);
    FI_rand(:,curtag) = FI{1,1};
end
clear FI;

FIm_rand = max(FI_rand);
figure;
plot(FI_rand);
grid on;