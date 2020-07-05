function [FI,rlz] = DuctileFractureModel(es_hist_tot,ss_hist_tot, ...
    c_mono,c_symm,c_cycl,Es,esu,k1,k2,db,b1,b2,num_rlz,ceiltag)
    if nargin == 13
        ceiltag = 1;
    end
    % mean and covariance matrices in DFM
    mu_1 = log([c_mono;c_cycl]);
    mu_2 = log([k1;k2]);
    mu_3 = log([b1;b2]);
    sigma_1 = [0.0361 -0.0171; -0.0171 0.1296]; % c_mono & c_cycl
    sigma_2 = [0.023 0.025; 0.025 0.044]; % k1 & k2
    sigma_3 = [0.0953 -0.0874; -0.0874 0.0881]; % b1 & b2
    % generate realizations (include the median)
    temprlz = mvnrnd(mu_1,sigma_1,num_rlz-1);
    rlz.c_mono = [c_mono;exp(temprlz(:,1))];
    rlz.c_cycl = [c_cycl;exp(temprlz(:,2))];
    temprlz = mvnrnd(mu_2,sigma_2,num_rlz-1);
    rlz.k1 = [k1;exp(temprlz(:,1))];
    rlz.k2 = [k2;exp(temprlz(:,2))];
    temprlz = mvnrnd(mu_3,sigma_3,num_rlz-1);
    rlz.b1 = [b1;exp(temprlz(:,1))];
    rlz.b2 = [b2;exp(temprlz(:,2))];
%     rlz.k1 = k1*ones(num_rlz,1);
%     rlz.k2 = k2*ones(num_rlz,1);
%     rlz.b1 = b1*ones(num_rlz,1);
%     rlz.b2 = b2*ones(num_rlz,1);
    % compute fracture index
    for tag_hist = 1:1:length(es_hist_tot(1,:))
        es_gage = es_hist_tot(:,tag_hist);
        ss_gage = ss_hist_tot(:,tag_hist);
        FI_temp = [];
        for tag_rlz = 1:1:num_rlz
            cur_c_mono = rlz.c_mono(tag_rlz);
            cur_c_cycl = rlz.c_cycl(tag_rlz);
            cur_k1 = rlz.k1(tag_rlz);
            cur_k2 = rlz.k2(tag_rlz);
            cur_b1 = rlz.b1(tag_rlz);
            cur_b2 = rlz.b2(tag_rlz);
            % compute the reference strain history
            es_max = [];
            for i = 1:1:length(es_gage)
                es_max(i,1) = max(es_gage(1:i));
            end
            % compute strain memory factor
            es_min = [];
            for i = 1:1:length(es_gage)
                es_min(i,1) = min(es_gage(1:i));
            end
            e_memo = min([ones(length(es_min),1),(es_max-es_min)/0.05],[],2);
            % necking amplification model
            for i = 1:1:length(es_gage)
                if es_gage(i,1) > esu
                    es_hist(i,1) = esu+cur_k1*(es_gage(i,1)-esu);
                    T(i,1) = 0.33+cur_k2*(es_gage(i,1)-esu);
                else
                    es_hist(i,1) = es_gage(i,1);
                    T(i,1) = 0.33;
                end
            end
            % buckle adjustment model
            phi = cur_b1*sinh((es_max-es_gage)/cur_b2);
            es_hist = -phi*db/2+es_hist;
            % compute damage indices
            for i = 1:1:length(es_hist)
                if abs(es_hist(i,1)) < abs(ss_gage(i,1)/Es)
                    ep_hist(i,1) = 0;
                else
                    ep_hist(i,1) = es_hist(i,1)-ss_gage(i,1)/Es;
                end
            end
            dep = diff(ep_hist);
            tag_nega = dep<0;
            dep_nega = tag_nega.*abs(dep);
            cep_nega = 0;
            dDI = cur_c_mono.*(((c_symm-1.0)*e_memo(2:end)+1.0).*exp(1.3*sign(dep).*T(2:end))- ...
                exp(-1.3*sign(dep).*T(2:end))).*abs(dep);
            DI_vgm = 0;
            for k = 1:1:length(dep_nega)
                DI_vgm(k+1,1) = max(0,DI_vgm(k,1)+dDI(k));
                cep_nega(k+1,1) = cep_nega(k,1)+dep_nega(k);
            end
            FI_temp(:,tag_rlz) = exp(cur_c_cycl*e_memo.*cep_nega).*DI_vgm;
            if max(FI_temp(:,tag_rlz)) >= 1.0 && ceiltag == 1
                FI_temp(min(find(FI_temp(:,tag_rlz)>=1.0)):end,tag_rlz) = 1.0;
            end
            % perturbation for median model
            if tag_rlz == 1
                % c_mono
                DI_vgm = 0;
                cep_nega = 0;
                ptb_c_mono = cur_c_mono*1.001;
                ptb_dDI = ptb_c_mono*(((c_symm-1.0)*e_memo(2:end)+1.0).*exp(1.3*sign(dep).*T(2:end))- ...
                    exp(-1.3*sign(dep).*T(2:end))).*abs(dep);
                for k = 1:1:length(dep_nega)
                    DI_vgm(k+1,1) = max(0,DI_vgm(k,1)+ptb_dDI(k));
                    cep_nega(k+1,1) = cep_nega(k,1)+dep_nega(k);
                end
                ptb_FI = exp(cur_c_cycl*e_memo.*cep_nega).*DI_vgm;
                rlz.sFI.c_mono{tag_hist,1} = (ptb_FI-FI_temp(:,tag_rlz))/ ...
                    cur_c_mono/0.001*0.19;
                % c_cycl
                DI_vgm = 0;
                cep_nega = 0;
                ptb_c_cycl = cur_c_cycl*1.001;
                dDI = cur_c_mono*(((c_symm-1.0)*e_memo(2:end)+1.0).*exp(1.3*sign(dep).*T(2:end))- ...
                exp(-1.3*sign(dep).*T(2:end))).*abs(dep);
                for k = 1:1:length(dep_nega)
                    DI_vgm(k+1,1) = max(0,DI_vgm(k,1)+dDI(k));
                    cep_nega(k+1,1) = cep_nega(k,1)+dep_nega(k);
                end
                ptb_FI = exp(ptb_c_cycl*e_memo.*cep_nega).*DI_vgm;
                rlz.sFI.c_cycl{tag_hist,1} = (ptb_FI-FI_temp(:,tag_rlz))/ ...
                    cur_c_cycl/0.001*0.36;
                % b1
                DI_vgm = 0;
                cep_nega = 0;
                ptb_b1 = cur_b1*1.001;
                % necking amplification model
                for i = 1:1:length(es_gage)
                    if es_gage(i,1) > esu
                        es_hist(i,1) = esu+cur_k1*(es_gage(i,1)-esu);
                        T(i,1) = 0.33+cur_k2*(es_gage(i,1)-esu);
                    else
                        es_hist(i,1) = es_gage(i,1);
                        T(i,1) = 0.33;
                    end
                end
                % buckle adjustment model
                phi = ptb_b1*sinh((es_max-es_gage)/cur_b2);
                es_hist = -phi*db/2+es_hist;
                % compute damage indices
                for i = 1:1:length(es_hist)
                    if abs(es_hist(i,1)) < abs(ss_gage(i,1)/Es)
                        ep_hist(i,1) = 0;
                    else
                        ep_hist(i,1) = es_hist(i,1)-ss_gage(i,1)/Es;
                    end
                end
                dep = diff(ep_hist);
                tag_nega = dep<0;
                dep_nega = tag_nega.*abs(dep);
                cep_nega = 0;
                dDI = cur_c_mono*(((c_symm-1.0)*e_memo(2:end)+1.0).*exp(1.3*sign(dep).*T(2:end))- ...
                    exp(-1.3*sign(dep).*T(2:end))).*abs(dep);
                DI_vgm = 0;
                for k = 1:1:length(dep_nega)
                    DI_vgm(k+1,1) = max(0,DI_vgm(k,1)+dDI(k));
                    cep_nega(k+1,1) = cep_nega(k,1)+dep_nega(k);
                end
                ptb_FI = exp(cur_c_cycl*e_memo.*cep_nega).*DI_vgm;
                rlz.sFI.b1{tag_hist,1} = (ptb_FI-FI_temp(:,tag_rlz))/ ...
                    cur_b1/0.001*0.31;
            end
            display(['Complete ',num2str(((tag_hist-1)* ...
                num_rlz+tag_rlz)/num_rlz/ ...
                length(es_hist_tot(1,:))*100),'%']);
        end
        FI{tag_hist,1} = FI_temp;
    end
end