pdb_array= {'q20_monomer' 'Abeta420' 'Abeta42_3' 'Abeta40_3' 'Abeta42_40_3'};
%pdb_array={'IkBa'};
Tf0=430;
Tf1=410;
id_protein1=1; 
id_protein2=1;
labels='q20';
Tmin=410; Tmax=450; dT=2;
qo_flag =0;
curve_shift_flag=1; q0_shift=0.55;

curve_color={'m' 'g' 'b' 'm' 'g' 'b' };
n_T=(Tmax-Tmin)/dT; T_plot = Tmin:dT:(Tmax-dT);     

for protein_i=id_protein1:id_protein2
    fsize=20; tsize=16; mr=2; mc=1; 
    label = labels(protein_i);             pdbID_upper = pdb_array{protein_i};    
    
    path = sprintf('/Users/mingchenchen/Documents/Q40_Qbias/Q20_hexamer_Q400_430K_QW_separatestart_norg_500_nomem_pbc_rgoff_notseparate_allbeta');    

    scrnsize = get(0,'ScreenSize'); figure('position', [1 scrnsize(4) 0.25*scrnsize(3) scrnsize(4)]);
    figure(1); grid on, hold on, set(gca, 'fontsize', fsize); xlabel('T (K)', 'fontsize', fsize); ylabel('Cv', 'fontsize', fsize); title([pdbID_upper, ' sim ', num2str(label(1))], 'fontsize', fsize);
    for i_label=1:length(label)
        sim_label = label(i_label);
        %load q and qo
        filename = sprintf('%s/p_total',path); q = load(filename);
        if qo_flag == 1
            filename = sprintf('%s/qo_%d',path, sim_label); qo = load(filename); 
        end
        Nsample = length(q);

        %load energy
        filename = sprintf('%s/e_total',path); E = load(filename);
        cv=zeros(n_T,1); cm=colormap(jet(n_T));
        for i_T=1:n_T
            %load pmf file and calculate pi_sample for cv calculation
            T_i=Tmin+(i_T-1)*(Tmax-Tmin)/n_T;
            %filename=sprintf('%s/%d_pmf.dat',path, sim_label);
            filename=sprintf('%s/q20_%d_pmf.dat',path, T_i);        
            FF=load(filename); qx=FF(:,1);  Fy = FF(:,2); nbin=length(qx);
            dq=qx(2)-qx(1); qmin=qx(1)-dq/2; qmax= qx(nbin)+dq/2;
            Py=exp(-Fy/(0.001987*T_i)); P_norm = sum(Py); Py=Py/P_norm;
            pi_sample = zeros(Nsample,1); ni_sample = zeros(nbin, 1);
            %calculate pi_sample
            for i_bin= 1:nbin
                qi_min = qmin + (i_bin-1)*dq; qi_max= qi_min + dq;
                ids = find( q >= qi_min & q < qi_max ) ;    
                ni_sample(i_bin) = length(ids);        
                if ni_sample(i_bin) > 0
                    pi_sample(ids) = Py(i_bin)/ni_sample(i_bin);
                end
            end    
            %fprintf('probability = %.3f\n', sum(pi_sample));        
            %Calculate <E^2> and (<E>)^2
            E_sq_mean = 0; E_mean = 0;
            for i_sample=1:Nsample
                E_i = E(i_sample,1);
                E_sq_mean = E_sq_mean + E_i^2*pi_sample(i_sample);
                E_mean = E_mean + E_i*pi_sample(i_sample);
            end
            cv(i_T) = ( E_sq_mean - E_mean^2 ) / 0.001987/ T_i^2;	
        end
        T_list = Tmin:dT:(Tmax-dT);
        cv_max = max(cv); id = find(cv==cv_max); Tf=Tf0;%Tf =T_list(id);
        %Tf = 370;
        fprintf('Tf=%f\n', Tf); T_list_normalized = T_list./Tf ;
        plot(T_list, cv, 'k', 'linewidth', 4, 'Marker', 'o');
        if curve_shift_flag==0
            Fmin = min(FF(:,2)); id_shift = find( FF(:,2) == Fmin );    
        else
            [~,id_shift]=min(abs(FF(:,1)-q0_shift));
        end
        %plot F Vs. Qo, shift free energy curve
        %filename=sprintf('%s/%d_%d_pmf.dat',path, sim_label, Tf); 
        %FF=load(filename); 
        scrnsize = get(0,'ScreenSize'); figure('position', [1 scrnsize(4) 0.25*scrnsize(3) scrnsize(4)]);
        figure(2); grid on, hold on, set(gca,'fontsize', fsize), ylabel('F (kcal/mol)')  , legend(['Tmin ', num2str(Tmin)], ['Tmax ', num2str(Tmax-dT)]);
        
        %xlim([0 1]);
        for i_plot = 1:length(T_plot)
            T_i_plot = T_plot(i_plot) ;
            filename=sprintf('%s/q20_%d_pmf.dat',path, T_i_plot); 
            FF=load(filename); 
            if qo_flag==0
                plot(FF(:,1), FF(:,2)-FF(id_shift(1),2) , 'color', cm(i_plot,:), 'linewidth', 3); 
            end
            if qo_flag==1                            
                qx=FF(:,1);  Fy = FF(:,2); nbin=length(qx);
                dq=qx(2)-qx(1); qmin=qx(1)-dq/2; %qmax= qx(nbin)+dq/2;
                Py=exp(-Fy/(0.001987*T_i)); P_norm = sum(Py); Py=Py/P_norm;
                pi_sample = zeros(Nsample,1); ni_sample = zeros(nbin, 1);
                %calculate pi_sample
                for i_bin= 1:nbin
                    qi_min = qmin + (i_bin-1)*dq; qi_max= qi_min + dq;
                    ids = find( q >= qi_min & q < qi_max ) ;    
                    ni_sample(i_bin) = length(ids);        
                    if ni_sample(i_bin) > 0
                        pi_sample(ids) = Py(i_bin)/ni_sample(i_bin);
                    end
                end    
                q_plot=linspace(0.1,0.9,nbin); dq=(0.9-0.1)/nbin; F_qo=ones(nbin,1);
                sum_p=0;
                for i_bin=1:nbin
                    qi_min = 0.1 + (i_bin-1)*dq; qi_max= qi_min + dq;
                    ids = find( qo >= qi_min & qo < qi_max ) ;
                    pi=sum(pi_sample(ids)); sum_p=sum_p + pi;
                    F_qo(i_bin)=-log(pi)*0.001987*T_i_plot;
                end
                fprintf('check total probability: sum_p=%d\n', sum_p);
                plot(q_plot, F_qo, 'color', cm(i_plot,:), 'linewidth', 3); 
            end
        end
        if qo_flag == 1
               xlabel('Q_{O}')
        else
               %xlabel('Q_{w}')
               xlabel('Q')
        end
    end
end
%}
legend('250K','260K','270K','280K','290K','300K','310K','320K','330K','340K','350K','360K','370K','380K','390K');title('ccQL Free energy plot')
scrnsize = get(0,'ScreenSize'); figure('position', [1 scrnsize(4) 0.25*scrnsize(3) scrnsize(4)]);
figure(3), hold on,set(gca,'fontsize', fsize), ylabel('F (kcal/mol)')
for  protein_i=id_protein1:id_protein2
        Tf=Tf0;
        %label = labels(protein_i);
        label = labels;
        pdbID_upper = pdb_array{protein_i};
        path = sprintf('/Users/mingchenchen/Documents/Q40_Qbias/Q20_hexamer_Q400_430K_QW_separatestart_norg_500_nomem_pbc_rgoff_notseparate_allbeta');
        %path = sprintf('~/work/mdbox/qbias/%s',      pdbID_upper);    
        filename=sprintf('%s/%s_%d_pmf.dat',path, label, Tf); FF=load(filename);
        if curve_shift_flag==0
            Fmin = min(FF(:,2)); id_shift = find( FF(:,2) == Fmin );    
        else
            [~,id_shift]=min(abs(FF(:,1)-q0_shift));
        end
        if qo_flag==0
            plot(FF(:,1), FF(:,2)-FF(id_shift(1),2), 'linewidth', 3, 'color',curve_color{protein_i}); grid on;
        end
        if qo_flag==1
                qx=FF(:,1);  Fy = FF(:,2); nbin=length(qx);
                dq=qx(2)-qx(1); qmin=qx(1)-dq/2; qmax= qx(nbin)+dq/2;
                Py=exp(-Fy/(0.001987*T_i)); P_norm = sum(Py); Py=Py/P_norm;
                pi_sample = zeros(Nsample,1); ni_sample = zeros(nbin, 1);
                %calculate pi_sample
                for i_bin= 1:nbin
                    qi_min = qmin + (i_bin-1)*dq; qi_max= qi_min + dq;
                    ids = find( q >= qi_min & q < qi_max ) ;    
                    ni_sample(i_bin) = length(ids);        
                    if ni_sample(i_bin) > 0
                        pi_sample(ids) = Py(i_bin)/ni_sample(i_bin);
                    end
                end    
                q_plot=linspace(0.1,0.9,nbin); dq=(0.9-0.1)/nbin; F_qo=ones(nbin,1);
                sum_p=0;
                for i_bin=1:nbin
                    qi_min = 0.1 + (i_bin-1)*dq; qi_max= qi_min + dq;
                    ids = find( qo >= qi_min & qo < qi_max ) ;
                    pi=sum(pi_sample(ids)); sum_p=sum_p + pi;
                    F_qo(i_bin)=-log(pi)*0.001987*T_i_plot;
                end
                %fprintf('check total probability: sum_p=%d\n', sum_p);
                plot(q_plot, F_qo, 'color', curve_color{protein_i}, 'linewidth', 3);
        end
        title(['T_f =', num2str(Tf),'K'],'fontsize', fsize);   
end


legend('250K','260K','270K','280K','290K','300K','310K','320K','330K','340K','350K','360K','370K','380K','390K');%title('ccQL Free energy plot 390K')
scrnsize = get(0,'ScreenSize'); figure('position', [1 scrnsize(4) 0.25*scrnsize(3) scrnsize(4)]);
figure(4), hold on,set(gca,'fontsize', fsize), ylabel('F (kcal/mol)')
for  protein_i=id_protein1:id_protein2
        Tf=Tf0;
        %label = labels(protein_i);
        label = labels;
        pdbID_upper = pdb_array{protein_i};
        path = sprintf('/Users/mingchenchen/Documents/Q40_Qbias/Q20_hexamer_Q400_430K_QW_separatestart_norg_500_nomem_pbc_rgoff_notseparate_allbeta');
        %path = sprintf('~/work/mdbox/qbias/%s',      pdbID_upper);    
        filename=sprintf('%s/%s_%d_pmf.dat',path, label, Tf1); FF=load(filename);
        if curve_shift_flag==0
            Fmin = min(FF(:,2)); id_shift = find( FF(:,2) == Fmin );    
        else
            [~,id_shift]=min(abs(FF(:,1)-q0_shift));
        end
        if qo_flag==0
            plot(FF(:,1), FF(:,2)-FF(id_shift(1),2), 'linewidth', 3, 'color',curve_color{protein_i}); grid on;
        end
        if qo_flag==1
                qx=FF(:,1);  Fy = FF(:,2); nbin=length(qx);
                dq=qx(2)-qx(1); qmin=qx(1)-dq/2; qmax= qx(nbin)+dq/2;
                Py=exp(-Fy/(0.001987*T_i)); P_norm = sum(Py); Py=Py/P_norm;
                pi_sample = zeros(Nsample,1); ni_sample = zeros(nbin, 1);
                %calculate pi_sample
                for i_bin= 1:nbin
                    qi_min = qmin + (i_bin-1)*dq; qi_max= qi_min + dq;
                    ids = find( q >= qi_min & q < qi_max ) ;    
                    ni_sample(i_bin) = length(ids);        
                    if ni_sample(i_bin) > 0
                        pi_sample(ids) = Py(i_bin)/ni_sample(i_bin);
                    end
                end    
                q_plot=linspace(0.1,0.9,nbin); dq=(0.9-0.1)/nbin; F_qo=ones(nbin,1);
                sum_p=0;
                for i_bin=1:nbin
                    qi_min = 0.1 + (i_bin-1)*dq; qi_max= qi_min + dq;
                    ids = find( qo >= qi_min & qo < qi_max ) ;
                    pi=sum(pi_sample(ids)); sum_p=sum_p + pi;
                    F_qo(i_bin)=-log(pi)*0.001987*T_i_plot;
                end
                %fprintf('check total probability: sum_p=%d\n', sum_p);
                plot(q_plot, F_qo, 'color', curve_color{protein_i}, 'linewidth', 3);
        end
        title(['T_f =', num2str(Tf),'K'],'fontsize', fsize);   
end


%plot F Vs. Qw
%{
binN_q = 50;
subplot(mr,mc,3), grid on, hold on, set(gca,'fontsize', fsize), xlabel('Q_{w}', 'fontsize', fsize); ylabel('F (kcal/mol)', 'fontsize', fsize)  , legend(['Tmin ', num2str(Tmin)], ['Tmax ', num2str(Tmax-dT)]);
q_lin=linspace(min(q), max(q), binN_q); [~, bin_q] = histc(q, q_lin);
for i_T=1:n_T
    T_i=Tmin+(i_T-1)*dT;    
    filename=sprintf('%s/%d_%d_pmf.dat',path, sim_label, T_i); FF=load(filename);
    qx=FF(:,1);  Fy = FF(:,2); nbin=length(qx);
    dq=qx(2)-qx(1); qmin=qx(1)-dq/2; qmax= qx(nbin)+dq/2;
    Py=exp(-Fy/(0.001987*T_i)); P_norm = sum(Py); Py=Py/P_norm;
    pi_sample = zeros(Nsample,1); ni_sample = zeros(nbin, 1);
    %calculate pi_sample
    for i_bin= 1:nbin
        qi_min = qmin + (i_bin-1)*dq; qi_max= qi_min + dq;
        ids = find( q >= qi_min & q < qi_max ) ;    
        ni_sample(i_bin) = length(ids);        
        if ni_sample(i_bin) > 0
            pi_sample(ids) = Py(i_bin)/ni_sample(i_bin);
        end
    end     
    p_q = zeros(binN_q, 1);
    for i_sample=1:Nsample
        pi=pi_sample(i_sample);     binq = bin_q(i_sample);
        p_q(binq) = p_q(binq) + pi ;
    end
    F_q = -0.001987*T_i*log(p_q);
    %if T_i == Tf
    %    Fmin = min(F_q); id_shift = find( F_q == Fmin );
    %end
    plot(q_lin, F_q-F_q(id_shift(1)), color(i_T), 'linewidth', 3); %hold on
    %filename=sprintf('%s/%d_%d_pmf.dat',path, sim_label, T_i); FF=load(filename);
    %plot(FF(:,1), FF(:,2)-FF(id_shift(1),2) , color(i_T+1), 'linewidth', 3); 
end
%}
