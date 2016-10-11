%pdb_array= { 'Abeta42_1' 'Abeta40_1'};
pdb_array= { 'q20' 'Abeta40_2' 'Abeta40_3' 'Abeta40_6_nonfibril'};
sim_labels = [14 7 6 4];

qa_name ='tc_total';%new rxn co
T=400;
binN=20;
color={'r','m','g','b'};
dih_range=16:22;

curve_shift_flag=1; q0_shift=-32;
%; qmono_flag =0;

%if qmono_flag == 1
%    path = sprintf('~/work/mdbox/qbias_mono/%s', pdbID_upper);
%else
%    path = sprintf('~/work/mdbox/qbias/%s',      pdbID_upper);
%end
cutoff=12;
scrnsize = get(0,'ScreenSize'); 
figure('position', [1 scrnsize(4) 0.3*scrnsize(3) 0.4*scrnsize(4)]), hold on
fsize=25; tsize=16; %mr=3; mc=2;
%for i_list = 1:length(T_list)
    %T = T_list(i_list);
    %subplot(mr, mc, i_list),
    %grid on, hold on, set(gca, 'fontsize', tsize); xlabel('Q', 'fontsize', fsize); ylabel('F (kcal/mol)', 'fontsize', fsize); 
    for i_label=1
        sim_label = sim_labels(i_label);
        pdbID_upper = pdb_array{i_label};
        path = sprintf('/Users/weilu/Documents/Research/server/top7_t400_q100');
        %title([pdbID_upper, ' ', num2str(T), 'K'], 'fontsize', fsize);
        
        filename = sprintf('%s/%s',path, qa_name); qa = load(filename);        
        if strcmp(qa_name,'dih')
            qa=(mean(qa(:,dih_range)'))';
        end
        
        %load q
        filename = sprintf('%s/p_total',path); q = load(filename);
        %if qo_flag == 1
        %    filename = sprintf('%s/qo_%d',path, sim_label); qo = load(filename); 
        %end
        Nsample = length(q);    

        %load pmf file and calculate pi_sample
        filename=sprintf('%s/q20_%d_pmf.dat',path, T);        

        FF=load(filename); qx=FF(:,1);  Fy = FF(:,2); nbin=length(qx);
        dq=qx(2)-qx(1); qmin=qx(1)-dq/2; qmax= qx(nbin)+dq/2;
        Py=exp(-Fy/(0.001987*T)); P_norm = sum(Py); Py=Py/P_norm;
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
        fprintf('probability = %.3f\n', sum(pi_sample));                        

        qa_lin=linspace(min(qa), max(qa),binN);
        count_qa=zeros(binN,1);
        [~,bin_index_x]=histc(qa, qa_lin);
        for i_sample=1:Nsample
            x=bin_index_x(i_sample);
            count_qa(x) = count_qa(x) + pi_sample(i_sample);
        end
        fprintf('Total probability for new coordinate is %.3f\n',sum(count_qa));
        F_qa=-0.001987*T*log(count_qa); ids = (F_qa>= cutoff); F_qa(ids) = cutoff;
        if curve_shift_flag==0
            Fmin = min(F_qa); id_shift = find( F_qa == Fmin );
        else
            [~,id_shift]=min(abs(qa_lin-q0_shift));
        end

        
        %plot F Vs. Qo, shift free energy curve    
        %Fmin = min(FF(:,2)); id_shift = find( FF(:,2) == Fmin );    
        %xlim([0 1]);
        %id_shift(1)=35;
        %plot(FF(:,1), FF(:,2)-FF(id_shift(1),2), 'linewidth', 3); 
        %if i_list == 1
        %    legend('T400','T420','T440','T450','fontsize',fsize);
        %end
        plot(qa_lin, F_qa-F_qa(id_shift(1)), color{i_label}, 'linewidth', 3);   
        fsize=30; tsize=16;
        xlabel('Number of Total Contacts', 'fontsize', fsize); ylabel('Free energy (kcal/mol)', 'fontsize', fsize); title('Free energy with contact number', 'fontsize', fsize);
   
    end