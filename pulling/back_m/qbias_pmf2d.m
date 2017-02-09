%pdb_array= {'Abeta42_2' 'Abeta420' 'Abeta40_2'
pdb_array= {'q30'};
%sim_labels = [4 4 4];
sim_labels = [12];

% what data to load
q_name={  'qw_ga_total' ,'qw_total'};
T=300;
%xy_limit=1; % whose range to use for plotting, 1 or 2?

n_contour=20; % # of contour lines
fsize=25; tsize=20; mr=1; mc=length(sim_labels);
cutoff=20;scrnsize = get(0,'ScreenSize');
figure('position', [1 scrnsize(4) 0.25*mc*scrnsize(3) 0.35*scrnsize(4)]);

for i_label=1:length(sim_labels)
    protein_i = i_label; pdbID_upper = pdb_array{protein_i};
    path = sprintf('.');
    sim_label = sim_labels(i_label);
    qa_name=q_name{1}; qb_name=q_name{2};

    filename = sprintf('%s/%s',path, qa_name); qa = load(filename);
    filename = sprintf('%s/%s',path, qb_name); qb = load(filename);
    if strcmp(q_name{1},'dih')
        qa=(mean(qa'))';
    end
    if strcmp(q_name{2},'dih')
        qb=(mean(qb'))';
    end
    if i_label==1
        qa_min=min(qa); qa_max=max(qa);
        qb_min=min(qb); qb_max=max(qb);
    end
    Nsample=size(qa,1);
    assert(Nsample==size(qb,1));
    %load pmf file and calculate pi_sample
    T_i=T; T=T_i;
    filename = sprintf('%s/p_total',path); q = load(filename);
    filename=sprintf('%s/2lhd_%d_pmf.dat',path, T_i);
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
    fprintf('probability = %.3f\n', sum(pi_sample));

    binN=20;ObinN=20;
    qa_lin=linspace(min(qa), max(qa),ObinN); qb_lin=linspace(min(qb), max(qb),binN); H=zeros(ObinN,binN);
    [~, bin_index_x] = histc(qa, qa_lin); [~, bin_index_y] = histc(qb, qb_lin);
    for i_sample = 1:Nsample
        x=bin_index_x(i_sample); y=bin_index_y(i_sample);
        H(x,y) = H(x,y) + pi_sample(i_sample);
    end
    H=H'; fprintf('sum(sum(H))=%.3f\n', sum(sum(H)));
    F=-0.001987*T*log(H); ids = (F>= cutoff); F(ids) = cutoff;
    subplot(mr,mc,i_label)
    F
    [~,h] = contourf(qa_lin, qb_lin,F,n_contour); shading flat,
    colormap(jet), colorbar,
    title([num2str(pdbID_upper), ' T= ', num2str(T)],'fontsize', fsize);
    %%%fill the top area white
    ccc = get(h,'children'); max_cdata = -inf; cdata_list=zeros(size(ccc,1), 1);

    for k=1:size(ccc,1)
        cd1 = get(ccc(k), 'cdata');
        if cd1 > max_cdata
            max_cdata = cd1 ;
        end
        cdata_list(k) = get(ccc(k),'cdata');
    end
    id = find(cdata_list == max_cdata);
    disp(ccc(id));
    for k=1:size(id,1)
        set(ccc(id(k)), 'facecolor', 'white');
    end

    xlabel(q_name{1}, 'fontsize', fsize),
    ylabel(q_name{2}, 'fontsize', fsize);
    set(gca, 'FontSize', fsize);
    xlim([qa_min, qa_max])
    ylim([qb_min, qb_max])
end
