pdb_array= {'q30' 'IkBa' 'Abeta42_3' 'Abeta40_3' 'Abeta42_40_3'};
repeat = 18; k_bias=200;
labels = [10 9 7 14 13 10 4 4 ];

for protein_i=1:1
    fsize=20; tsize=14; mr=3; mc=1;
    color=colormap(jet(repeat));
    label = labels(protein_i); pdbID_upper = pdb_array{protein_i};
    path = sprintf('.');

    %load q
    for i_label=1:length(label)
        scrnsize = get(0,'ScreenSize'); figure('position', [1 scrnsize(4) 0.5*scrnsize(3) scrnsize(4)]);
        filename = sprintf('%s/p_total',path); q = load(filename);
        Nline= length(q)/repeat;
        %cat 1e_total 2e_total 3e_total 4e_total 5e_total 6e_total 7e_total 8e_total 9e_total 10e_total >> e_total

        %plot q
        subplot(mr,mc,1), grid on, hold on, set(gca, 'fontsize', tsize);
        plot(q,'k'); xlabel('frames', 'fontsize', fsize);ylabel('Q', 'fontsize', fsize);
        %title([pdbID_upper, ' sim ', num2str(label(i_label)), '  N-qbias=',num2str(repeat), ' k-bias=',num2str(k_bias)], 'fontsize', fsize);
        title([pdbID_upper, '  N-qbias=',num2str(repeat), ' k-bias=',num2str(k_bias)], 'fontsize', fsize);
        %load energy
        for i_label=1:length(label)
            filename=sprintf('%s/e_total', path);
            E = load(filename);
        end

        %plot Q distribution
        subplot(mr,mc,2), grid on, hold on, set(gca, 'fontsize', tsize); xlabel('Q', 'fontsize', fsize);ylabel('P(Q)', 'fontsize', fsize);
        binN_Qi=50;
        for i_repeat=1:repeat
            i_start = 1 + Nline*(i_repeat -1); i_end = i_start + Nline - 1;
            qwindow=q(i_start:i_end,1);
            [p_q,qgrid]=hist_curve(qwindow, min(qwindow),max(qwindow), binN_Qi,0);
            fprintf('%d %f %f %f\n', i_repeat, min(qwindow), max(qwindow), mean(qwindow));
            if i_repeat<repeat % plot the last one outside of the loop
                plot(qgrid, p_q, 'color', color(i_repeat,:), 'linewidth', 1.5);
            end
        end
        [p_q_total, qgrid_total]=hist_curve(q, min(q),max(q), 50,0);
        %plot(qgrid_total, p_q_total, '--k', 'linewidth', 3);
        [AX,H1,H2]=plotyy(qgrid, p_q, qgrid_total, p_q_total);
        set(H1,'color', color(repeat,:), 'linewidth',3);
        set(H2,'linestyle','--', 'linewidth',4);
        set(AX(2),'fontsize', tsize);
        set(get(AX(2),'ylabel'), 'string','P(Q)','fontsize', fsize);
        %plot E distribution
        subplot(mr,mc,3), grid, hold on, set(gca, 'fontsize', tsize); xlabel('E (kcal/mol)', 'fontsize', fsize);ylabel('P(E)', 'fontsize', fsize);
        binN_Ei = 100;
        for i_repeat=1:repeat
            i_start = 1+(i_repeat-1)*Nline;
            i_end   = i_repeat*Nline;

            Ewindow=E(i_start:i_end, 1);
            [E_q,Egrid]=hist_curve(Ewindow, min(Ewindow),max(Ewindow), binN_Ei,0);
            plot(Egrid, E_q, 'color', color(i_repeat,:), 'linewidth', 1.5);
        end
%         [p_E_total, Egrid_total]=hist_curve(E, min(E),max(E), 100,0);
%         %plot(qgrid_total, p_q_total, '--k', 'linewidth', 3);
%         [AX,H1,H2]=plotyy(Egrid, E_q, Egrid_total, p_E_total);
%         set(H1,'color', color(repeat,:), 'linewidth',3);
%         set(H2,'linestyle','--', 'linewidth',3);
%         set(AX(2),'fontsize', tsize);
%         set(get(AX(2),'ylabel'), 'string','P(Q)','fontsize', fsize);
    end
end
