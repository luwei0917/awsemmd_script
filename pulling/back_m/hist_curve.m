function [p_x, x_grid] = hist_curve(x, minx, maxx, binN, plot_flag)
x_grid = linspace(minx, maxx, binN);
x_grid = x_grid' ;
delta_x = (maxx-minx)/binN;
x_bin_count = histc(x, x_grid);
p_x = x_bin_count / length(x) /delta_x ;
if plot_flag==1
    figure,
    plot(x_grid, p_x);
end
   