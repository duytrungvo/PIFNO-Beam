function plot_domain(x, area_range, func_name)

figure

area(x, area_range, FaceAlpha=0.5, LineStyle="none")
newcolors = [1 1 1; 0.7 0.7 0.7];
% newcolors = [1 1 1; 0.3010 0.7450 0.9330];
% newcolors = [1 1 1; 0 0 1];
% newcolors = [1 1 1; 0 0.4470 0.7410]
colororder(newcolors)
xlabel('$x$', 'Interpreter','latex')
ylabel(func_name, 'Interpreter','latex')
% ylim([0,2])
set(findobj(gcf,'type','axes'),"FontSize", 18)

end