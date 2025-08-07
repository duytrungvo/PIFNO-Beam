function plot_domain_v2(x, y, alpha, area_range, func_name)

figure
% nexttile
% plot(x, y)
% hold on
% nexttile
area(x, area_range, FaceAlpha=0.5, LineStyle="none")
newcolors = [1 1 1; 0.7 0.7 0.7];
% newcolors = [1 1 1; 0.3010 0.7450 0.9330];
% newcolors = [1 1 1; 0 0 1];
% newcolors = [1 1 1; 0 0.4470 0.7410];
% newcolors = [0.7 0.7 0.7; 1 1 1];
colororder(newcolors)
hold on
plot(x, y, 'k-.')
% for parabolic symmetric varying in depth 
% text(x(end-515)*ones(length(alpha),1), y(:,end-515), num2str(round(alpha,2)))
% for sinusoidal in load
% text(x(end-528)*ones(length(alpha),1), y(:,end-528), num2str(round(alpha,2)))
text(x(end-200)*ones(length(alpha),1), y(:,end-200), num2str(round(alpha,2)))
hold off
xlabel('$x$', 'Interpreter','latex')
ylabel(func_name, 'Interpreter','latex')
% ylim([0,2])
set(findobj(gcf,'type','axes'),"FontSize", 14)
% hold off
end