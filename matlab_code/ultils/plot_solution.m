function plot_solution(x, input, output, shape, key, alpha, beta, I0, q0, L, non_dim_w, non_dim_m)
    if shape == "rectangular"
        figure
        subplot(4,1,1)
        plot(x/L, (input(key,:,1)/I0).^(1/3))
        xlabel('$x/L$', Interpreter='latex')
        ylabel('$h/h_0$', Interpreter='latex')
        title(['$\beta$ = ', num2str(beta(key))], Interpreter='latex')
        ylim([-0.1, 1.1])
    elseif shape == "circular"
        figure
        subplot(4,1,1)
        plot(x/L, (input(key,:,1)/I0).^(1/4))
        xlabel('$x/L$', Interpreter='latex')
        ylabel('$h/h_0$', Interpreter='latex')
        title(['$\beta$ = ', num2str(beta(key))], Interpreter='latex')
        ylim([-0.1, 1.1])
    end
    subplot(4,1,2)
    yqmax = max(input(key,:,2)/q0);
    plot(x/L, input(key,:,2)/q0)
    xlabel('$x/L$', Interpreter='latex')
    ylabel('$q/q_0$', Interpreter='latex')
    title(['$\alpha$ = ', num2str(alpha(key))], Interpreter='latex')
    ylim([-0.1, yqmax+0.1])
    
    subplot(4,1,3)
    plot(x/L, output(key,:,1)*non_dim_w, 'b')
    xlabel('$x/L$', Interpreter='latex')
    ylabel('$w$', Interpreter='latex')
    grid on
    
    subplot(4,1,4)
    plot(x/L, output(key,:,2)*non_dim_m, 'b')
    xlabel('$x/L$', Interpreter='latex')
    ylabel('$M$', Interpreter='latex')
    grid on
    
    set(findobj(gcf,'type','axes'),"FontSize", 14)

end