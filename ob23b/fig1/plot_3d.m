function [fig] = plot_3d(y, fitness, range, num_eval)
    
    y = reshape(y,1,length(y));
    %num_eval = 100;
    
    x_1 = linspace(y(1)-range, y(1)+range, num_eval);
    x_2 = linspace(y(2)-range, y(2)+range, num_eval);

    [X_1, X_2] = meshgrid(x_1, x_2);
    X_0_mesh = X_1;
    X_1_mesh = X_2;
    X_1 = X_1(:)';
    X_2 = X_2(:)';
    X_i = y(3:end)';
    X_i = repmat(X_i,1,num_eval^2);

    X = [X_1;X_2;X_i];

    F = fitness(X);
    F_mesh = reshape(F, [num_eval, num_eval]);
    
    fig = figure; 
    
    hold on; 
    xlabel('$y_1$', 'interpreter','latex'); 
    ylabel('$y_2$', 'interpreter','latex'); 
    zlabel('$f(\textbf{\emph y})$', 'interpreter','latex');
    
    surf(X_0_mesh, X_1_mesh, F_mesh); 

        
end
 

