clear
load('data.mat');

P_target = 0.99;

num_var = length(data);
try data.N;
    str_var = '$N$';
    VAR = [data(:).N];
end
try data.A;
    str_var = '$A$';
    VAR = [data(:).A];
end
try data.ALPHA;
    str_var = '$\alpha$';
    VAR = [data(:).ALPHA];
end   

figure; hold on; legend;
for i=1:num_var
    plot(data(i).MU, data(i).PS, 'DisplayName', [str_var, '=', num2str(VAR(i))]);
end
legend('location','NorthWest')
set(gca, 'XScale', 'log');
xlabel('$\mu$');
ylabel('$P_S$'); yticks([0:0.2:1]);
myfigstyle(gcf, 8, 4.5, 9, 9);

MU_givenP = nan*zeros(num_var, 1);
for i=1:num_var
    ids = find(data(i).PS >= P_target);
    if ~isempty(ids)
        MU_givenP(i) = data(i).MU(ids(1));
    end
end

get_scaling_sim(VAR, MU_givenP, P_target, 'firstCall', str_var, 'log');
% get lines
get_scaling_lines(VAR, 2);
% plot again for better visibility
get_scaling_sim(VAR, MU_givenP, P_target, 'secondCall', str_var, 'log');

xticks([30,50,100,200,300])

myfigstyle_scaling(gcf, 8, 4.5, 9, 9);
legend('off');
xlim([VAR(1),VAR(end)])
% xticks([VAR(1):2:VAR(end)])
% yticks([0:200:1000])
% thicker sim line
Line = findall(gcf, 'type', 'line');
Line(1).LineWidth=1.5;