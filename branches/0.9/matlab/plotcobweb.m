function plotcobweb(x)

x = x(:)';
x = repmat(x, [2 1]);
x = [x(:) x(:)];
figure; hold on;
plot(x(2:end-1,1),x(3:end,2),'k')
plot([min(x(:,1)) max(x(:,1))], [min(x(:,1)) max(x(:,1))], 'Color', [.75 .75 .75], 'LineWidth', 2);
