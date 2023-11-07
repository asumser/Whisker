function print_sum_stat(values)

fprintf('mean: %.3f std: %.3f sem: %.3f median: %.3f iqr: %.3f n = %u\n',...
    mean(values,1,'omitnan'),std(values,0,1,'omitnan'),std(values,0,1,'omitnan')/sqrt(sum(~isnan(values))),median(values,1,'omitnan'),iqr(values,1),sum(~isnan(values)))