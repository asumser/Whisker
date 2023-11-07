function print_avg_comp(values,vpm,pom)

vpm_values=values(vpm,:);
vpm_values=vpm_values(all(~isnan(vpm_values),2),:);
p_vpm=signrank(vpm_values(:,1),vpm_values(:,2));
pom_values=values(pom,:);
pom_values=pom_values(all(~isnan(pom_values),2),:);
p_pom=signrank(pom_values(:,1),pom_values(:,2));

p_comp=cat(2,ranksum(vpm_values(:,1),pom_values(:,1)),ranksum(vpm_values(:,2),pom_values(:,2)));
fprintf('vpm: %.3f Hz +/- %.3f => %.3f Hz +/- %.3f p=%.6f n= %u\n pom: %.3f Hz +/- %.3f => %.3f Hz +/- %.3f p=%.6f n= %u\n base_comp_p = %.6f resp_comp_p = %.6f\n',...
    mean(vpm_values(:,1),1,'omitnan'),std(vpm_values(:,1),0,1,'omitnan')/sqrt(size(vpm_values,1)),mean(vpm_values(:,2),1,'omitnan'),std(vpm_values(:,2),0,1,'omitnan')/sqrt(size(vpm_values,1)),p_vpm,size(vpm_values,1),...
    mean(pom_values(:,1),1,'omitnan'),std(pom_values(:,1),0,1,'omitnan')/sqrt(size(pom_values,1)),mean(pom_values(:,2),1,'omitnan'),std(pom_values(:,2),0,1,'omitnan')/sqrt(size(pom_values,1)),p_pom,size(pom_values,1),...
    p_comp(1),p_comp(2))