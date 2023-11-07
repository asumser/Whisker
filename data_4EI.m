Puff_PSTH=P_PSTH_instrate([vpmr;pomr],:);
Touch_PSTH=T_PSTH_instrate([vpmr;pomr],:);
Bins_PSTH=-psth_time_before_trig:psth_binsize:psth_time_after_trig;
Bins_PSTH=Bins_PSTH(1:end-1)+psth_binsize/2;
Puff_Raster=RasterHP([vpmr;pomr]);
Touch_Raster=RasterHT([vpmr;pomr]);
Raster_column_key={'relative spike time','trial no'};
Cell_numbers=[vpmr;pomr];
is_vpm=false(size(Cell_numbers));is_vpm(1:numel(vpmr))=true;
is_pom=false(size(Cell_numbers));is_pom(numel(vpmr)+1:end)=true;
save('vpm_pom_cluster_data_for_Emilio_Oct_23.mat','Puff_PSTH','Touch_PSTH',...
    'Bins_PSTH','Puff_Raster','Touch_Raster','Raster_column_key','Cell_numbers','is_vpm','is_pom')