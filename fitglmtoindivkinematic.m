split_time_before_trig=50;%ms
split_time_after_trig=50;%ms
split_safety_margin=150;%time_after_trig;
isi=10;
Wdname='R:\Whisker\Wdec\Wdec_by_field.mat';
split_binsize=1;%ms
smoothwidth=milliseconds(7.5);%ms
meth='gaussian';
bins=-split_time_before_trig:split_binsize:split_time_after_trig;
pbins=bins(1:end-1);
pbins=pbins+mean(diff(pbins)/2);
binsR=[-50 0 50];
splits={[0 100/3 200/3 100]};
splitname={'Tertile'};


% splits={[0 25 50 75 100];[0 100/3 200/3 100];[0 50 100]};
% splitname={'Quartile';'Tertile';'median'};


%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));clear temp*
[RasterHP,P_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,split_time_before_trig,split_time_after_trig,exclude,split_safety_margin,split_binsize,isi,[]);
[RasterPP,~]=get_raster_pop_no_self_safety(parameter,P_trig, DiscreteData,triglength_limit,ppms,split_time_before_trig,split_time_after_trig,exclude,split_safety_margin,split_binsize,isi,[]);

parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[RasterHT,T_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,split_time_before_trig,split_time_after_trig,exclude,split_safety_margin,split_binsize,isi,[]);


selc={vpmr;pomr};
selN={'VPM';'POm'};
selS={'Puff';'Touch'};
%%
split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Amp';'Angle';'Phase'};
mkabs=false;
normbefore=false;
fign_add='';
%%
split_var_nameC={'Acceleration';'Curvature';'Velocity';'Setp';'Angle';'Phase'};
mkabs=true;
normbefore=false;
fign_add='Abs_';
%%
split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Phase'};
mkabs=true;
normbefore=true;
fign_add='Abs_Norm';
%%
split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Amp';'Angle';'Phase'};
mkabs=false;
normbefore=true;
fign_add='Norm_';
%%
close all
sp_tr_Ta=cell(size(T_trig,1),numel(split_var_nameC));
sp_tr_Pa=cell(size(P_trig,1),numel(split_var_nameC));
for SN=1:numel(split_var_nameC)
    split_var_name=split_var_nameC{SN};
    % smoothker1=[1:1000 999:-1:0];
    % smoothkerLU=linspace(-smoothwidth/2,smoothwidth/2,numel(smoothker1));
    % smoothker2=interp1(smoothkerLU,smoothker1,bins,'linear',0);
    % smoothker2=smoothker2./sum(smoothker2);

    %


    switch split_var_name
        case 'Interval'
            sp_tr_T=cell(size(T_trig));
            for n=1:numel(RasterHT)
                if ~isempty(T_trig{n})
                    sp_tr_T{n}=cat(1,T_trig{n}(1),diff(T_trig{n}))/20000;
                end
            end
            sp_tr_P=cellfun(@(x) cat(1,x(1),diff(x))/20000,P_trig,'UniformOutput',0);

        case {'Angle','Setp','Phase','Curvature','Acceleration','Velocity','Amp'}
            switch split_var_name
                case {'Angle','Setp','Amp'}
                    eval(sprintf('a=%s;',split_var_name));
                    switch split_var_name
                        case {'Angle','Setp'}
                            a=cellfun(@(x) x-median(x,'omitnan'),a,'UniformOutput',0);
                    end
                case {'Acceleration','Velocity'}
                    a=load(Wdname, split_var_name);
                    a=a.(split_var_name);
                case 'Phase'
                    a=Phase;
                    for x=1:numel(a)
                        % a{x}(smoothdata(Amp{x},2,'movmedian',.5,'SamplePoints',milliseconds(1:numel(Amp{x})))<3)=nan;%cellfun(@(x) x-median(x),a,'UniformOutput',0);
                        a{x}(Amp{x}<3)=nan;%cellfun(@(x) x-median(x),a,'UniformOutput',0);
                    end
                case 'Curvature'
                    a=load(Wdname, split_var_name);
                    a=a.(split_var_name);
                    a=cellfun(@(x) x-median(x,'omitnan'),a,'UniformOutput',0);
            end
            pbinsx=pbins(pbins<35 & pbins>-35);
            [~,sp_tr_P]=get_trig_cont_pop(P_trig,pbinsx,a,20,DiscreteData);
            [~,sp_tr_T]=get_trig_cont_pop(T_trig,pbinsx,a,20,DiscreteData);
            xT=cellfun(@isempty,sp_tr_T);
            xP=cellfun(@isempty,sp_tr_P);
            sp_tr_P(xP)={nan(size(pbinsx))};
            sp_tr_T(xT)={nan(size(pbinsx))};
            if mkabs
                sp_tr_P=cellfun(@abs,sp_tr_P,'UniformOutput',0);
                sp_tr_T=cellfun(@abs,sp_tr_T,'UniformOutput',0);
            end
            if normbefore
                sp_tr_P=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan')-mean(x(:,pbinsx<0),2,'omitnan'),sp_tr_P,'UniformOutput',0);
                sp_tr_T=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan')-mean(x(:,pbinsx<0),2,'omitnan'),sp_tr_T,'UniformOutput',0);
            else
                sp_tr_P=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan'),sp_tr_P,'UniformOutput',0);
                sp_tr_T=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan'),sp_tr_T,'UniformOutput',0);
            end
            sp_tr_P(xP)={[]};
            sp_tr_T(xT)={[]};
    end

    sp_tr_Ta(:,SN)=cellfun(@(x) reshape(x,[],1),sp_tr_T,'UniformOutput',0);
    sp_tr_Pa(:,SN)=cellfun(@(x) reshape(x,[],1),sp_tr_P,'UniformOutput',0);
    sp_tr_Tall=cat(1,sp_tr_T{cat(1,selc{:})});
    sp_tr_Pall=cat(1,sp_tr_P{cat(1,selc{:})});

    %Ilim=prctile(cat(2,sp_tr_Pall,sp_tr_Tall),[0 100/3 200/3 100]);
    %Ilim=prctile(cat(2,sp_tr_Pall,sp_tr_Tall),[0 25 50 75 100]);
    %Ilim=prctile(cat(2,sp_tr_Pall,sp_tr_Tall),[0 100/2 100]);

    %
    %Ilim=[0 0.25 1.1 inf];
    %

    %Ilim=[0 inf];
    %%
end
%%
trial_rate_a=cell(numel(T_trig),2);
for n=1:numel(T_trig)
    if ~isempty(RasterHP{n}) && ~isempty(sp_tr_Pa{n,1}) && ~any(n==[80 81])
        tempR=zeros(numel(sp_tr_Pa{n,1}),2);
        for x=1:numel(sp_tr_Pa{n,1})
            tempR(x,:)= histcounts(RasterHP{n}(RasterHP{n}(:,2)==x,1),binsR)./(diff(binsR)/1000);
        end
        trial_rate_a{n,1}=tempR;
    end
    if ~isempty(RasterHT{n}) && ~isempty(sp_tr_Ta{n,1}) && ~any(n==[80 81])
        tempR=zeros(numel(sp_tr_Ta{n,1}),2);
        for x=1:numel(sp_tr_Ta{n,1})
            tempR(x,:)= histcounts(RasterHT{n}(RasterHT{n}(:,2)==x,1),binsR)./(diff(binsR)/1000);
        end
        trial_rate_a{n,2}=tempR;
    end
end
%%
includePred=[1 2 4];
test_prop=.75;
fit_reps=100;
incPredNames=split_var_nameC(includePred);
split_var_nameA=cat(1,split_var_nameC,{'isTouch'});
nx=cat(1,vpmr,pomr);
linCoef=cell(size(trial_rate_a,1),5);
predR=nan(size(trial_rate_a,1),4,5);
w=waitbar(0);
interval_I=find(includePred==3);
ftype='linear';
glms=linCoef;
for ny=1:numel(nx)
    n=nx(ny);

    %Puff
    pred=cat(2,sp_tr_Pa{n,:});
    pred(pred(:,interval_I)>1.25,interval_I)=1.25;
    resp=trial_rate_a{n,1}(:,2);%-mean(trial_rate_a{n,1}(:,1));
    resp=resp.*(diff(binsR(1:2))/1000);

    %resp=resp./max(resp);
    %pred=abs(pred);
    pred=pred-mean(pred,'omitnan');
    pred=pred./var(pred,0,1,'omitnan');
    [glms{n,1},linCoef{n,1},varExp,Corr,CorrP,varExpVar]=kinematic_glm(pred(:,includePred),resp,test_prop,split_var_nameA(includePred),fit_reps,[],ftype);
    predR(n,1,1)=Corr;
    predR(n,2,1)=CorrP;
  %  predR(n,3,1)=varExp;
   % predR(n,4,1)=var(varExpVar);

    %Touch
    pred=cat(2,sp_tr_Ta{n,:});
    pred(pred(:,interval_I)>1.25,interval_I)=1.25;
    resp=trial_rate_a{n,2}(:,2);%-mean(trial_rate_a{n,1}(:,1));
    resp=resp.*(diff(binsR(1:2))/1000);
    %resp=resp./max(resp);
    %pred=abs(pred);
     pred=pred-mean(pred,'omitnan');
    pred=pred./var(pred,0,1,'omitnan');
    [~,linCoef{n,2},varExp,Corr,CorrP,varExpVar]=kinematic_glm(pred(:,includePred),resp,test_prop,split_var_nameC(includePred),fit_reps,[],ftype);
    predR(n,1,2)=Corr;
    predR(n,2,2)=CorrP;
 %   predR(n,3,2)=varExp;
  %  predR(n,4,2)=var(varExpVar);

    %Puff/Touch joint
    pred=cat(1,cat(2,sp_tr_Pa{n,:}),cat(2,sp_tr_Ta{n,:}));
      pred(pred(:,interval_I)>1.25,interval_I)=1.25;
    resp=cat(1,trial_rate_a{n,1}(:,2),trial_rate_a{n,2}(:,2));
    resp=resp.*(diff(binsR(1:2))/1000);
     pred=pred-mean(pred,'omitnan');
    pred=pred./var(pred,0,1,'omitnan');
    [~,linCoef{n,3},varExp,Corr,CorrP,varExpVar]=kinematic_glm(pred(:,includePred),resp,test_prop,split_var_nameA(includePred),fit_reps,[],ftype);
    predR(n,1,3)=Corr;
    predR(n,2,3)=CorrP;
 %   predR(n,3,3)=varExp;
 %   predR(n,4,3)=var(varExpVar);

    %Puff/Touch joint abs
    pred=cat(1,cat(2,sp_tr_Pa{n,:}),cat(2,sp_tr_Ta{n,:}));
   pred(pred(:,interval_I)>1.25,interval_I)=1.25;
    resp=cat(1,trial_rate_a{n,1}(:,2),trial_rate_a{n,2}(:,2));
    resp=resp.*(diff(binsR(1:2))/1000);
    
     pred=pred-mean(pred,'omitnan');
    pred=pred./var(pred,0,1,'omitnan');
    pred=abs(pred);
    [~,linCoef{n,4},varExp,Corr,CorrP,varExpVar]=kinematic_glm(pred(:,includePred),resp,test_prop,split_var_nameA(includePred),fit_reps,[],ftype);
    predR(n,1,4)=Corr;
    predR(n,2,4)=CorrP;
  %  predR(n,3,4)=varExp;
  %  predR(n,4,4)=var(varExpVar);

    %Puff/Touch joint with identity

   
    pred=cat(2,pred,zeros(size(pred,1),1));
    pred(size(trial_rate_a{n,1}(:,2),1)+1:end,end)=1;

    includePredS=[includePred size(pred,2)];
    [~,linCoef{n,5},varExp,Corr,CorrP,varExpVar]=kinematic_glm(pred(:,includePredS),resp,test_prop,split_var_nameA(includePredS),fit_reps,numel(includePredS),ftype);
    predR(n,1,5)=Corr;
    predR(n,2,5)=CorrP;
  %  predR(n,3,5)=varExp;
  %  predR(n,4,5)=var(varExpVar);

    waitbar(ny/numel(nx),w)
end
close(w)
%%
% figure, scatter(resp,predict(a,pred))
test=cat(3,linCoef{pomr,5});
test=mean(squeeze(test(2:end,2,:)<.05),2)
incPredNames
%%
figure,
tiledlayout(2,3)
for n=1:5
    nexttile;
    histogram(predR(vpmr,1,n).^2,0:.025:.5);hold on;
    histogram(predR(pomr,1,n).^2,0:.025:.5)
end
%%
figure,
tiledlayout(1,2)
nexttile;
plot(1:2,permute(predR(vpmr,1,[1 4]),[1 3 2]).^2)
nexttile;
plot(1:2,permute(predR(pomr,1,[1 4]),[1 3 2]).^2)
%%
n=vpmr(1);
pred=cat(1,cat(2,sp_tr_Pa{n,:}),cat(2,sp_tr_Ta{n,:}));
resp=cat(1,trial_rate_a{n,1}(:,2),trial_rate_a{n,2}(:,2));
resp=resp.*(diff(binsR(1:2))/1000);
pred=abs(pred);
pred=pred-mean(pred,'omitnan');
pred=pred./std(pred,0,1,'omitnan');
pred=cat(2,pred,zeros(size(pred,1),1));
pred(size(trial_rate_a{n,1}(:,2),1)+1:end,end)=1;
[a]=fitglm(pred,resp,'interactions','VarNames',cat(1,split_var_nameC,{'isTouch'},{'rate'}),'CategoricalVars',9,'Distribution','poisson','PredictorVars',[1 2 4 3 9]);
[c,b]=corr(resp,predict(a,pred));
linCoef{n,3}=cat(2,a.Coefficients{:,'Estimate'},a.Coefficients{:,'pValue'});
predR(n,1,3)=c;
predR(n,2,3)=b;



%%
pred=cat(1,cat(2,sp_tr_Pa{n,:}),cat(2,sp_tr_Ta{n,:}));
resp=cat(1,trial_rate_a{n,1}(:,2)-mean(trial_rate_a{n,1}(:,1)),trial_rate_a{n,2}(:,2)-mean(trial_rate_a{n,2}(:,1)));

% resp=cat(1,trial_rate_a{n,1}(:,2),trial_rate_a{n,2}(:,2));
% resp=resp.*(diff(binsR(1:2))/1000);

%pred=abs(pred);
%pred=pred-mean(pred,'omitnan');
pred=pred./std(pred,0,1,'omitnan');
[a]=fitglm(pred,resp,'purequadratic','VarNames',cat(1,split_var_nameC,{'rate'}),'Distribution','normal','PredictorVars',[1:4 7])
%%
poi=nan(numel(nx),1);
for ny=1:numel(nx)
    n=nx(ny);

    %Puff
    pred=cat(2,sp_tr_Pa{n,:});
    pred(pred(:,interval_I)>1.25,interval_I)=1.25;
    resp=trial_rate_a{n,1}(:,2);%-mean(trial_rate_a{n,1}(:,1));
    resp=resp.*(diff(binsR(1:2))/1000);
    poi(ny)=var(resp)/mean(resp);
end