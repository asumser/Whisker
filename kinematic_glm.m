function [glms,linCoefs,varExp,Corr,CorrP,varExpVar]=kinematic_glm(pred,resp,test_prop,split_var_nameC,fit_reps,catVars,ftype)

glms=cell(fit_reps+1,1);
    CorrP=nan(fit_reps,1);
    Corr=nan(fit_reps,1);
    varExp=nan(fit_reps,1);
    varExpVar=nan(fit_reps,1);
    linCoefs=nan(size(pred,2)+1,2);
    distr='poisson';
    re=any(isnan(pred),2);
    resp(re)=[];
    pred(re,:)=[];
 if sum(resp>0)/numel(resp)>.1
    parfor r=1:fit_reps
        
    trainsel=randperm(numel(resp),round(test_prop*numel(resp)));
    testsel=1:numel(resp);
    testsel(trainsel)=[];
    
    [glms{r}]=fitglm(pred(trainsel,:),resp(trainsel),ftype,'VarNames',cat(1,split_var_nameC,{'rate'}),'Distribution',distr,'CategoricalVars',catVars);
%     [test]=stepwiseglm(pred(trainsel,:),resp(trainsel),'linear','VarNames',cat(1,split_var_nameC,{'rate'}),'Distribution','normal','CategoricalVars',catVars);
%     [test,FitInfo] = lassoglm(pred,resp,'poisson','CV',10);
%figure, scatter(resp(testsel),predict(glms{r},pred(testsel,:)))
    [Corr(r),CorrP(r)]=corr(resp(testsel),predict(glms{r},pred(testsel,:)));
   % title(sprintf('rho=%.2f',Corr(r)))
    %
   % var(resp(testsel)-predict(glms{r},pred(testsel,:)))./var(resp(testsel))
    varExp(r)=1-var(resp(testsel)-predict(glms{r},pred(testsel,:)))./var(resp(testsel));
    end
  
     [glms{fit_reps+1}]=fitglm(pred,resp,'linear','VarNames',cat(1,split_var_nameC,{'rate'}),'Distribution',distr,'CategoricalVars',catVars);
    linCoefs=cat(2,glms{end}.Coefficients{:,'Estimate'},glms{end}.Coefficients{:,'pValue'});

    varExpVar=varExp;
    varExp=mean(varExp,'omitnan');
%     figure, histogram(Corr)
    Corr=mean(Corr,'omitnan');
    CorrP=mean(CorrP,'omitnan');
 else
     CorrP=nan(1,1);
    Corr=nan(1,1);
    varExp=nan(1,1);
 end