function [glms,linCoefs,varExp,Corr,CorrP,varExpVar]=kinematic_glm(pred,resp,test_prop,split_var_nameC,fit_reps,catVars,ftype)

glms=cell(fit_reps+1,1);
    CorrP=nan(fit_reps,1);
    Corr=nan(fit_reps,1);
    varExp=nan(fit_reps,1);
 
    parfor r=1:fit_reps
    trainsel=randperm(numel(resp),round(test_prop*numel(resp)));
    testsel=1:numel(resp);
    testsel(trainsel)=[];
    [glms{r}]=fitglm(pred(trainsel,:),resp(trainsel),ftype,'VarNames',cat(1,split_var_nameC,{'rate'}),'Distribution','poisson','CategoricalVars',catVars);
%     [test]=stepwiseglm(pred(trainsel,:),resp(trainsel),'linear','VarNames',cat(1,split_var_nameC,{'rate'}),'Distribution','normal','CategoricalVars',catVars);
%     [test,FitInfo] = lassoglm(pred,resp,'poisson','CV',10);
    [Corr(r),CorrP(r)]=corr(resp(testsel),predict(glms{r},pred(testsel,:)));
    varExp(r)=1-var(resp(testsel)-predict(glms{r},pred(testsel,:)))./var(resp(testsel));
    end
  
     [glms{fit_reps+1}]=fitglm(pred,resp,'linear','VarNames',cat(1,split_var_nameC,{'rate'}),'Distribution','poisson','CategoricalVars',catVars);
    linCoefs=cat(2,glms{end}.Coefficients{:,'Estimate'},glms{end}.Coefficients{:,'pValue'});

    varExpVar=varExp;
    varExp=mean(varExp);
    Corr=mean(Corr);
    CorrP=mean(CorrP);