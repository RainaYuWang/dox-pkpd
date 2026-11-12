function results = run_pkpd_ci(mode)
% RUN_PKPD_CI Run the unchanged quick/full fit and check its saved outputs.
if nargin<1,mode='quick';end
mode=validatestring(mode,{'quick','full','sensitivity','dose-beta'});
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root);
out=fullfile(root,'joint_results');
% Never mix a previous run's files with the current run or delete user results.
if exist(out,'dir')
 error('PKPD:ExistingResults','Move existing joint_results aside before a new CI run.');
end
mkdir(out);
timer=tic;
doSensitivity=strcmp(mode,'sensitivity');
doDoseBeta=strcmp(mode,'dose-beta');
nStarts=2;if ~strcmp(mode,'quick'),nStarts=8;end
status=struct('status','running','mode',mode,'nStarts',nStarts,...
 'seed',42,'matlab_version',version,'platform',computer,...
 'commit',getenv('GITHUB_SHA'),'run_id',getenv('GITHUB_RUN_ID'),...
 'started_utc',char(datetime('now','TimeZone','UTC')),'elapsed_seconds',0);
status.command=sprintf('run_joint_pkpd(%d,%d,%d)',nStarts,doSensitivity,doDoseBeta);
status.rng_algorithm='twister';
writeStatus(out,status);
try
 results=run_joint_pkpd(nStarts,doSensitivity,doDoseBeta);
 % These are runtime/output checks, not tests of predictive scientific validity.
 required={'01_cell_fits.png','02_residuals.png','03_PK_reference.png',...
  '04_equivalent_parameters.png','representative_parameters.csv',...
  'growth_parameters.csv','equivalent_parameter_sets.csv',...
  'reversible_combinations.csv','old_pk_multistart.csv',...
  'joint_irreversible_multistart.csv','joint_reversible_multistart.csv',...
  'predictions.csv','residuals.csv','fit_by_dose.csv','search_bounds.csv',...
  'joint_fit.mat','run_inputs.mat','run.log','model_d_fit.mat',...
  'model_d_fixed_parameters.csv','model_d_multistart.csv',...
  'model_d_kbf_profile.csv','model_d_kbf_profile.mat','05_model_d_kbf_profile.png'};
 for j=1:numel(required)
  info=dir(fullfile(out,required{j}));
  assert(~isempty(info) && info.bytes>0,'PKPD:MissingOutput',...
   'Missing or empty output: %s',required{j});
 end
 expectedStarts=2*min(4,nStarts)+3*nStarts+1;
 checkpoints=dir(fullfile(out,'checkpoints','*_start_*.mat'));
 assert(numel(checkpoints)==expectedStarts,'PKPD:Checkpoints',...
  'Expected %d completed starts; found %d.',expectedStarts,numel(checkpoints));
 summary=readtable(fullfile(out,'representative_parameters.csv'));
 assert(height(summary)==4 && all(isfinite(summary.scaled_SSE)) && ...
  all(summary.scaled_SSE>=0) && all(summary.scaled_SSE<1e9),...
  'PKPD:InvalidFit','Invalid fit summary.');
 predictions=readtable(fullfile(out,'predictions.csv'));
 values=predictions{:,{'CF','CB','N'}};
 assert(height(predictions)>0 && all(isfinite(values(:))) && ...
  all(predictions.N>0),'PKPD:InvalidPredictions','Invalid predictions.');
 saved=load(fullfile(out,'joint_fit.mat'),'results');
 assert(isequal(saved.results.summary,results.summary),...
  'PKPD:SavedResults','Saved MAT summary differs from returned result.');
 profile=readtable(fullfile(out,'model_d_kbf_profile.csv'));
 assert(height(profile)>=100 && all(isfinite(profile.scaled_SSE)) && ...
  all(profile.scaled_SSE<1e9),'PKPD:Profile','Incomplete/invalid D profile.');
 assert(abs(min(profile.scaled_SSE)-summary.scaled_SSE(4))<1e-4,...
  'PKPD:ProfileOptimum','D joint fit and beta profile optimum disagree.');
 assert(isequal(results.model_d_growth,[0.031094959444983555,35143.615284616222]),...
  'PKPD:GrowthReference','D must use the archived full growth fit.');
 if doSensitivity
  sens=readtable(fullfile(out,'model_d_sensitivity_summary.csv'));
  starts=readtable(fullfile(out,'model_d_sensitivity_starts.csv'));
  curves=readtable(fullfile(out,'model_d_sensitivity_predictions.csv'));
  assert(height(sens)==9 && height(starts)==9*nStarts,...
   'PKPD:SensitivityCount','Missing sensitivity fits.');
  assert(all(isfinite(sens.scaled_SSE)) && all(sens.scaled_SSE<1e9),...
   'PKPD:SensitivityFit','Invalid sensitivity fit.');
  vals=curves{:,{'CF','CB','N','beta_CB'}};
  assert(all(isfinite(vals(:))) && all(curves.N>0),...
   'PKPD:SensitivityCurves','Invalid sensitivity trajectories.');
  reference=[3.5e-6,1.25,0.21];names={'kEF','kFE','kFB'};
  for k=1:9
   expected=reference;which=find(strcmp(names,sens.varied_parameter{k}));
   expected(which)=expected(which)*sens.multiplier(k);
   assert(isequal(sens{k,{'fixed_kEF','fixed_kFE','fixed_kFB'}},expected),...
    'PKPD:SensitivityFixed','More than one reference parameter changed.');
  end
  refs=sens(sens.multiplier==1,:);
  assert(max(refs.scaled_SSE)-min(refs.scaled_SSE)<1e-10,...
   'PKPD:SensitivityRepeat','Identical reference fits disagree.');
  status.sensitivity_cases=height(sens);status.sensitivity_starts=height(starts);
 end
 if doDoseBeta
  diagnostic=load(fullfile(out,'model_d_dose_beta.mat'));
  assert(diagnostic.kBF==0.0049373889028291049 && ...
   diagnostic.sharedBeta==0.00048479871747565857,...
   'PKPD:DoseBetaReference','Archived D values must remain fixed.');
  assert(isequal(diagnostic.g,results.model_d_growth) && ...
   isequal(diagnostic.cfg,results.config),...
   'PKPD:DoseBetaFixed','D reference PK and growth changed.');
  doseSummary=diagnostic.summary;
  assert(height(doseSummary)==3 && size(diagnostic.allRuns,1)==3*nStarts,...
   'PKPD:DoseBetaCount','Missing dose-specific beta fits.');
  savedInputs=load(fullfile(out,'run_inputs.mat'),'cs');
  assert(isequal(diagnostic.cs,savedInputs.cs),'PKPD:DoseBetaData','Data changed.');
  assert(abs(sum(doseSummary.shared_scaled_SSE)-0.9939692549524157)<1e-9,...
   'PKPD:DoseBetaBaseline','Archived shared-beta SSE was not reproduced.');
  assert(all(doseSummary.dose_scaled_SSE<=doseSummary.shared_scaled_SSE+1e-8),...
   'PKPD:DoseBetaNested','Dose beta fit is worse than its shared-beta starting point.');
  assert(all(abs(doseSummary.dose_scaled_SSE-doseSummary.scan_scaled_SSE)<1e-7),...
   'PKPD:DoseBetaOptimum','Multistart and independent beta scan disagree.');
  vals=diagnostic.predictions(:,4:end);
  assert(all(isfinite(vals(:))) && all(diagnostic.predictions(:,6)>0),...
   'PKPD:DoseBetaCurves','Invalid diagnostic trajectories.');
  for doseId=1:3
   for modelId=1:2
    mask=diagnostic.residuals(:,1)==modelId & ...
     diagnostic.residuals(:,2)==doseSummary.dose_nM(doseId);
    expected=doseSummary.shared_scaled_SSE(doseId);
    if modelId==2,expected=doseSummary.dose_scaled_SSE(doseId);end
    assert(abs(sum(diagnostic.residuals(mask,7).^2)-expected)<1e-10,...
     'PKPD:DoseBetaResiduals','Residuals and objective disagree.');
   end
  end
  checkpointsDose=dir(fullfile(out,'dose_beta_checkpoints','dose_*','*_start_*.mat'));
  assert(numel(checkpointsDose)==3*nStarts,'PKPD:DoseBetaCheckpoints','Missing start checkpoints.');
  status.dose_beta_starts=size(diagnostic.allRuns,1);
  status.dose_beta_shared_SSE=sum(doseSummary.shared_scaled_SSE);
  status.dose_beta_separate_SSE=sum(doseSummary.dose_scaled_SSE);
 end
 status.profile_points=height(profile);
 status.status='success';
 status.completed_starts=numel(checkpoints);
 status.scaled_SSE=summary.scaled_SSE';
 status.max_equivalent_N_difference=results.max_equivalent_N_difference;
 status.elapsed_seconds=toc(timer);
 writeStatus(out,status);
 fprintf('CI OUTPUT CHECKS PASSED: %d starts; %s; %.1f s.\n',...
  expectedStarts,mode,status.elapsed_seconds);
catch err
 status.status='failed';status.elapsed_seconds=toc(timer);
 status.error_identifier=err.identifier;
 status.error_report=getReport(err,'extended','hyperlinks','off');
 writeStatus(out,status);
 fprintf(2,'%s\n',status.error_report);
 rethrow(err);
end
end

function writeStatus(out,status)
fid=fopen(fullfile(out,'ci_status.json'),'w');
assert(fid~=-1,'PKPD:StatusWrite','Cannot write CI status.');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(status));
end
