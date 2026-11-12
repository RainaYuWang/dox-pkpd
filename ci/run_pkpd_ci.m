function results = run_pkpd_ci(mode)
% RUN_PKPD_CI Run the unchanged quick/full fit and check its saved outputs.
if nargin<1,mode='quick';end
mode=validatestring(mode,{'quick','full','sensitivity','dose-beta','saturation'});
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
doSaturation=strcmp(mode,'saturation');
doDoseBeta=strcmp(mode,'dose-beta') || doSaturation;
nStarts=2;if ~strcmp(mode,'quick'),nStarts=8;end
status=struct('status','running','mode',mode,'nStarts',nStarts,...
 'seed',42,'matlab_version',version,'platform',computer,...
 'commit',getenv('GITHUB_SHA'),'run_id',getenv('GITHUB_RUN_ID'),...
 'started_utc',char(datetime('now','TimeZone','UTC')),'elapsed_seconds',0);
status.command=sprintf('run_joint_pkpd(%d,%d,%d,%d)',nStarts,doSensitivity,doDoseBeta,doSaturation);
status.rng_algorithm='twister';
writeStatus(out,status);
try
 results=run_joint_pkpd(nStarts,doSensitivity,doDoseBeta,doSaturation);
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
  % Output agreement tolerance, NOT a change to ODE or optimization tolerances.
  % The first run found a 5.52e-7 SSE difference in the shallow low-dose minimum.
  comparisonTolerance=1e-6;
  assert(all(abs(doseSummary.dose_scaled_SSE-doseSummary.scan_scaled_SSE)<comparisonTolerance) && ...
   all(abs(log(doseSummary.dose_beta./doseSummary.scan_beta))<0.001),...
   'PKPD:DoseBetaOptimum','Multistart and independent beta scan disagree.');
  assert(all(doseSummary.dose_scaled_SSE-doseSummary.local_audit_scaled_SSE<comparisonTolerance),...
   'PKPD:DoseBetaLocalAudit','Local objective audit found a materially lower SSE.');
  status.dose_beta_comparison_tolerance=comparisonTolerance;
  status.dose_beta_scan_SSE_difference=abs(doseSummary.dose_scaled_SSE-doseSummary.scan_scaled_SSE)';
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
 if doSaturation
  sat=load(fullfile(out,'model_d_saturation.mat'));
  assert(sat.kBF==diagnostic.kBF && sat.sharedBeta==diagnostic.sharedBeta && ...
   isequal(sat.g,diagnostic.g) && isequal(sat.cfg,diagnostic.cfg) && isequal(sat.cs,diagnostic.cs),...
   'PKPD:SaturationFixed','Reference PK, growth, kBF or data changed.');
  assert(size(sat.runs,1)==nStarts+1 && height(sat.profile)>=100,...
   'PKPD:SaturationCount','Incomplete multistart/profile.');
  checkpointsSat=dir(fullfile(out,'saturation_checkpoints','*_start_*.mat'));
  assert(numel(checkpointsSat)==nStarts+1,'PKPD:SaturationCheckpoints','Missing start checkpoints.');
  assert(all(isfinite(sat.profile.scaled_SSE)) && all(sat.profile.scaled_SSE<1e9),...
   'PKPD:SaturationProfile','Invalid profile values.');
  assert(abs(min(sat.profile.scaled_SSE)-sat.parameters.scaled_SSE)<1e-5,...
   'PKPD:SaturationOptimum','Joint fit and Ks profile disagree.');
  assert(sat.parameters.scaled_SSE<=status.dose_beta_shared_SSE+1e-5,...
   'PKPD:SaturationNested','Saturation fit is worse than linear limit.');
  assert(abs(sat.summary.scaled_SSE(1)-status.dose_beta_shared_SSE)<1e-10 && ...
   abs(sat.summary.scaled_SSE(2)-status.dose_beta_separate_SSE)<1e-10,...
   'PKPD:SaturationComparators','Original comparators changed.');
  vals=sat.predictions(:,4:end);
  assert(all(isfinite(vals(:))) && all(sat.predictions(:,6)>0),...
   'PKPD:SaturationCurves','Invalid saturation predictions.');
  for m=1:3
   mask=sat.residuals(:,1)==m;
   assert(abs(sum(sat.residuals(mask,7).^2)-sat.summary.scaled_SSE(m))<1e-10,...
    'PKPD:SaturationResiduals','SSE and residuals disagree.');
  end
  exact=sat.limits(isinf(sat.limits(:,1)),:);
  large=sat.limits(sat.limits(:,1)==1e8,:);
  assert(all(exact(:,4)<1e-10) && all(exact(:,5)<1e-10),...
   'PKPD:SaturationExactLimit','Ks=Inf does not recover the original linear ODE solution.');
  assert(all(large(:,4)<1e-4) && all(large(:,5)<1e-5),...
   'PKPD:SaturationLargeLimit','Large finite Ks has not approached the linear solution.');
  linearProfile=sat.profile.scaled_SSE(isinf(sat.profile.Ks_nM));
  assert(abs(linearProfile-status.dose_beta_shared_SSE)<1e-5,...
   'PKPD:SaturationLinearProfile','Re-optimized beta at Ks=Inf disagrees with linear D.');
  status.saturation_starts=size(sat.runs,1);
  status.saturation_profile_points=height(sat.profile);
  status.saturation_beta=sat.parameters.beta;status.saturation_Ks=sat.parameters.Ks_nM;
  status.saturation_SSE=sat.parameters.scaled_SSE;
  status.saturation_large_Ks_max_N_relative_difference=max(large(:,4));
  status.saturation_exact_limit_max_N_relative_difference=max(exact(:,4));
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
