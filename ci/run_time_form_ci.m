function results=run_time_form_ci(mode)
% Separate entry: only new A/B diagnostic fits, plus archived comparisons.
mode=validatestring(mode,{'time-quick','time-full'});
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
out=fullfile(root,'joint_results');
assert(~exist(out,'dir'),'PKPD:ExistingResults','Move existing joint_results aside before running.');
mkdir(out);timer=tic;nStarts=2;if strcmp(mode,'time-full'),nStarts=8;end
status=struct('status','running','mode',mode,'nStarts',nStarts,...
 'matlab_version',version,'platform',computer,'commit',getenv('GITHUB_SHA'),...
 'run_id',getenv('GITHUB_RUN_ID'),'started_utc',char(datetime('now','TimeZone','UTC')),...
 'elapsed_seconds',0,'command',sprintf('run_time_form_diagnostics(%d)',nStarts),...
 'rng_algorithm','twister','old_models_refitted',false,'seed_rule','57724+100*form+dose_index');
writeStatus(out,status);
try
 results=run_time_form_diagnostics(nStarts);
 summary=results.summary;params=results.parameters;starts=results.starts;
 assert(height(params)==6 && height(starts)==6*nStarts,'PKPD:TimeCounts','Incomplete diagnostic fits.');
 checkpoints=dir(fullfile(out,'time_form_checkpoints','form_*','*_start_*.mat'));
 assert(numel(checkpoints)==6*nStarts,'PKPD:TimeCheckpoints','Missing per-start checkpoints.');
 assert(all(isfinite(params.scaled_SSE)) && all(params.scaled_SSE<1e9) && ...
  all(params.scaled_SSE<=params.zero_killing_SSE+1e-10),...
  'PKPD:TimeObjective','Invalid fit or worse than included zero-killing candidate.');
 assert(isequal(summary.parameters_refitted_this_comparison,[0;0;3;6;0;0]),...
  'PKPD:TimeFreeParameters','Old models must not be re-fitted.');
 baseline=readtable(fullfile(out,'baseline_delay_summary.csv'));
 for m=[1,2,5,6]
  assert(abs(summary.scaled_SSE(m)-baseline.scaled_SSE(baseline.model_id==m))<1e-10,...
   'PKPD:TimeBaseline','Archived baseline SSE changed.');
 end
 curves=readtable(fullfile(out,'time_form_predictions.csv'));
 residuals=readtable(fullfile(out,'time_form_residuals.csv'));
 assert(height(residuals)==32*6 && all(isfinite(curves.N)) && all(curves.N>0) && ...
  all(isfinite(curves.killing_rate)) && all(curves.killing_rate>=0),...
  'PKPD:TimeOutputs','Missing observations or invalid/nonpositive output.');
 for m=1:6
  assert(abs(sum(residuals.scaled_residual(residuals.model_id==m).^2)-summary.scaled_SSE(m))<1e-10,...
   'PKPD:TimeResiduals','SSE disagrees with saved residuals.');
 end
 new=ismember(curves.model_id,[3,4]);
 assert(all(curves.killing_rate(new & curves.tau_h<0)==0),...
  'PKPD:TimeBeforeOnset','Killing before the verified onset is not zero.');
 for file={'time_form_comparison.png','time_form_killing.png','time_form_numerical_checks.csv',...
   'time_form_peak_checks.csv','time_form_clock_checks.csv','time_form_fit.mat','time_form_inputs.mat'}
  info=dir(fullfile(out,file{1}));assert(~isempty(info) && info.bytes>0,'PKPD:TimeMissing','Missing %s',file{1});
 end
 status.status='success';status.scaled_SSE=summary.scaled_SSE';
 status.time_form_completed_starts=height(starts);status.paper_definitions_verified=results.paper.verified;
 status.baseline_run_id=results.baseline.run_id;status.time_form_numerical_checks='passed';
 status.elapsed_seconds=toc(timer);writeStatus(out,status);
 fprintf('TIME-FORM CI PASSED: %s; %d starts; %.1f s.\n',mode,height(starts),toc(timer));
catch err
 status.status='failed';status.elapsed_seconds=toc(timer);status.error_identifier=err.identifier;
 status.error_report=getReport(err,'extended','hyperlinks','off');writeStatus(out,status);
 fprintf(2,'%s\n',status.error_report);rethrow(err);
end
end

function writeStatus(out,status)
fid=fopen(fullfile(out,'ci_status.json'),'w');assert(fid~=-1,'PKPD:StatusWrite','Cannot write status.');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(status));
end
