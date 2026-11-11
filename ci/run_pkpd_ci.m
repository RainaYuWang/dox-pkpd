function results = run_pkpd_ci(mode)
% RUN_PKPD_CI Run the unchanged quick/full fit and check its saved outputs.
if nargin<1,mode='quick';end
mode=validatestring(mode,{'quick','full'});
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root);
out=fullfile(root,'joint_results');
% Never mix a previous run's files with the current run or delete user results.
if exist(out,'dir')
 error('PKPD:ExistingResults','Move existing joint_results aside before a new CI run.');
end
mkdir(out);
timer=tic;
nStarts=2;if strcmp(mode,'full'),nStarts=8;end
status=struct('status','running','mode',mode,'nStarts',nStarts,...
 'seed',42,'matlab_version',version,'platform',computer,...
 'commit',getenv('GITHUB_SHA'),'run_id',getenv('GITHUB_RUN_ID'),...
 'started_utc',char(datetime('now','TimeZone','UTC')),'elapsed_seconds',0);
status.command=sprintf('run_joint_pkpd(%d)',nStarts);
status.rng_algorithm='twister';
writeStatus(out,status);
try
 results=run_joint_pkpd(nStarts);
 % These are runtime/output checks, not tests of predictive scientific validity.
 required={'01_cell_fits.png','02_residuals.png','03_PK_reference.png',...
  '04_equivalent_parameters.png','representative_parameters.csv',...
  'growth_parameters.csv','equivalent_parameter_sets.csv',...
  'reversible_combinations.csv','old_pk_multistart.csv',...
  'joint_irreversible_multistart.csv','joint_reversible_multistart.csv',...
  'predictions.csv','residuals.csv','fit_by_dose.csv','search_bounds.csv',...
  'joint_fit.mat','run_inputs.mat','run.log'};
 for j=1:numel(required)
  info=dir(fullfile(out,required{j}));
  assert(~isempty(info) && info.bytes>0,'PKPD:MissingOutput',...
   'Missing or empty output: %s',required{j});
 end
 expectedStarts=2*min(4,nStarts)+2*nStarts+1;
 checkpoints=dir(fullfile(out,'checkpoints','*_start_*.mat'));
 assert(numel(checkpoints)==expectedStarts,'PKPD:Checkpoints',...
  'Expected %d completed starts; found %d.',expectedStarts,numel(checkpoints));
 summary=readtable(fullfile(out,'representative_parameters.csv'));
 assert(height(summary)==3 && all(isfinite(summary.scaled_SSE)) && ...
  all(summary.scaled_SSE>=0) && all(summary.scaled_SSE<1e9),...
  'PKPD:InvalidFit','Invalid fit summary.');
 predictions=readtable(fullfile(out,'predictions.csv'));
 values=predictions{:,{'CF','CB','N'}};
 assert(height(predictions)>0 && all(isfinite(values(:))) && ...
  all(predictions.N>0),'PKPD:InvalidPredictions','Invalid predictions.');
 saved=load(fullfile(out,'joint_fit.mat'),'results');
 assert(isequal(saved.results.summary,results.summary),...
  'PKPD:SavedResults','Saved MAT summary differs from returned result.');
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
