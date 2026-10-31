% Run only synthetic experiments from this repository.
projectRoot = fileparts(mfilename('fullpath'));
originalDirectory = pwd;
restoreDirectory = onCleanup(@() cd(originalDirectory));
cd(projectRoot);
if ~exist('fig', 'dir')
    mkdir('fig');
end
rng(1, 'twister');
run('task1_3.m');
run('task7.m');
run('advanced.m');
fprintf('Synthetic experiments complete. Figures saved in fig/.\n');
