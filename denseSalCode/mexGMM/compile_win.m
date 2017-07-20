% make the mex file
% Original: Jianming Zhang
% 3/22/2016

function compile_win()

% set the values
opts.opencv_include_path    =   'D:\opencv249\opencv\build\include'; % OpenCV include path
opts.opencv_lib_path        =   'D:\opencv249\opencv\build\x64\vc12\lib'; % OpenCV lib path
%please turn to "https://github.com/MingMingCheng" for the required lib
opts.CmCode_include_path    =   'E:\lab\C_C++\saliency-detection\code\SuZhuo\ChengMingming';
opts.CmCode_Lib_path        =   'E:\lab\C_C++\saliency-detection\code\SuZhuo\ChengMingming\Lib';

opts.clean                  =   false; % clean mode
opts.dryrun                 =   false; % dry run mode
opts.verbose                =   1; % output verbosity
opts.debug                  =   false; % enable debug symbols in MEX-files


% Clean
clear functions

if opts.clean
    if opts.verbose > 0
        fprintf('Cleaning all generated files...\n');
    end

    cmd = fullfile('../',['*.' mexext]);
    if opts.verbose > 0, disp(cmd); end
    if ~opts.dryrun, delete(cmd); end

    cmd = fullfile('*.obj');
    if opts.verbose > 0, disp(cmd); end
    if ~opts.dryrun, delete(cmd); end

    return;
end

% compile flags
[cv_cflags,cv_libs] = pkg_config(opts);
mex_flags = sprintf('%s %s', cv_cflags, cv_libs);
if opts.verbose > 1
    mex_flags = ['-v ' mex_flags];    % verbose mex output
end
if opts.debug
    mex_flags = ['-g ' mex_flags];    % debug vs. optimized builds
end
compstr = computer;
is64bit = strcmp(compstr(end-1:end),'64');
if (is64bit)
  mex_flags = ['-largeArrayDims ' mex_flags];
end

% Compile MxArray and BMS
src = 'MxArray.cpp';
   
cmd = sprintf('mex %s -c %s', mex_flags, src);
if opts.verbose > 0, disp(cmd); end
if ~opts.dryrun, eval(cmd); end

cpplis = {'DenseCRF.cpp'};
src = {'GetFilteredLogFg.cpp'};
% src = {'GetFilteredFg.cpp'};
% src = {'GetFBprobFromGMM.cpp'};
%src = {'TestGMMA.cpp'};
% Compile the mex file
for i = 1:numel(src)
    obj = 'MxArray.obj';
    cmd = sprintf('mex %s %s %s %s -outdir ./', mex_flags, src{i},cpplis{1}, obj);
    if opts.verbose > 0, disp(cmd); end
    if ~opts.dryrun, eval(cmd); end
end

end

%
% Helper functions for windows
%
function [cflags,libs] = pkg_config(opts)
    %PKG_CONFIG  constructs OpenCV-related option flags 
    I_path = opts.opencv_include_path;
    L_path = opts.opencv_lib_path;
    l_options = strcat({' -l'}, lib_names(L_path));
    
    %if opts.debug
    %    l_options = strcat(l_options,'d');    % link against debug binaries
    %end
    l_options = [l_options{:}];

    if ~exist(I_path,'dir')
        error('OpenCV include path not found: %s', I_path);
    end
    if ~exist(L_path,'dir')
        error('OpenCV library path not found: %s', L_path);
    end
    
        %%
    %CmCode
    I_pathCm = opts.CmCode_include_path;
    L_pathCm = opts.CmCode_Lib_path;
    libsCm = cell(2,1);
    libsCm{1} = 'CmLib';
    libsCm{2} = 'CmLibd';
    l_optionsCm = strcat({' -l'}, libsCm);
    l_optionsCm = [l_optionsCm{:}];
    
    if ~exist(I_pathCm,'dir')
        error('CmCode include path not found: %s', I_pathCm);
    end
    if ~exist(L_pathCm,'dir')
        error('CmCode library path not found: %s', L_pathCm);
    end
    %%

    cflags = sprintf('-I''%s''', I_path);
    libs = sprintf('-L''%s'' %s', L_path, l_options);
    
    cflagsCm = sprintf('-I''%s''', I_pathCm);
    libsCm = sprintf('-L''%s'' %s', L_pathCm, l_optionsCm);
    
    cflags = [cflags, ' ', cflagsCm];
    libs = [libs, ' ', libsCm];
end

function l = lib_names(L_path)
    %LIB_NAMES  return library names
    d = dir( fullfile(L_path,'opencv_*.lib') );
    l = regexp({d.name}, '(opencv_core.+)\.lib|(opencv_imgproc.+)\.lib|(opencv_highgui.+)\.lib', 'tokens', 'once');
    l = [l{:}];
end