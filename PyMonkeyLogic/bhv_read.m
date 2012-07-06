function BHV = bhv_read(varargin)
% SYNTAX:
%          BHV = bhvread(datafile)
%
% Reads behavior data files created by MonkeyLogic.
%
% File format is:
%
%   BHV.MagicNumber                'uint32' x 1
%   BHV.FileHeader                 'uchar'  x 64
%   BHV.FileVersion                'double' x 1
%   BHV.StartTime                  'uchar'  x 32
%   BHV.ExperimentName             'uchar'  x 128
%   BHV.SubjectName                'uchar'  x 128
%   BHV.ComputerName               'uchar'  x 128
%   BHV.ConditionsFile             'uchar'  x 128
%   BHV.NumConds                   'uint16' x 1
%   BHV.ObjectsPerCond             'uint16' x 1
%   BHV.TaskObject                 'uchar'  x 64 x BHV.NumConds x BHV.ObjectsPerCond
%   BHV.NumTimingFiles             'uint8'  x 1
%   BHV.TimingFiles                'uchar'  x 128  x BHV.NumTimingFiles
%   BHV.ErrorLogic                 'uchar'  x 64
%   BHV.BlockLogic                 'uchar'  x 64
%   BHV.CondLogic                  'uchar'  x 64
%   BHV.BlockSelectFunction        'uchar'  x 64
%   BHV.CondSelectFunction         'uchar'  x 64
%   BHV.VideoRefreshRate           'double' x 1
%   BHV.VideoBufferPages           'uint16' x 1
%   BHV.ScreenXresolution          'uint16' x 1
%   BHV.ScreenYresolution          'uint16' x 1
%   BHV.ViewingDistance            'double' x 1
%   BHV.PixelsPerDegree            'double' x 1
%   BHV.AnalogInputType            'uchar'  x 32
%   BHV.AnalogInputFrequency       'double' x 1
%   BHV.AnalogInputDuplication     'uchar'  x 32
%   BHV.InputCalibrationMethod     'uchar'  x 32
%   BHV.PhotoDiode                 'uchar'  x 12
%   BHV.ScreenBackgroundColor      'double' x 3
%   BHV.EyeTraceColor              'double' x 3
%   BHV.JoyTraceColor              'double' x 3
%   BHV.Stimuli:
%       BHV.Stimuli.NumPics        'uint16' x 1
%       BHV.Stimuli.PIC.Name       'uchar'  x 128  x BHV.Stimuli.NumPics
%       BHV.Stimuli.PIC.Size       'uint16' x 3    x BHV.Stimuli.NumPics
%       BHV.Stimuli.PIC.Data       'uint8'  x sum-over-i(prod(BHV.Stimuli.PIC.Size(i)))
%       BHV.Stimuli.NumMovs        'uint16' x 1
%       BHV.Stimuli.MOV.Name       'uchar'  x 128  x BHV.Stimuli.NumMovs
%       BHV.Stimuli.MOV.Size       'uint16' x 3    x BHV.Stimuli.NumMovs
%       BHV.Stimuli.MOV.NumFrames  'uint16' x 1    x BHV.Stimuli.NumMovs
%       BHV.Stimuli.MOV.Data       'uint8'  x sum-over-i(prod(BHV.Stimuli.MOV.Size(i)))
%   BHV.Padding                    'uint8'  x 1024
%   BHV.NumTrials                  'uint16' x 1
%       .--------------------------------------------------------------
%       | BHV.TrialNumber(trial)                 'uint16' x 1
%       | n (length of following field)          'uint8' x 1
%       | BHV.AbsoluteTrialStartTime(trial, 1:n) 'double' x 1 x n
%       | BHV.BlockNumber(trial)                 'uint16' x 1
%       | BHV.CondNumber(trial)                  'uint16' x 1
%       | BHV.TrialError(trial)                  'uint16' x 1
%       | BHV.CycleRate(trial)                   'uint16' x 1
%       | BHV.NumCodes(trial)                    'uint16' x 1
%       | BHV.CodeNumbers{trial}                 'uint16' x 1 x BHV.NumCodes
%       | BHV.CodeTimes{trial}                   'uint16' x 1 x BHV.NumCodes
%       | BHV.NumXEyePoints(trial)               'uint32' x 1
%       | BHV.AnalogData{trial}.EyeSignal[x]     'float32' x BHV.NumXEyePoints
%       | BHV.NumYEyePoints(trial)               'uint32' x 1
%       | BHV.AnalogData{trial}.EyeSignal[y]     'float32' x BHV.NumYEyePoints
%       | BHV.NumXJoyPoints(trial)               'uint32' x 1
%       | BHV.AnalogData{trial}.Joystick [x]     'float32' x BHV.NumXJoyPoints
%       | BHV.NumYJoyPoints(trial)               'uint32' x 1
%       | BHV.AnalogData{trial}.Joystick [y]     'float32' x BHV.NumYJoyPoints
%       | BHV.ReactionTime(trial)                'uint16' x 1
%       | BHV.ObjectStatusRecord(trial).Number   'uint32' x 1
%       | BHV.ObjectStatusRecord(trial).Status   'ubit1' x ...Number
%       | BHV.ObjectStatusRecord(trial).Time     'uint32' x 1
%       .--------------------------------------------------------------
%   BHV.NumBehavioralCodesUsed     'uint16' x 1
%   BHV.CodeNumbersUsed            'uint16' x 1   x BHV.NumBehavioralCodes
%   BHV.CodeNamesUsed              'uchar'  x 64  x BHV.NumBehavioralCodes
%   BHV.FinishTime                 'uchar'  x 32
%
% Created by WA 7/06
% Modified 5/22/08 -WA (added BlockIndex)
% Modified 8/13/08 -WA (added Movies)

BHV = struct;
if ~ispref('MonkeyLogic', 'Directories'),
    success = set_ml_preferences;
    if ~success,
        return
    end
end
MLDATA.Directories = getpref('MonkeyLogic', 'Directories');
expdir = MLDATA.Directories.ExperimentDirectory;

if isempty(varargin),
    [datafile pathname] = uigetfile([expdir '*.bhv'], 'Choose BHV file');
    if datafile == 0,
        return
    end
    datafile = [pathname datafile];
else
    datafile = varargin{1};
end

[pname fname ext] = fileparts(datafile);
if isempty(pname),
    pname = expdir;
else
    pname = [pname filesep];
end
if isempty(ext),
    ext = '.bhv';
end
datafile = [pname fname ext];
BHV.DataFileName = fname;
BHV.FullDataFile = datafile;

fidbhv = fopen(datafile, 'r');
if fidbhv == -1,
    error(sprintf('*** Unable to open data file: %s', datafile));
end

BHV.MagicNumber = fread(fidbhv, 1, 'uint32');
if BHV.MagicNumber ~= 2837160,
    error('*** %s is not recognized as a "BHV" file ***', datafile);
end
BHV.FileHeader = deblank(char(fread(fidbhv, 64, 'uchar')'));
BHV.FileVersion = fread(fidbhv, 1, 'double');
if BHV.FileVersion < 1.5,
    error('*** BHV files of version < 1.5 are no longer supported ***');
end
BHV.StartTime = deblank(char(fread(fidbhv, 32, 'uchar')'));
BHV.ExperimentName = deblank(char(fread(fidbhv, 128, 'uchar')'));
if BHV.FileVersion > 1.5,
    BHV.Investigator = deblank(char(fread(fidbhv, 128, 'uchar')'));
end
BHV.SubjectName = deblank(char(fread(fidbhv, 128, 'uchar')'));
if BHV.FileVersion > 2.1,
    BHV.ComputerName = deblank(char(fread(fidbhv, 128, 'uchar')'));
end
BHV.ConditionsFile = deblank(char(fread(fidbhv, 128, 'uchar')'));
BHV.NumConds = fread(fidbhv, 1, 'uint16');
BHV.ObjectsPerCond = fread(fidbhv, 1, 'uint16');
BHV.TaskObject = fread(fidbhv, BHV.NumConds * BHV.ObjectsPerCond * 64, 'uchar');
BHV.TaskObject = deblank(cellstr(char(reshape(BHV.TaskObject, 64, BHV.NumConds * BHV.ObjectsPerCond)')));
BHV.TaskObject = reshape(BHV.TaskObject, BHV.NumConds, BHV.ObjectsPerCond);
if BHV.FileVersion >= 2.65,
    BHV.TimingFileByCond = fread(fidbhv, BHV.NumConds * 128, 'uchar');
    BHV.TimingFileByCond = deblank(cellstr(char(reshape(BHV.TimingFileByCond, 128, BHV.NumConds)')));
    if BHV.FileVersion >= 2.71,
        maxblocks = fread(fidbhv, 1, 'uint8');
        bbc = fread(fidbhv, [BHV.NumConds maxblocks], 'uint8');
        for i = 1:BHV.NumConds,
            bbci = bbc(i,:);
            bbci = bbci(bbci ~= 0);
            BHV.BlockByCond{i} = bbci;
        end
    else
        BHV.BlockByCond = fread(fidbhv, BHV.NumConds, 'uint8');
    end
    BHV.InfoByCond = fread(fidbhv, BHV.NumConds * 128, 'uchar');
    BHV.InfoByCond = deblank(cellstr(char(reshape(BHV.InfoByCond, 128, BHV.NumConds)')));
    for i = 1:length(BHV.InfoByCond),
        BHV.InfoByCond{i} = eval(sprintf('struct(%s)',BHV.InfoByCond{i}));
    end
end
    
BHV.NumTimingFiles = fread(fidbhv, 1, 'uint8');
for i = 1:BHV.NumTimingFiles,
    BHV.TimingFiles{i} = deblank(char(fread(fidbhv, 128, 'uchar')'));
end
BHV.ErrorLogic = deblank(char(fread(fidbhv, 64, 'uchar')'));
BHV.BlockLogic = deblank(char(fread(fidbhv, 64, 'uchar')'));
BHV.CondLogic = deblank(char(fread(fidbhv, 64, 'uchar')'));
BHV.BlockSelectFunction = deblank(char(fread(fidbhv, 64, 'uchar')'));
BHV.CondSelectFunction = deblank(char(fread(fidbhv, 64, 'uchar')'));
if BHV.FileVersion > 2.0,
    BHV.VideoRefreshRate = fread(fidbhv, 1, 'double');
    BHV.VideoBufferPages = fread(fidbhv, 1, 'uint16');
end
BHV.ScreenXresolution = fread(fidbhv, 1, 'uint16');
BHV.ScreenYresolution = fread(fidbhv, 1, 'uint16');
BHV.ViewingDistance = fread(fidbhv, 1, 'double');
BHV.PixelsPerDegree = fread(fidbhv, 1, 'double');
if BHV.FileVersion > 2.01,
    BHV.AnalogInputType = deblank(char(fread(fidbhv, 32, 'uchar')));
    BHV.AnalogInputFrequency = fread(fidbhv, 1, 'double');
end
if BHV.FileVersion > 2.0,
    BHV.AnalogInputDuplication = deblank(char(fread(fidbhv, 32, 'uchar')));
end
BHV.EyeSignalCalibrationMethod = deblank(char(fread(fidbhv, 32, 'uchar')'));
tmatrix = fread(fidbhv, 1, 'uint8');
if tmatrix,
    BHV.EyeTransform.ndims_in = fread(fidbhv, 1, 'uint16');
    BHV.EyeTransform.ndims_out = fread(fidbhv, 1, 'uint16');
    fwdfcn = deblank(fread(fidbhv, 64, '*char'))';
    eval(['BHV.EyeTransform.forward_fcn = ' fwdfcn ';']);
    invfcn = deblank(fread(fidbhv, 64, '*char'))';
    eval(['BHV.EyeTransform.inverse_fcn = ' invfcn ';']);
    tsize = fread(fidbhv, 1, 'uint16');
    tsqrt = sqrt(tsize);
    BHV.EyeTransform.tdata.T = reshape(fread(fidbhv, tsize, 'double'), tsqrt, tsqrt);
    BHV.EyeTransform.tdata.Tinv = reshape(fread(fidbhv, tsize, 'double'), tsqrt, tsqrt);
else
    BHV.EyeTransform = [];
end
BHV.JoystickCalibrationMethod = deblank(char(fread(fidbhv, 32, 'uchar')'));
tmatrix = fread(fidbhv, 1, 'uint8');
if tmatrix,
    BHV.JoyTransform.ndims_in = fread(fidbhv, 1, 'uint16');
    BHV.JoyTransform.ndims_out = fread(fidbhv, 1, 'uint16');
    fwdfcn = deblank(fread(fidbhv, 64, '*char'))';
    eval(['BHV.JoyTransform.forward_fcn = ' fwdfcn ';']);
    invfcn = deblank(fread(fidbhv, 64, '*char'))';
    eval(['BHV.JoyTransform.inverse_fcn = ' invfcn ';']);
    tsize = fread(fidbhv, 1, 'uint16');
    tsqrt = sqrt(tsize);
    BHV.JoyTransform.tdata.T = reshape(fread(fidbhv, tsize, 'double'), tsqrt, tsqrt);
    BHV.JoyTransform.tdata.Tinv = reshape(fread(fidbhv, tsize, 'double'), tsqrt, tsqrt);
else
    BHV.JoyTransform = [];
end
BHV.PhotoDiodePosition = deblank(char(fread(fidbhv, 12, 'uchar')'));

if BHV.FileVersion >= 1.9,
   BHV.ScreenBackgroundColor = fread(fidbhv, 3, 'double');
   BHV.EyeTraceColor = fread(fidbhv, 3, 'double');
   BHV.JoyTraceColor = fread(fidbhv, 3, 'double');
end

BHV.Stimuli.NumPics = fread(fidbhv, 1, 'uint16');
if BHV.Stimuli.NumPics > 0,
    for i = 1:BHV.Stimuli.NumPics,
        PIC(i).Name = deblank(char(fread(fidbhv, 128, 'uchar')'));
    end
    for i = 1:BHV.Stimuli.NumPics,
        PIC(i).Size = fread(fidbhv, 3, 'uint16');
    end
    for i = 1:BHV.Stimuli.NumPics,
        sz = PIC(i).Size;
        d = fread(fidbhv, prod(sz), 'uint8');
        PIC(i).Data = reshape(d, sz(1), sz(2), sz(3));
    end
    BHV.Stimuli.PIC = PIC;
else
    BHV.Stimuli.PIC = [];
end

if BHV.FileVersion >= 2.5,
    BHV.Stimuli.NumMovs = fread(fidbhv, 1, 'uint16');
    if BHV.Stimuli.NumMovs > 0,
        for i = 1:BHV.Stimuli.NumMovs,
            MOV(i).Name = deblank(char(fread(fidbhv, 128, 'uchar')'));
        end
        for i = 1:BHV.Stimuli.NumMovs,
            sz = fread(fidbhv, 4, 'uint16');
            MOV(i).Size = sz(1:3);
            MOV(i).NumFrames = sz(4);
        end
        for i = 1:BHV.Stimuli.NumMovs,
            sz = MOV(i).Size;
            numframes = MOV(i).NumFrames;
            for fnum = 1:numframes,
                d = fread(fidbhv, prod(sz), 'uint8');
                mframe = reshape(d, sz(1), sz(2), sz(3));
                mframe = flipdim(mframe, 1);
                MOV(i).Data{fnum} = mframe;
            end
        end
        BHV.Stimuli.MOV = MOV;
    else
        BHV.Stimuli.MOV = [];
    end
end

BHV.Padding = fread(fidbhv, 1024, 'uint8');
BHV.NumTrials = fread(fidbhv, 1, 'uint16');
for trial = 1:BHV.NumTrials,
    BHV.TrialNumber(trial, 1) = fread(fidbhv, 1, 'uint16');
    if BHV.FileVersion > 2.2,
        numc = fread(fidbhv, 1, 'uint8');
        BHV.AbsoluteTrialStartTime(trial, 1:numc) = fread(fidbhv, numc, 'double')';
    end
    BHV.BlockNumber(trial, 1) = fread(fidbhv, 1, 'uint16');
    % Create "BlockIndex" - ordinal position of each trial's block:
    bnum = BHV.BlockNumber;
    dblock = find(diff(bnum)) + 1;
    dblock(2:length(dblock)+1) = dblock;
    dblock(1) = 1;
    dblock(length(dblock)+1) = BHV.NumTrials+1;
    BHV.BlockIndex = zeros(BHV.NumTrials, 1);
    for i = 1:length(dblock)-1,
        x1 = dblock(i);
        x2 = dblock(i+1)-1;
        BHV.BlockIndex(x1:x2) = i;
    end
    BHV.ConditionNumber(trial, 1) = fread(fidbhv, 1, 'uint16');
    BHV.TrialError(trial, 1) = fread(fidbhv, 1, 'uint16');
    if BHV.FileVersion >= 2.05,
        BHV.CycleRate(trial, 1) = fread(fidbhv, 1, 'uint16');
        if BHV.FileVersion >= 2.72,
            BHV.MinCycleRate(trial, 1) = fread(fidbhv, 1, 'uint16');
        end
    end
    BHV.NumCodes(trial, 1) = fread(fidbhv, 1, 'uint16');
    BHV.CodeNumbers{trial} = fread(fidbhv, BHV.NumCodes(trial), 'uint16');
    BHV.CodeTimes{trial} = fread(fidbhv, BHV.NumCodes(trial), 'uint16');
    if BHV.FileVersion >= 1.5,
        BHV.NumXEyePoints = fread(fidbhv, 1, 'uint32');
        if BHV.NumXEyePoints > 0,
            if BHV.FileVersion > 1.6,
                xeye = fread(fidbhv, BHV.NumXEyePoints, 'float32');
            else
                xeye = fread(fidbhv, BHV.NumXEyePoints, 'double');
            end
        end
        BHV.NumYEyePoints = fread(fidbhv, 1, 'uint32');
        if BHV.NumYEyePoints > 0,
            if BHV.FileVersion > 1.6,
                yeye = fread(fidbhv, BHV.NumYEyePoints, 'float32');
            else
                yeye = fread(fidbhv, BHV.NumYEyePoints, 'double');
            end
        end
        if BHV.NumXEyePoints == BHV.NumYEyePoints && BHV.NumXEyePoints > 0,
            BHV.AnalogData{trial}.EyeSignal = [xeye yeye];
        elseif BHV.NumXEyePoints > BHV.NumYEyePoints, %only recorded one of the two
            BHV.AnalogData{trial}.EyeSignal = xeye;
        end
        
        BHV.NumXJoyPoints = fread(fidbhv, 1, 'uint32');
        if BHV.NumXJoyPoints > 0,
            if BHV.FileVersion > 1.6,
                xjoy = fread(fidbhv, BHV.NumXJoyPoints, 'float32');
            else
                xjoy = fread(fidbhv, BHV.NumXJoyPoints, 'double');
            end
        end
        BHV.NumYJoyPoints = fread(fidbhv, 1, 'uint32');
        if BHV.NumYJoyPoints > 0,
            if BHV.FileVersion > 1.6,
                yjoy = fread(fidbhv, BHV.NumYJoyPoints, 'float32');
            else
                yjoy = fread(fidbhv, BHV.NumYJoyPoints, 'double');
            end
        end
        if BHV.NumXJoyPoints == BHV.NumYJoyPoints && BHV.NumXJoyPoints > 0,
            BHV.AnalogData{trial}.Joystick = [xjoy yjoy];
        elseif BHV.NumXJoyPoints > BHV.NumYJoyPoints, %only recorded one of the two
            BHV.AnalogData{trial}.Joystick = xjoy;
        end
        
        if BHV.FileVersion > 2.5,
            for i = 1:9
                gname = sprintf('Gen%i', i);
                BHV.NumGenPoints = fread(fidbhv, 1, 'uint32');
                gen = fread(fidbhv, BHV.NumGenPoints, 'float32');
                BHV.AnalogData{trial}.General.(gname) = gen;
            end
        end
        
        if BHV.FileVersion >= 1.8,
            BHV.NumPhotoDiodePoints = fread(fidbhv, 1, 'uint32');
            if BHV.NumPhotoDiodePoints > 0,
                pd = fread(fidbhv, BHV.NumPhotoDiodePoints, 'float32');
            else
                pd = [];
            end
            BHV.AnalogData{trial}.PhotoDiode = pd;
        end
        BHV.ReactionTime(trial) = fread(fidbhv, 1, 'int16');
    end
    BHV.ReactionTime(BHV.ReactionTime == -1) = NaN;
    if BHV.FileVersion >= 1.9,
        if BHV.FileVersion >= 2.00,
            numstat = fread(fidbhv, 1, 'uint32');
            for i = 1:numstat,
                numbits = fread(fidbhv, 1, 'uint32');
                BHV.ObjectStatusRecord(trial).Status{i} = fread(fidbhv, numbits, 'uint8')';
                BHV.ObjectStatusRecord(trial).Time(i) = fread(fidbhv, 1, 'uint32');
                if any(BHV.ObjectStatusRecord(trial).Status{i} > 1),
                    numfields = fread(fidbhv, 1, 'uint8');
                    if numfields,
                        for fnum = 1:numfields,
                            datacount = fread(fidbhv, 1, 'uint32');
                            BHV.ObjectStatusRecord(trial).Data{i, fnum} = fread(fidbhv, datacount, 'double');
                        end
                    else
                        BHV.ObjectStatusRecord(trial).Data{i} = NaN;
                    end
                    %[Status] is a vector, one element per object, that indicates if each object is currently visible (1),
                    %invisible (0), or altered in some other way (2 = static change in position; 3 = Movie frame / translation change).
                    %[Data] contains information regarding movie frame number, on-line changes in object position 
                    %(or possibly orientation, etc, later)
                end
            end
        else
            numstat = fread(fidbhv, 1, 'uint32');
            for i = 1:numstat,
                numbits = fread(fidbhv, 1, 'uint32');
                BHV.ObjectStatusRecord(trial).Status{i} = fread(fidbhv, numbits, 'ubit1')';
                BHV.ObjectStatusRecord(trial).Time(i) = fread(fidbhv, 1, 'uint32');
                BHV.ObjectStatusRecord(trial).Data = {};
            end
        end
    else
        BHV.ObjectStatusRecord = [];
    end
    if BHV.FileVersion >= 1.95,
        numreward = fread(fidbhv, 1, 'uint32');
        BHV.RewardRecord(trial).RewardOnTime = fread(fidbhv, numreward, 'uint32');
        BHV.RewardRecord(trial).RewardOffTime = fread(fidbhv, numreward, 'uint32');
    else
        BHV.RewardRecord = [];
    end
    
    if BHV.FileVersion >= 2.7,
        numUserVars = fread(fidbhv, 1, 'uint8');
        for i = 1:numUserVars,
            varn = deblank(char(fread(fidbhv, 32, 'uchar')'));
            type = char(fread(fidbhv, 1, 'uchar'));
            if type == 'd',
                lenv = fread(fidbhv, 1, 'uint8');
                varv = fread(fidbhv, lenv, 'double');
            elseif type == 'c',
                varv = deblank(char(fread(fidbhv, 128, 'uchar')'));
            else
                varv = [];
            end
            BHV.UserVars(trial).(varn) = varv;
        end
    end
end
bo = cat(1, 0, BHV.BlockNumber);
BHV.BlockOrder = bo(find(diff(bo))+1);
BHV.NumBehavioralCodesUsed = fread(fidbhv, 1, 'uint16');
BHV.CodeNumbersUsed = fread(fidbhv, BHV.NumBehavioralCodesUsed, 'uint16');
for i = 1:BHV.NumBehavioralCodesUsed,
    BHV.CodeNamesUsed{i, 1} = deblank(char(fread(fidbhv, 64, 'uchar')'));
end

if BHV.FileVersion > 2.05,
    VarChanges = struct;
    numf = fread(fidbhv, 1, 'uint16');
    for i = 1:numf,
        fn = deblank(fread(fidbhv, 64, '*char')');
        n = fread(fidbhv, 1, 'uint16');
        VarChanges.(fn).Trial = fread(fidbhv, n, 'uint16');
        VarChanges.(fn).Value = fread(fidbhv, n, 'double');
    end
    BHV.VariableChanges = VarChanges;
end

BHV.FinishTime = deblank(char(fread(fidbhv, 32, 'uchar')'));

fclose(fidbhv);

rmlist = {'NumTimingFiles' 'NumStimuli' 'NumTrials' 'NumBehavioralCodesUsed' 'NumConds' 'ObjectsPerCond' 'Padding' 'NumXEyePoints' 'NumYEyePoints' 'NumXJoyPoints' 'NumYJoyPoints' 'NumGenPoints'};
for i = 1:length(rmlist),
    if isfield(BHV, rmlist{i}),
        BHV = rmfield(BHV, rmlist{i});
    end
end
