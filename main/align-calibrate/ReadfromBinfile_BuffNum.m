function [I,n,VarName,VarValue,StartShift,Ver] = ReadfromBinfile_BuffNum(Data_Path,LinesPerScan,frame,NoOfAuxData)

if nargin<4
    NoOfAuxData = 0;
end

if nargin<3
    frame = 0;
end

if nargin<2
    LinesPerScan = inf;
end

fid = fopen(Data_Path,'r');
numc = fread(fid,1,'uint32');
Ver = 'unidentified';
if numc < 50
    Ver = char(fread(fid,numc,'char')');
    numc = fread(fid,1,'uint32');
    type = char(fread(fid,numc,'char')');
end

% % % if strcmpi(Ver,'CloseLoop_Ver11')|| strcmpi(Ver,'OCT_Ver11')
% % %     BeginingPos = ftell(fid);
% % %     HrzOffsetDetected = 0;
% % %     txt1 = ', Horiz. Offset:';
% % %     txt2 = ', Vert. Offset:';
% % %     while (~HrzOffsetDetected)
% % %         txt1Detected = 0;
% % %         txt2Detected = 0;
% % %         i = 0;
% % %         numc = length(txt1);
% % %         while (~txt1Detected)
% % %             fseek(fid,BeginingPos+i,'bof');
% % %             tmp = char(fread(fid,numc,'char')');
% % %             if strcmpi(tmp,txt1)
% % %                 txt1Detected = 1;
% % %                 Pos1 = ftell(fid) -  length(txt1);
% % %             end
% % %             i = i + 1;
% % %         end
% % %         while (~txt2Detected)
% % %             fseek(fid,BeginingPos+i,'bof');
% % %             tmp = char(fread(fid,numc,'char')');
% % %             if strcmpi(tmp,txt2)
% % %                 txt2Detected = 1;
% % %                 Pos2 = ftell(fid) -  length(txt2);
% % %             end
% % %             i = i + 1;
% % %         end
% % %         if (Pos2-Pos1)>8
% % %             Ver = [Ver,'_MultiVess'];
% % %         end
% % %         HrzOffsetDetected = 1;
% % %     end
% % % end        
% % % 
% % % fseek(fid,BeginingPos,'bof');
if strcmpi(Ver,'CloseLoop_Ver11')|| strcmpi(Ver,'OCT_Ver11')
    numc = fread(fid,1,'uint32');
    temp = char(fread(fid,numc,'char')');
    TotalVars = fread(fid,1,'double');
    %     disp([temp,' ' ,num2str(TotalVars)]);
    for VarCounter = 1:TotalVars
        numc = fread(fid,1,'uint32');
        VarName{VarCounter} = char(fread(fid,numc,'char')');
        VarName{VarCounter} = VarName{VarCounter}(3:end-1);
        VarValue{VarCounter} = fread(fid,1,'double');
        %         disp([VarName{VarCounter},' ' ,num2str(VarValue{VarCounter})])
    end
    numc = fread(fid,1,'uint32');
    DataNames = char(fread(fid,numc,'char')');
    numc = fread(fid,1,'uint32');
    DataBegins= char(fread(fid,numc,'char')');
    StartShift = ftell(fid);
    
elseif strcmpi(Ver,'OCT_Ver15')
    
    numc = fread(fid,1,'uint32');
    temp = char(fread(fid,numc,'char')');
    TotalVars = fread(fid,1,'double');
    %     disp([temp,' ' ,num2str(TotalVars)]);
    for VarCounter = 1:TotalVars
        numc = fread(fid,1,'uint32');
        VarName{VarCounter} = char(fread(fid,numc,'char')');
        VarName{VarCounter} = VarName{VarCounter}(3:end-1);
        if strcmpi(VarName{VarCounter},'Horiz. Offset') || strcmpi(VarName{VarCounter},'Vert. Offset') || strcmpi(VarName{VarCounter},'Angle') 
            VarValue{VarCounter} = fread(fid,1,'int32');
            VarValue{VarCounter} = fread(fid,VarValue{VarCounter},'double');
        else
            VarValue{VarCounter} = fread(fid,1,'double');
        end
        %         disp([VarName{VarCounter},' ' ,num2str(VarValue{VarCounter})])
    end
    numc = fread(fid,1,'uint32');
    DataNames = char(fread(fid,numc,'char')');
    numc = fread(fid,1,'uint32');
    DataBegins= char(fread(fid,numc,'char')');
    StartShift = ftell(fid);
    
elseif strcmpi(Ver,'OCT_Ver17')
    
    numc = fread(fid,1,'uint32');
    temp = char(fread(fid,numc,'char')');
    TotalVars = fread(fid,1,'double');
    %     disp([temp,' ' ,num2str(TotalVars)]);
    for VarCounter = 1:TotalVars
        numc = fread(fid,1,'uint32');
        VarName{VarCounter} = char(fread(fid,numc,'char')');
        VarName{VarCounter} = VarName{VarCounter}(3:end-1);
        if strcmpi(VarName{VarCounter},'Horiz. Offset') || strcmpi(VarName{VarCounter},'Vert. Offset') || strcmpi(VarName{VarCounter},'Angle') 
            VarValue{VarCounter} = fread(fid,1,'int32');
            VarValue{VarCounter} = fread(fid,VarValue{VarCounter},'double');
        else
            VarValue{VarCounter} = fread(fid,1,'double');
        end
        %         disp([VarName{VarCounter},' ' ,num2str(VarValue{VarCounter})])
    end
    numc = fread(fid,1,'uint32');
    DataNames = char(fread(fid,numc,'char')');
    numc = fread(fid,1,'uint32');
    DataBegins= char(fread(fid,numc,'char')');
    StartShift = ftell(fid);
    
elseif strcmpi(Ver(1:8),'oct_ver6')
    
    numc = fread(fid,1,'uint32');
    temp = char(fread(fid,numc,'char')');
    TotalVars = fread(fid,1,'double');
    disp([temp,' ' ,num2str(TotalVars)]);
    for VarCounter = 1:TotalVars
        numc = fread(fid,1,'uint32');
        VarName{VarCounter} = char(fread(fid,numc,'char')');
        VarValue{VarCounter} = fread(fid,1,'double');
%         disp([VarName{VarCounter},' ' ,num2str(VarValue{VarCounter})])
    end
    numc = fread(fid,1,'uint32');
    DataBegins= char(fread(fid,numc,'char')');
    StartShift = ftell(fid);
else
    Start = 0;
    StartShift = 0;
    fseek(fid,Start,'bof');
    VarName = [];
    VarValue = [];
end

CamPix = 1024;
Start = (LinesPerScan*(1024*2)+NoOfAuxData*8)*(frame-1) + StartShift;
End = Start + LinesPerScan*(1024*2);
fseek(fid,Start,'bof');
I = [];
n = [];
if isinf(End)
    I = fread(fid,'int16');
    n = [];
else
    I = fread(fid,(End-Start)/2,'int16');
    for nCounter = 1:NoOfAuxData
        n{nCounter} = fread(fid,1,'double');
    end
end
fclose(fid);
I = reshape(I,CamPix,[])';