%variables:
studyNames={'DeepeyeANT','EEGManyLabs','AOC', 'GCP'};
studySaveFolders={'DeepeyeANT/data','EEGManyLabs/data', 'AOC/data', 'GCP/data'};
study_IDs={'DE','EML', 'AOC', 'GCP'};
sourceFolders={'sourceDE',...
    'sourceML',...
    'sourceAOC',...
    'sourceGCP'}; % %TODO - CORRECT ME: place where actual data is saved (../study/data)

tarFolder='/home/methlab/transfer/'; %TODO- CORRECT ME

studyfolder={'/home/methlab/Desktop/eyepredict_ANT/All_Subjects',...
    '/home/methlab/Desktop/EEGManyLabs/data',...
    '/home/methlab/Desktop/AOC_data',...
    '/home/methlab/Desktop/GCP_data'};

%% make sure the internet conenction is established
waitfor(msgbox('Make sure the Stimulus Computer is connected to the internet (not to EYELINK), then press ok'));
connected=0;
cnt=1;
while connected==0

    try
        url ='https://www.google.de';
        code=webread(url);
        connected=1;
    catch
        connected=0;
    end

    if connected==0
        waitfor(msgbox(['Still no connection, try connecting again, then press ok. ('...
            num2str(3-cnt) ' more attempts)' ]));
    else
        break;
    end

    if cnt==3
        break;
    end
    cnt=cnt+1;

end

%% mount the network drive:
result=32;
cnt=0;
while not(result==0)
    if cnt==3
        waitfor(msgbox('Connection could not be established, wrong password? Please rerun "SaveData" script','Error','error'));
        return;
    end
    dtmp=dir(tarFolder);
    if length(dtmp)>3
        result=0;
        %do nothin, already mounted
    else
        shortname=inputdlg('Please enter your UZH shortname:                       :' ,'Establishing connection to Server');
        %pw=inputdlg('Bitte das UZH passwort angeben','Verbingung zu Server wird hergestellt');
        waitfor(msgbox('Please enter passwords in MATLAB command window'));

        %order= ['sudo mount -t cifs //psyger-stor02.d.uzh.ch/methlab_xfr/ /home/methlab/transfer/ -o "username=' shortname  ',workgroup=UZH"'];
        order= ['sudo mount -t cifs -o username=' shortname  ',workgroup=UZH,uid=$(id -u),gid=$(id -g) //psyger-stor02.d.uzh.ch/methlab_xfr/ /home/methlab/transfer/ -o "username=' shortname  ',workgroup=UZH"'];
        order=cell2mat(order);
        result=system(order);
        cnt=cnt+1;
        if not(result==0)
            if cnt==3
            else
                waitfor(msgbox('Conenction failed, please try again.'));
            end
        end
    end
end
%% get study name:

[index,tf]=listdlg('PromptString',{'Which Study did you perform?'},...
    'SelectionMode','single',...
    'ListString',studyNames,...
    'ListSize',[600, 200] );

performedStudy=studyNames{index};

%% get subjectID

correctID=0;
cnt=0;
while correctID==0

    prompt = {'Enter subject ID:','Repeat subject ID:'};
    if cnt==0
        dlgtitle=' ';
    else
        dlgtitle = 'MISSMATCH - try again';
    end
    dims = [1 50];
    definput = {'',''};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    correctID=strcmp(answer{1},answer{2});
    cnt=cnt+1;
end

ID=answer{1};

%check if there is a folder for the entered subject:

d=dir(cell2mat([studyfolder(index)])); % EXCHANGE WITH STUDYFOLDER - studyfolders{index}

if not(any(strcmp({d.name},ID)))
    waitfor(msgbox('No Folder found for this subject, makle sure you are using the corret ID and rerun SaveData script','Error','error'));
    return
end

%% Transfer data (will only be executed if Folder of subject was found!)

% % % % %mount: methlab_data/transfer
% % % %
% % % % %TODO : CHECK WHAT HAPPENS IF I TRY TO MOUNT AN ALERADY MOUTNED POINT
% % % % %(IN SAME DIR) - or check before if moutned already!
% % % % if not exist([tarFolder filesep study_IDs{index}])
% % % %     system('mount()') %replace with proper command, inclduing hardcoded usrname&pw of dummy account
% % % % end
% % % % %retry once:
% % % % if not exist([tarFolder filesep study_IDs{index}])
% % % %     system('dir') %replace with proper command, inclduing hardcoded usrname&pw of dummy account
% % % % end
% % % % % still not existing? - then something is wrong, break!
% % % % if not exist([tarFolder filesep study_IDs{index}])
% % % %    waitfor(msgbox('Connecting to server failed, Data was NOT saved! Please retry.','Error','error'));
% % % %    return
% % % % end

%copy to mounted directory from soruce.
%waitfor() needed?
dtmp2= dir([tarFolder studySaveFolders{index} filesep ID]);
if not(any(endsWith({dtmp2.name},'.mat')))
    %system(['sudo chmod a+rwx ' tarFolder])
    % system(['sudo chmod -R 777 ' tarFolder])
    if not(exist([tarFolder studySaveFolders{index} filesep ID]))
        system([ 'sudo mkdir ' tarFolder studySaveFolders{index} filesep ID])
    end
    cp_order= ['sudo cp -RT ' studyfolder{index} filesep ID ' ' tarFolder studySaveFolders{index} filesep ID];
    fail=system(cp_order);
    % [SUCCESS,MESSAGE,MESSAGEID] =copyfile([studyfolder{index} filesep ID], [tarFolder studySaveFolders{index} filesep ID],'f')
    %msgbox('Transfer in Progress, please wait');
else
    [indAns,tf]=listdlg('PromptString',{'(Some) Subject Data is already existing on methlab_xfr.' , 'Are sou sure you want to continue?'},...
        'SelectionMode','single',...
        'ListString',{'Yes', 'No'},...
        'ListSize',[600, 200] );

    if indAns==1
        if not(exist([tarFolder studySaveFolders{index} filesep ID]))
            system([ 'mkdir ' tarFolder studySaveFolders{index} filesep ID])
        end
        cp_order= ['cp -RT ' studyfolder{index} filesep ID ' ' tarFolder studySaveFolders{index} filesep ID];
        fail=system(cp_order);
        % [SUCCESS,MESSAGE,MESSAGEID] =copyfile([studyfolder{index} filesep ID], [tarFolder studySaveFolders{index} filesep ID],'f')
        %msgbox('Transfer in Progress, please wait');
    else
        waitfor(msgbox('Saving aborted!','Error','error'));
        return;
    end


end


if not(fail==0)
    waitfor(msgbox('Saving did not work, please try again carefully. If still no sucess, contact Dawid.','Error','error'));
else
    msgbox('Transfer done succesfully! You can now shut down the Stimulus Computer.');
end
