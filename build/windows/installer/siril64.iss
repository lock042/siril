; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

#include "version.isi"

#define MyAppName "Siril"
#define MyAppExeName "siril.exe"
#define RootDir ROOTDIR

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={ADA3C347-68C3-4EAA-92B3-C1BDBD836EDB}
AppName=Siril
AppVersion={#MAJOR}.{#MINOR}.{#MICRO}
AppPublisher=Free-Astro
AppPublisherURL=https://www.siril.org/
AppSupportURL=https://www.siril.org/
AppUpdatesURL=https://www.siril.org/
DefaultDirName={commonpf}\Siril
DefaultGroupName=Siril
OutputDir={#OUTPUT}
OutputBaseFilename=siril-{#MAJOR}.{#MINOR}.{#MICRO}-setup
Compression=lzma
SolidCompression=yes
ChangesAssociations=yes
ArchitecturesInstallIn64BitMode=x64
DisableDirPage=yes

WizardImageFile=windows-installer-intro-big.bmp
WizardImageStretch=yes
WizardSmallImageFile=siril.bmp

LicenseFile=gpl-3.0.rtf

UninstallDisplayIcon={app}\bin\{#MyAppExeName}

[Languages]
Name: "en"; MessagesFile: "compiler:Default.isl";
Name: "fr"; MessagesFile: "compiler:Languages\French.isl";

[Tasks]
Name: desktopicon; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "{#RootDir}\siril\bin\siril.exe"; DestDir: "{app}\bin"; Flags: ignoreversion
Source: "{#RootDir}\siril\bin\*"; DestDir: "{app}\bin"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "{#RootDir}\siril\lib\*"; DestDir: "{app}\lib"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "{#RootDir}\siril\share\*"; DestDir: "{app}\share"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "{#RootDir}\scripts\*.ssf"; DestDir: "{app}\scripts"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "{#RootDir}\3rdparty\scripts\en\*.ssf"; DestDir: "{app}\scripts"; Flags: ignoreversion recursesubdirs createallsubdirs

; NOTE: Don't use "Flags: ignoreversion" on any shared system files

[Registry]
Root: HKCR; Subkey: ".seq";                             ValueData: "{#MyAppName}";          Flags: uninsdeletevalue; ValueType: string;  ValueName: ""
Root: HKCR; Subkey: "{#MyAppName}";                     ValueData: "Program {#MyAppName}";  Flags: uninsdeletekey;   ValueType: string;  ValueName: ""
Root: HKCR; Subkey: "{#MyAppName}\DefaultIcon";         ValueData: "{app}\bin\{#MyAppExeName},1";               ValueType: string;  ValueName: ""
Root: HKCR; Subkey: "{#MyAppName}\shell\open\command";  ValueData: """{app}\bin\{#MyAppExeName}"" ""%1""";  ValueType: string;  ValueName: ""

[Icons]
Name: "{group}\Siril"; Filename: "{app}\bin\{#MyAppExeName}";
Name: "{group}\{cm:UninstallProgram,Siril}"; Filename: "{uninstallexe}"
Name: "{commondesktop}\Siril"; Filename: "{app}\bin\{#MyAppExeName}"; Tasks: desktopicon;

[Run]
Filename: "{app}\bin\siril.exe"; Description: "{cm:LaunchProgram,Siril}"; Flags: postinstall waituntilidle skipifsilent
Filename: "{code:getFSURL}"; Description: "{code:getFSstring}"; Flags: postinstall nowait shellexec unchecked

[Code]
function getFSURL(s : string) : string;
    var langage : string;
begin
    case ActiveLanguage() of  
        'en' : langage := 'https://siril.org/tutorials/first-steps/';
        'fr' : langage := 'https://siril.org/fr/tutorials/first-steps/';
    end;
    Result := langage;
end;

function getFSstring(s : string) : string;
    var langage : string;
begin
    case ActiveLanguage() of  
        'en' : langage := 'Visit our First Steps page';
        'fr' : langage := 'Ouvrir la page Premiers Pas';
    end;
    Result := langage;
end;

procedure OpenBrowser(Url: string);
var
  ErrorCode: Integer;
begin
  ShellExec('open', Url, '', '', SW_SHOWNORMAL, ewNoWait, ErrorCode);
end;

procedure DonateClick(Sender: TObject);
begin
    case ActiveLanguage() of  
        'en' : OpenBrowser('https://siril.org/donate/');
        'fr' : OpenBrowser('https://siril.org/fr/donate/');
    end;
end;

procedure GetNewsGenericURL(Sender: TObject);
begin
    case ActiveLanguage() of  
        'en' : OpenBrowser('https://siril.org/download/' + '{#SetupSetting("AppVersion")}' + '/');
        'fr' : OpenBrowser('https://siril.org/fr/download/' + '{#SetupSetting("AppVersion")}' + '/');
    end;
end;

procedure CurPageChanged(CurPageID: Integer);
var
  Button: TButton;
  Button2: TButton;
begin
  if CurPageID = wpFinished then
    begin
      Button := TButton.Create(WizardForm);
      Button.Parent := WizardForm;
      Button.Left := ScaleX(16);
      Button.Top := WizardForm.NextButton.Top;
      Button.Width := WizardForm.NextButton.Width;
      Button.Height := WizardForm.NextButton.Height;
      case ActiveLanguage() of  
          'en' : Button.Caption := 'Donate';
          'fr' : Button.Caption := 'Dons';
      end;
      Button.OnClick := @DonateClick;

      Button2 := TButton.Create(WizardForm);
      Button2.Parent := WizardForm;
      Button2.Left := Button.Left + Button.Width + ScaleX(5);
      Button2.Top := WizardForm.NextButton.Top;
      Button2.Width := WizardForm.NextButton.Width;
      Button2.Height := WizardForm.NextButton.Height;
      case ActiveLanguage() of  
          'en' : Button2.Caption := 'News';
          'fr' : Button2.Caption := 'News';
      end;
      Button2.OnClick := @GetNewsGenericURL;
    end;

procedure CurStepChanged(CurStep: TSetupStep);
var
  ResultCode: Integer;
  Uninstall: String;
begin
  if (CurStep = ssInstall) then
    begin
      if RegQueryStringValue(HKLM, 'SOFTWARE\Microsoft\Windows\CurrentVersion\Uninstall\{#AppId}_is1', 'UninstallString', Uninstall) then
        begin
          MsgBox('Warning: an old version of {#AppName} is installed! Now the old one will be removed and the new one installed!', mbInformation, MB_OK);
          Exec(RemoveQuotes(Uninstall), ' /SILENT', '', SW_SHOWNORMAL, ewWaitUntilTerminated, ResultCode);
        end;
    end;
end;
