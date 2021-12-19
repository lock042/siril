# Fetching version number
$VersionPattern = "version : '(\d+).(\d+).(.+)'," 
Get-ChildItem meson.build |
    Select-String -Pattern $VersionPattern | select-object -First 1 |
    Foreach-Object {
        $MAJOR_VERSION, $MINOR_VERSION, $MICRO_VERSION = $_.Matches[0].Groups[1..3].Value
    }

$VERSIONSTR=$MAJOR_VERSION+'.'+$MINOR_VERSION+'.'+$MICRO_VERSION
Write-Output $VERSIONSTR

#storing cwd for later use
$RootDir=(Get-Item .).FullName
Write-Output $RootDir

# Preparing installer with Innoset
cd build\windows\installer
$INNOPATH="c:\program files (x86)\inno setup 6\iscc.exe"
Write-Output $INNOPATH
#Check for existence of icss.exe, otherwise, throw erroe message and exit
if (-not(Test-Path -Path $INNOPATH -PathType Leaf)) {
    Write-Host "Innosetup exe was not found at $INNOPATH"
    Write-Host "Aborting"
    exit 1
}

# Running inno setup with parameters
$Output="..\..\..\WinInstaller"  #location to store installer

$Param1="-DVERSION="+$VERSIONSTR
$Param2="-DOUTPUT="+$Output
$Param3="-DROOTDIR="+$RootDir

&$INNOPATH $Param1 $Param2 $Param3 siril64.iss

# Test if the installer was created and return success/failure
cd $Output
Write-Output (Get-Item .).FullName
$EXE_ROOT = 'siril-'+$VERSIONSTR+'-setup'
$EXE_NAME = $EXE_ROOT+'.exe'
$SHA256_NAME = $EXE_ROOT+'.SHA256SUMS'
$SHA512_NAME = $EXE_ROOT+'.SHA512SUMS'
If (Test-Path -Path $EXE_NAME ) {
    Get-FileHash -Path $EXE_NAME -Algorithm SHA256 | Out-File $SHA256_NAME
    Get-FileHash -Path $EXE_NAME -Algorithm SHA512 | Out-File $SHA512_NAME
    exit 0
} else {
    exit 1
}

