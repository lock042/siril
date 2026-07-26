# Fetching version number
$VersionPattern = "version : '(\d+).(\d+).(.+)'," 
Get-ChildItem meson.build |
    Select-String -Pattern $VersionPattern | select-object -First 1 |
    Foreach-Object {
        $MAJOR_VERSION, $MINOR_VERSION, $MICRO_VERSION = $_.Matches[0].Groups[1..3].Value
    }

$VERSIONSTR=$MAJOR_VERSION+'.'+$MINOR_VERSION+'.'+$MICRO_VERSION
Write-Output ('Siril Version: '+$VERSIONSTR)

#storing cwd for later use
$RootDir=(Get-Item .).FullName
Write-Output ('Root Directory: '+$RootDir)

# Checking Inno-setup install
$a="path to executable: "
$INNOPATH=(iscc --shimgen-noop | Select-String $a) -split $a | select-object -Last 1 | ForEach-Object Trim
Write-Output ('Inno-setup executable: '+$INNOPATH)
#Check for existence of iscc.exe, otherwise, throw error message and exit
if (!$INNOPATH) {
    Write-Host "Inno-setup exe was not found on your system"
    Write-Host "Aborting"
    exit 1
}

# Running inno setup with parameters
cd build\windows\installer
$Output="..\..\..\WinInstaller"  #location to store installer

$Param1="-DVERSION="+$VERSIONSTR
$Param2="-DOUTPUT="+$Output
$Param3="-DROOTDIR="+$RootDir
$Param4="-DMSYSTEM="+$Env:MSYSTEM

&$INNOPATH $Param1 $Param2 $Param3 $Param4 siril64.iss

# Test if the installer was created and return success/failure
cd $Output
Write-Output ('Installer package directory: '+(Get-Item .).FullName)
$EXE_ROOT = 'siril-'+$VERSIONSTR+'-'+$Env:MSYSTEM.ToLower()+'-setup'
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


