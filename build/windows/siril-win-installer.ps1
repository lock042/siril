# Fetching version number
$VersionPattern = "version : '(\d+).(\d+).(.+)'," 
Get-ChildItem meson.build |
    Select-String -Pattern $VersionPattern | select-object -First 1 |
    Foreach-Object {
        $MAJOR_VERSION, $MINOR_VERSION, $MICRO_VERSION = $_.Matches[0].Groups[1..3].Value
    }

$VERSIONSTR=$MAJOR_VERSION+'.'+$MINOR_VERSION+'.'+$MICRO_VERSION

# Preparing installer with Innoset
cd build\windows\installer
.\compile.bat $VERSIONSTR ..\..\.. siril-w64 ..\..\.. siril-w64

# Test if the installer was created and return success/failure
cd _Output
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

