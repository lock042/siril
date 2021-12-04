@echo off
if [%1]==[] goto help
if [%2]==[] goto help
set VER=%~1
set SIRIL_BASE=%~2
set SIRIL64=%~3
set DEPS_BASE=%~4
set DEPS64=%~5

set INNOPATH=c:\program files (x86)\inno setup 6
echo "%INNOPATH%\iscc.exe"
if not exist "%INNOPATH%\iscc.exe" goto noinno

::i'd use %*, but shift has no effect on it
shift
shift
shift
shift
shift
shift
shift
set PARAMS=
:doparams
if "%1"=="" goto paramsdone
set PARAMS=%PARAMS% %1
shift
goto doparams
:paramsdone

"%INNOPATH%\iscc.exe" -DVERSION="%VER%" -DSIRIL_DIR="%SIRIL_BASE%" -DDIR64="%SIRIL64%" -DDEPS_DIR="%DEPS_BASE%" -DDDIR64="%DEPS64%" -DDEBUG_SYMBOLS -DPYTHON -DLUA %PARAMS% siril64.iss
goto :eof

:help
echo Usage: %~n0%~x0 ver.si.on siril_base_dir siril_x64_dir deps_base_dir deps_x64_dir [iscc_parameters]
echo Example: %~n0%~x0 2.9.4 X:\siril-output\2.9-dev x64 x:\siril-deps x64 -DPYTHON -DDEBUG_SYMBOLS
goto :eof
:noinno
echo Inno Setup path could not be read from Registry - install Inno Setup or set INNOPATH environment variable pointing at it's
echo install directory
goto :eof
