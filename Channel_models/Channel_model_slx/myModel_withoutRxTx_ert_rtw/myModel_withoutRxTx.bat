set "VSCMD_START_DIR=%CD%"

call "setup_msvc150.bat"

call "%VS150COMNTOOLS%..\..\VC\Auxiliary\Build\VCVARSALL.BAT" AMD64

cd .

if "%1"=="" (nmake  -f myModel_withoutRxTx.mk all) else (nmake  -f myModel_withoutRxTx.mk %1)
@if errorlevel 1 goto error_exit

exit /B 0

:error_exit
echo The make command returned an error of %errorlevel%
An_error_occurred_during_the_call_to_make
