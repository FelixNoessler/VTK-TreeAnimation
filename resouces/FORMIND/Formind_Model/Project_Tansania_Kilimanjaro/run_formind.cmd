@echo off
set var1=FORMIND is currently running
set var2=Please wait a minute ..
set var3=Simulation finished.
echo ################################
echo. 
echo %var1%.
echo %var2%.
echo.    
echo ################################

set olddir=%CD%

cd ..\bin

Formind-full.exe %olddir%\formind_parameters\KiLi_FLM3_PFT6.par 1> %olddir%\stout.txt 2> %olddir%\sterr.txt

echo.
echo %var3%.
echo.    
echo ################################

Pause