@echo off
setlocal
cd /d "%~dp0"
if not exist "%~dp0GCP_upload_data.ps1" (
    echo GCP_upload_data.ps1 is missing. Keep it in the same folder as this file.
    echo.
    pause
    exit /b 1
)
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0GCP_upload_data.ps1"
set "ERR=%ERRORLEVEL%"
echo.
pause
exit /b %ERR%
