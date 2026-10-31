@echo off
setlocal
cd /d "%~dp0"
echo Starting FazPrep 3D Conformer Generator...
where streamlit >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    streamlit run app.py
) else (
    python -m streamlit run app.py
)
if %ERRORLEVEL% NEQ 0 (
    pause
)
