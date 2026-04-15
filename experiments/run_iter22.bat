@echo off
cd /d C:\Users\mic23\prototype-embed\chameleon-research-loop
C:\Users\mic23\miniconda3\envs\chameleon\python.exe experiments\iter_22_runner.py > experiments\iter_22_output.txt 2>&1
echo Exit code: %ERRORLEVEL%
