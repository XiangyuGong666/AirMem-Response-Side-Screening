param([string]$Python = "python", [string]$OutputRoot = "$env:TEMP/AirMemResponseBuild")
$ErrorActionPreference = 'Stop'
New-Item -ItemType Directory -Force -Path $OutputRoot | Out-Null
$root = Split-Path -Parent $PSScriptRoot
$package = Join-Path $root 'src/airmem_response_screening'
& $Python -m PyInstaller --noconfirm --onefile --windowed `
  --name AirMem-Response-Side-Screening --paths (Join-Path $root src) `
  --exclude-module PyQt5 --exclude-module PyQt6 --exclude-module PySide2 --exclude-module PySide6 `
  --exclude-module IPython --exclude-module pytest --exclude-module torch --exclude-module tensorflow --exclude-module notebook --exclude-module jupyter `
  --add-data "$(Join-Path $package *.py):airmem_response_screening" `
  --add-data "$(Join-Path $package config):airmem_response_screening/config" `
  --add-data "$(Join-Path $package model_inputs):airmem_response_screening/model_inputs" `
  --collect-data matplotlib --hidden-import matplotlib.backends.backend_tkagg --hidden-import matplotlib.backends.backend_pdf `
  --distpath (Join-Path $OutputRoot dist) --workpath (Join-Path $OutputRoot build) `
  --specpath $OutputRoot (Join-Path $PSScriptRoot windows_entry.py)
if ($LASTEXITCODE) { exit $LASTEXITCODE }
