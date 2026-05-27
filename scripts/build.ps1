<#
.SYNOPSIS
    Build the manuscript on Windows (PowerShell 5.1+).

.DESCRIPTION
    1. Optionally runs ``python scripts/build_figures.py``.
    2. Compiles ``manuscript/main.tex`` with pdflatex + bibtex + pdflatex x2;
       PDF is written to ``build/main.pdf``.

.PARAMETER Figures
    Force figure regeneration before compiling.

.PARAMETER SkipFigures
    Skip the figure-generation step.

.PARAMETER Clean
    After a successful build, remove intermediate files in ``build/`` (optionally keep PDF).

.PARAMETER KeepPdf
    With ``-Clean``, keep ``build/main.pdf``.

.PARAMETER Engine
    LaTeX engine (default: pdflatex).

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File scripts\build.ps1 -Figures
#>

[CmdletBinding()]
param(
    [switch] $Figures,
    [switch] $SkipFigures,
    [switch] $Clean,
    [switch] $KeepPdf,
    [string] $Engine = "pdflatex"
)

$ErrorActionPreference = "Stop"

$RepoRoot   = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$Manuscript = Join-Path $RepoRoot "manuscript"
$BuildDir   = Join-Path $RepoRoot "build"
$MainTex    = "main.tex"
$Jobname    = "main"

$PythonExe = $null
foreach ($name in @("python", "py")) {
    $cmd = Get-Command $name -ErrorAction SilentlyContinue
    if ($cmd) { $PythonExe = $cmd.Source; break }
}

function Write-Step([string] $msg) {
    Write-Host ""
    Write-Host "==> $msg" -ForegroundColor Cyan
}

function Invoke-Native([string] $exe, [string[]] $arguments) {
    Write-Host "    $exe $($arguments -join ' ')" -ForegroundColor DarkGray
    & $exe @arguments
    if ($LASTEXITCODE -ne 0) {
        throw "command '$exe' exited with code $LASTEXITCODE"
    }
}

# ---------------------------------------------------------------------------
# 1. Figures
# ---------------------------------------------------------------------------
if (-not $SkipFigures) {
    $needFigures = $Figures.IsPresent
    $figDir = Join-Path $Manuscript "figures"
    $pdfs = @(Get-ChildItem -Path $figDir -Filter "*.pdf" -ErrorAction SilentlyContinue)
    if (-not (Test-Path $figDir) -or ($pdfs.Count -eq 0)) {
        $needFigures = $true
    }
    if ($needFigures) {
        Write-Step "Generating figures via build_figures.py"
        if (-not $PythonExe) {
            throw "python (or py) is required to build figures. Install Python 3.10+ and try again."
        }
        Invoke-Native $PythonExe @((Join-Path $RepoRoot "scripts\build_figures.py"))
    } else {
        Write-Step "Figures already present (use -Figures to force a rebuild)"
    }
}

# ---------------------------------------------------------------------------
# 2. LaTeX (run from manuscript/ so \input and \bibliography paths resolve)
# ---------------------------------------------------------------------------
New-Item -ItemType Directory -Force -Path $BuildDir | Out-Null
$BuildDirAbs = (Resolve-Path $BuildDir).Path

$prevBib = $env:BIBINPUTS
$env:BIBINPUTS = "$Manuscript" + $(if ($env:BIBINPUTS) { ";$env:BIBINPUTS" } else { "" })

Push-Location $Manuscript
try {
    $texArgs = @(
        "-interaction=nonstopmode",
        "-halt-on-error",
        "-output-directory=$BuildDirAbs",
        $MainTex
    )

    Write-Step "Running $Engine (pass 1)"
    Invoke-Native $Engine $texArgs

    Write-Step "Running bibtex"
    Push-Location $BuildDirAbs
    try {
        Invoke-Native bibtex @($Jobname)
    } catch {
        Write-Warning "bibtex returned a non-zero status; check $BuildDirAbs\$Jobname.blg"
    } finally {
        Pop-Location
    }

    Write-Step "Running $Engine (pass 2)"
    Invoke-Native $Engine $texArgs

    Write-Step "Running $Engine (pass 3, final cross-references)"
    Invoke-Native $Engine $texArgs
} finally {
    Pop-Location
    if ($null -ne $prevBib) { $env:BIBINPUTS = $prevBib } else { Remove-Item Env:BIBINPUTS -ErrorAction SilentlyContinue }
}

$Pdf = Join-Path $BuildDir "$Jobname.pdf"
if (-not (Test-Path $Pdf)) {
    throw "LaTeX finished but $Pdf was not produced. Check $BuildDir\$Jobname.log"
}

Write-Host ""
Write-Host "PDF built successfully: $Pdf" -ForegroundColor Green

# ---------------------------------------------------------------------------
# 3. Cleanup (optional)
# ---------------------------------------------------------------------------
if ($Clean) {
    if ($KeepPdf) {
        Get-ChildItem -Path $BuildDir -Force | Where-Object { $_.Name -ne "$Jobname.pdf" } |
            Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
        Write-Host "Cleaned intermediate files in $BuildDir (kept $Jobname.pdf)" -ForegroundColor Yellow
    } else {
        Remove-Item -Recurse -Force $BuildDir
        Write-Host "Removed $BuildDir" -ForegroundColor Yellow
    }
}
