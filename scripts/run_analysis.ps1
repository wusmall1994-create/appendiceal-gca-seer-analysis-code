$ErrorActionPreference = 'Stop'

$projectDir = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$rscriptCommand = Get-Command Rscript -ErrorAction SilentlyContinue

if (-not $rscriptCommand) {
    throw 'Rscript was not found on PATH. Install R and add its bin directory to PATH.'
}

Push-Location $projectDir
try {
    & $rscriptCommand.Source '.\R\00_install_packages.R'
    if ($LASTEXITCODE -ne 0) { throw "Package installation failed with exit code $LASTEXITCODE" }

    & $rscriptCommand.Source '.\R\01_main_analysis.R'
    if ($LASTEXITCODE -ne 0) { throw "Main analysis failed with exit code $LASTEXITCODE" }
}
finally {
    Pop-Location
}

Write-Host 'Analysis completed. Generated files are under outputs/formal_next_phase/.'

