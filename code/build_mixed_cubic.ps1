$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
$buildDir = Join-Path $projectRoot 'tmp/pdfs/susskind-mixed-cubic-build'
New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
Push-Location (Join-Path $projectRoot 'manuscript/susskind-mixed-cubic')
try {
    foreach ($pass in 1..2) {
        $ErrorActionPreference = 'Continue'
        & pdflatex -interaction=nonstopmode -halt-on-error "-output-directory=$buildDir" main.tex *> (Join-Path $buildDir "pass$pass.txt")
        $ErrorActionPreference = 'Stop'
        if ($LASTEXITCODE -ne 0) { throw "LaTeX pass $pass failed" }
    }
    $log = Get-Content -Raw (Join-Path $buildDir 'main.log')
    if ($log -match 'undefined|Rerun to get|Overfull|Fatal error' -or $log -notmatch 'Output written on') { throw 'Inspect mixed-cubic build log.' }
    Copy-Item -LiteralPath (Join-Path $buildDir 'main.pdf') -Destination (Join-Path $projectRoot 'output/pdf/de-sitter-susskind-mixed-cubic.pdf')
    Copy-Item -LiteralPath (Join-Path $buildDir 'main.log') -Destination (Join-Path $projectRoot 'results/susskind-mixed-cubic/main.log')
} finally { Pop-Location }
