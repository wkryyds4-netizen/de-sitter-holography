param([string]$PdfLaTeX = 'pdflatex', [string]$BibTeX = 'bibtex')
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path $PSScriptRoot -Parent
$buildDir = Join-Path $projectRoot 'tmp/pdfs/latex-build'
$outputDir = Join-Path $projectRoot 'output/pdf'
New-Item -ItemType Directory -Force -Path $buildDir, $outputDir | Out-Null
Push-Location (Join-Path $projectRoot 'manuscript')
try {
    & $PdfLaTeX -interaction=nonstopmode -halt-on-error "-output-directory=$buildDir" main.tex *> (Join-Path $buildDir 'pass1.txt')
    if ($LASTEXITCODE -ne 0) { throw "First LaTeX pass failed; see $buildDir/pass1.txt" }
    # BibTeX resolves the bibliography from this source directory.
    & $BibTeX (Join-Path $buildDir 'main') *> (Join-Path $buildDir 'bibtex.txt')
    if ($LASTEXITCODE -ne 0) { throw "BibTeX failed; see $buildDir/bibtex.txt" }
    foreach ($pass in 2..3) {
        & $PdfLaTeX -interaction=nonstopmode -halt-on-error "-output-directory=$buildDir" main.tex *> (Join-Path $buildDir "pass$pass.txt")
        if ($LASTEXITCODE -ne 0) { throw "LaTeX pass $pass failed" }
    }
    $log = Get-Content -Raw (Join-Path $buildDir 'main.log')
    if ($log -match 'undefined|Rerun to get|Overfull|^!|Fatal error' -or $log -notmatch 'Output written on') { throw 'LaTeX did not produce a clean output; inspect main.log.' }
    if (-not (Test-Path -LiteralPath (Join-Path $buildDir 'main.bbl'))) { throw 'Missing compiled bibliography.' }
    Copy-Item -LiteralPath (Join-Path $buildDir 'main.pdf') -Destination (Join-Path $outputDir 'de-sitter-observer-foundations.pdf') -Force
    Write-Output (Join-Path $outputDir 'de-sitter-observer-foundations.pdf')
} finally { Pop-Location }
