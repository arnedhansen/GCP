# Find which call throws "argument types do not match".
# Nothing is copied. The server is not mounted.
# Run GCP_upload_data_debug.bat and send back the SUMMARY block.

$ErrorActionPreference = 'Continue'
$localRoot = 'C:\Users\neuro\Documents\Eego\export'
$script:failNames = New-Object System.Collections.ArrayList

function Write-Probe {
    param(
        [Parameter(Mandatory = $true)][string]$Id,
        [Parameter(Mandatory = $true)][string]$Name,
        [Parameter(Mandatory = $true)][scriptblock]$Action
    )
    Write-Host ''
    Write-Host ("--- {0} {1}" -f $Id, $Name)
    try {
        & $Action
        Write-Host ("PASS {0}" -f $Id)
    } catch {
        $outer = $_.Exception.Message
        $inner = ''
        if ($null -ne $_.Exception.InnerException) {
            $inner = $_.Exception.InnerException.Message
        }
        Write-Host ("FAIL {0}" -f $Id)
        Write-Host ("OUTER {0}" -f $outer)
        if (-not [string]::IsNullOrWhiteSpace($inner)) {
            Write-Host ("INNER {0}" -f $inner)
        }
        Write-Host ("LINE {0}" -f $_.InvocationInfo.ScriptLineNumber)
        [void]$script:failNames.Add(("{0} {1}" -f $Id, $Name))
    }
}

Write-Host 'GCP upload debug. Nothing is copied. The server is not mounted.'
Write-Host ("PowerShell {0}" -f $PSVersionTable.PSVersion)
Write-Host ("CLR {0}" -f $PSVersionTable.CLRVersion)
Write-Host ("Local folder {0}" -f $localRoot)

$sample = $null
$gcpFile = $null

Write-Probe -Id '01' -Name 'export folder exists' -Action {
    if (-not (Test-Path -LiteralPath $localRoot -PathType Container)) {
        throw "Missing folder: $localRoot"
    }
    Write-Host 'folder exists'
}

Write-Probe -Id '02' -Name 'list cnt evt seg' -Action {
    $counts = @{}
    foreach ($filter in @('*.cnt', '*.evt', '*.seg')) {
        $items = @(Get-ChildItem -LiteralPath $localRoot -Recurse -File -Filter $filter -ErrorAction SilentlyContinue)
        $real = @($items | Where-Object { $null -ne $_ })
        $counts[$filter] = $real.Count
        if ($null -eq $script:sample -and $real.Count -gt 0) {
            $script:sample = $real[0]
        }
    }
    Write-Host ("cnt {0}  evt {1}  seg {2}" -f $counts['*.cnt'], $counts['*.evt'], $counts['*.seg'])
}

Write-Probe -Id '03' -Name 'GCP name and subject id' -Action {
    if ($null -eq $script:sample) {
        Write-Host 'no cnt evt or seg files, skipping name check'
        return
    }
    $name = $script:sample.Name
    Write-Host ("first file {0}" -f $name)
    if ($name -match '^(?<id>\d+)_GCP.*\.(cnt|evt|seg)$') {
        Write-Host ("subject id {0}" -f $Matches['id'])
        $script:gcpFile = $script:sample
    } else {
        Write-Host 'first file does not match ID_GCP'
    }
}

Write-Probe -Id '04' -Name 'file length as int64' -Action {
    $file = $script:gcpFile
    if ($null -eq $file) { $file = $script:sample }
    if ($null -eq $file) {
        Write-Host 'no sample file'
        return
    }
    $n = [int64]$file.Length
    Write-Host ("length {0}" -f $n)
}

Write-Probe -Id '05' -Name 'refresh file' -Action {
    $file = $script:gcpFile
    if ($null -eq $file) { $file = $script:sample }
    if ($null -eq $file) {
        Write-Host 'no sample file'
        return
    }
    $file.Refresh()
    Write-Host 'refresh returned'
}

Write-Probe -Id '06' -Name 'StartsWith with StringComparison' -Action {
    $ok = 'Z:\OCC\GCP\data\601\file.cnt'.StartsWith('Z:\OCC\GCP\data\', [System.StringComparison]::OrdinalIgnoreCase)
    Write-Host ("result {0}" -f $ok)
}

Write-Probe -Id '07' -Name 'StartsWith after ToLower' -Action {
    $path = 'Z:\OCC\GCP\data\601\file.cnt'.ToLowerInvariant()
    $root = 'Z:\OCC\GCP\data\'.ToLowerInvariant()
    $ok = $path.StartsWith($root)
    Write-Host ("result {0}" -f $ok)
}

Write-Probe -Id '08' -Name 'generic List object Add' -Action {
    $file = $script:gcpFile
    if ($null -eq $file) { $file = $script:sample }
    if ($null -eq $file) {
        $file = Get-Item -LiteralPath $localRoot
    }
    $list = New-Object System.Collections.Generic.List[object]
    $list.Add($file)
    Write-Host ("count {0}" -f $list.Count)
}

Write-Probe -Id '09' -Name 'ArrayList Add' -Action {
    $file = $script:gcpFile
    if ($null -eq $file) { $file = $script:sample }
    if ($null -eq $file) {
        $file = Get-Item -LiteralPath $localRoot
    }
    $list = New-Object System.Collections.ArrayList
    [void]$list.Add($file)
    Write-Host ("count {0}" -f $list.Count)
}

Write-Probe -Id '10' -Name 'generic List int64 Add' -Action {
    $list = New-Object System.Collections.Generic.List[int64]
    $list.Add([int64]1)
    Write-Host ("count {0}" -f $list.Count)
}

Write-Probe -Id '11' -Name 'group files the way the upload does' -Action {
    $byName = @{}
    $files = @()
    foreach ($filter in @('*.cnt', '*.evt', '*.seg')) {
        $files += @(Get-ChildItem -LiteralPath $localRoot -Recurse -File -Filter $filter -ErrorAction SilentlyContinue)
    }
    $n = 0
    foreach ($file in $files) {
        if ($null -eq $file) { continue }
        if ($file.Name -notmatch 'GCP') { continue }
        $key = $file.Name.ToLowerInvariant()
        if (-not $byName.ContainsKey($key)) {
            $byName[$key] = New-Object System.Collections.ArrayList
        }
        [void]$byName[$key].Add($file)
        $n++
        if ($n -ge 5) { break }
    }
    Write-Host ("grouped {0} files into {1} names" -f $n, $byName.Keys.Count)
}

Write-Host ''
Write-Host '======== SUMMARY ========'
Write-Host ("PowerShell {0}" -f $PSVersionTable.PSVersion)
Write-Host ("CLR {0}" -f $PSVersionTable.CLRVersion)
if ($script:failNames.Count -eq 0) {
    Write-Host 'Failed steps: none'
} else {
    Write-Host ("Failed steps: {0}" -f $script:failNames.Count)
    foreach ($name in $script:failNames) {
        Write-Host ("FAIL {0}" -f $name)
    }
}
Write-Host '======== END SUMMARY ========'
