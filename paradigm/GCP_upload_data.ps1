# Upload GCP EEG files from the eego export folder to methlab_data.
#
# Reads C:\Users\neuro\Documents\Eego\export for files named <ID>_GCP*.cnt,
# <ID>_GCP*.evt, and <ID>_GCP*.seg. Each file is copied to
# \\idnas37.d.uzh.ch\g_psyplafor_methlab_data$\OCC\GCP\data\<ID>\
# which is V:\OCC\GCP\data\<ID> on the analysis computers.
#
# A filename counts as already uploaded when it is in that subject folder
# or in archive, with the same size. It is not copied again.
# Recordings already on the server are never overwritten or deleted.
# Unfinished copies from this script use the suffix .uploading and may be
# replaced on a later run. Z: is disconnected before this script ends.
#
# Run GCP_upload_data.bat. Keep this file in the same folder.

$ErrorActionPreference = 'Continue'
$localRoot = 'C:\Users\neuro\Documents\Eego\export'
$shareUnc = '\\idnas37.d.uzh.ch\g_psyplafor_methlab_data$'
$drive = 'Z:'
$serverRel = 'OCC\GCP\data'
$script:shouldDisconnect = $false
$script:stillMounted = $false

function Get-MappedUnc {
    param([Parameter(Mandatory = $true)][string]$DeviceId)
    $disk = Get-CimInstance -ClassName Win32_LogicalDisk -Filter ("DeviceID='{0}'" -f $DeviceId) -ErrorAction SilentlyContinue
    if ($null -eq $disk) {
        return $null
    }
    return [string]$disk.ProviderName
}

function Test-MethlabDrive {
    $unc = Get-MappedUnc -DeviceId $drive
    return ($unc -like '*g_psyplafor_methlab_data*')
}

function Disconnect-FailedAttempt {
    & net.exe use $drive /delete /y *> $null
}

function ConvertTo-PlainText {
    param([Parameter(Mandatory = $true)][System.Security.SecureString]$Secure)
    $bstr = [Runtime.InteropServices.Marshal]::SecureStringToBSTR($Secure)
    try {
        return [Runtime.InteropServices.Marshal]::PtrToStringBSTR($bstr)
    } finally {
        [Runtime.InteropServices.Marshal]::ZeroFreeBSTR($bstr)
    }
}

function Format-WinArg {
    param([Parameter(Mandatory = $true)][AllowEmptyString()][string]$Value)
    $sb = New-Object System.Text.StringBuilder
    [void]$sb.Append('"')
    $slashes = 0
    foreach ($ch in $Value.ToCharArray()) {
        if ($ch -eq '\') {
            $slashes++
        } elseif ($ch -eq '"') {
            [void]$sb.Append('\', ($slashes * 2 + 1))
            [void]$sb.Append('"')
            $slashes = 0
        } else {
            if ($slashes -gt 0) {
                [void]$sb.Append('\', $slashes)
                $slashes = 0
            }
            [void]$sb.Append($ch)
        }
    }
    if ($slashes -gt 0) {
        [void]$sb.Append('\', ($slashes * 2))
    }
    [void]$sb.Append('"')
    return $sb.ToString()
}

function Connect-MethlabData {
    $unc = Get-MappedUnc -DeviceId $drive
    if ($unc -like '*g_psyplafor_methlab_data*') {
        if (Test-Path -LiteralPath ($drive + '\') -PathType Container) {
            Write-Host 'methlab_data is already connected on Z:. That drive letter will be disconnected at the end.'
        } else {
            Write-Host 'Z: points at methlab_data but cannot be read. Connecting again.'
            & net.exe use $drive /delete /y *> $null
            if (Test-MethlabDrive) {
                throw 'Z: could not be disconnected. Nothing was copied.'
            }
            $unc = $null
        }
    }

    if ($unc -notlike '*g_psyplafor_methlab_data*') {
        if ($null -eq $unc) {
            # Z: is free.
        } elseif ([string]::IsNullOrWhiteSpace($unc)) {
            throw 'Z: is a local disk. Nothing was copied. This script left that disk alone.'
        } elseif ($unc -like '*methlab_xfr*') {
            Write-Host 'Z: is connected to the transfer share. That drive letter will be disconnected. Files on the server stay where they are.'
            & net.exe use $drive /delete /y *> $null
            if ($null -ne (Get-MappedUnc -DeviceId $drive)) {
                throw 'Z: could not be disconnected from the transfer share. Nothing was copied.'
            }
        } else {
            throw "Z: is already connected to $unc. Nothing was copied. That connection was left in place."
        }

        $mapped = $false
        $attempt = 0
        while ($attempt -lt 3) {
            $short = Read-Host 'UZH shortname'
            if ([string]::IsNullOrWhiteSpace($short)) {
                Write-Host 'Connection cancelled. Nothing was copied.'
                return $null
            }
            $short = $short.Trim()
            if ($short -match '^(?i)UZH[\\/](.+)$') {
                $short = $Matches[1].Trim()
            }
            if ($short -notmatch '^[A-Za-z0-9]+$') {
                Write-Host 'The UZH shortname can only contain letters and digits. Type the shortname only, without UZH\ in front.'
                continue
            }

            $secure = Read-Host 'UZH password' -AsSecureString
            if ($null -eq $secure -or $secure.Length -eq 0) {
                Write-Host 'The UZH password was empty. Try again.'
                continue
            }

            $attempt++
            $plain = $null
            $argLine = $null
            try {
                $plain = ConvertTo-PlainText -Secure $secure
                $argLine = @(
                    'use',
                    (Format-WinArg $drive),
                    (Format-WinArg $shareUnc),
                    (Format-WinArg $plain),
                    (Format-WinArg ('/user:UZH\' + $short)),
                    '/persistent:no'
                ) -join ' '
                $psi = New-Object System.Diagnostics.ProcessStartInfo
                $psi.FileName = Join-Path $env:SystemRoot 'System32\net.exe'
                $psi.Arguments = $argLine
                $psi.UseShellExecute = $false
                $proc = [System.Diagnostics.Process]::Start($psi)
                $proc.WaitForExit()
            } finally {
                if (-not [string]::IsNullOrEmpty($plain)) {
                    $plain = ' ' * $plain.Length
                }
                $plain = $null
                $argLine = $null
            }
            if (Test-MethlabDrive) {
                $mapped = $true
                break
            }
            Disconnect-FailedAttempt
            if ($attempt -lt 3) {
                Write-Host 'Connection failed. Check the UZH shortname and the UZH password.'
            }
        }
        if (-not $mapped) {
            throw 'methlab_data was not connected. Nothing was copied. Run this script again.'
        }
    }

    $script:shouldDisconnect = $true
    $serverRoot = Join-Path ($drive + '\') $serverRel
    if (-not (Test-Path -LiteralPath $serverRoot -PathType Container)) {
        throw ('Server GCP folder not found: {0}. This is OCC\GCP\data on g_psyplafor_methlab_data$. Nothing was copied. The folder was not created.' -f $serverRoot)
    }
    return $serverRoot
}

function Disconnect-MethlabData {
    if ($script:shouldDisconnect) {
        Write-Host ''
        Write-Host 'Disconnecting Z:.'
        & net.exe use $drive /delete /y *> $null
    }
    Write-Host ''
    if (Test-MethlabDrive) {
        $script:stillMounted = $true
        $unc = Get-MappedUnc -DeviceId $drive
        Write-Host 'METHLAB_DATA IS STILL MOUNTED'
        if (-not [string]::IsNullOrWhiteSpace($unc)) {
            Write-Host $unc
        }
    } else {
        $script:stillMounted = $false
        Write-Host 'METHLAB_DATA IS UNMOUNTED'
    }
    Write-Host ''
}

function Get-DirectoryNames {
    param([Parameter(Mandatory = $true)][string]$Root)
    return @(Get-ChildItem -LiteralPath $Root -Directory -ErrorAction Stop | ForEach-Object { $_.Name })
}

function Get-ServerSizes {
    param([Parameter(Mandatory = $true)][string]$SubjectDir)
    $map = @{}
    $folders = @(
        $SubjectDir
        (Join-Path $SubjectDir 'archive')
    )
    foreach ($folder in $folders) {
        if (-not (Test-Path -LiteralPath $folder -PathType Container)) {
            continue
        }
        $items = @(Get-ChildItem -LiteralPath $folder -File -ErrorAction Stop)
        foreach ($item in $items) {
            $key = $item.Name.ToLowerInvariant()
            if (-not $map.ContainsKey($key)) {
                $map[$key] = New-Object System.Collections.ArrayList
            }
            [void]$map[$key].Add([int64]$item.Length)
        }
    }
    Write-Output $map -NoEnumerate
}

function Test-SameSizeOnServer {
    param($SizeMap, [string]$Name, [int64]$Length)
    $key = $Name.ToLowerInvariant()
    if (-not $SizeMap.ContainsKey($key)) {
        return $false
    }
    foreach ($size in @($SizeMap[$key])) {
        if ([int64]$size -eq $Length) {
            return $true
        }
    }
    return $false
}

function Test-AnyCopyOnServer {
    param($SizeMap, [string]$Name)
    return $SizeMap.ContainsKey($Name.ToLowerInvariant())
}

function Get-SubjectId {
    param([Parameter(Mandatory = $true)][string]$Name)
    if ($Name -match '^(?<id>\d+)_GCP.*\.(cnt|evt|seg)$') {
        return [string]$Matches['id']
    }
    return ''
}

function Test-UnderRoot {
    param([string]$Root, [string]$Path)
    $rootFull = [System.IO.Path]::GetFullPath(($Root.TrimEnd('\') + '\'))
    if (-not $rootFull.EndsWith('\')) {
        $rootFull = $rootFull + '\'
    }
    $pathFull = [System.IO.Path]::GetFullPath($Path)
    # Windows PowerShell 5.1 rejects StartsWith(string, StringComparison).
    return $pathFull.ToLowerInvariant().StartsWith($rootFull.ToLowerInvariant())
}

function Copy-ToTemp {
    param(
        [Parameter(Mandatory = $true)][string]$Source,
        [Parameter(Mandatory = $true)][string]$TempPath,
        [Parameter(Mandatory = $true)][string]$DisplayName
    )
    $sourceStream = [System.IO.File]::Open(
        $Source,
        [System.IO.FileMode]::Open,
        [System.IO.FileAccess]::Read,
        [System.IO.FileShare]::ReadWrite)
    $destStream = $null
    try {
        # Create replaces only a previous .uploading file from this script.
        $destStream = [System.IO.File]::Open(
            $TempPath,
            [System.IO.FileMode]::Create,
            [System.IO.FileAccess]::Write,
            [System.IO.FileShare]::None)
        $buffer = New-Object byte[] 8388608
        $total = $sourceStream.Length
        $done = [int64]0
        $lastPct = -1
        while (($read = $sourceStream.Read($buffer, 0, $buffer.Length)) -gt 0) {
            $destStream.Write($buffer, 0, $read)
            $done += $read
            if ($total -gt 0) {
                $pct = [int](($done * 100) / $total)
                if ($pct -ne $lastPct) {
                    Write-Host -NoNewline ("`r    {0}  {1,3}%" -f $DisplayName, $pct)
                    $lastPct = $pct
                }
            }
        }
        $destStream.Flush()
    } finally {
        Write-Host ''
        if ($null -ne $destStream) {
            $destStream.Dispose()
        }
        $sourceStream.Dispose()
    }
}

function Show-NameList {
    param([string[]]$Lines)
    if ($Lines.Count -eq 0) {
        Write-Host '    (none)'
        return
    }
    foreach ($line in $Lines) {
        Write-Host ("    {0}" -f $line)
    }
}

function Invoke-EegUpload {
    if ($env:OS -ne 'Windows_NT') {
        throw 'This script is for the Windows EEG laptop.'
    }

    Write-Host 'GCP EEG upload to methlab_data'
    Write-Host ("Local folder:  {0}" -f $localRoot)
    Write-Host ("Server folder: {0}\{1}\<subject ID>" -f $shareUnc, $serverRel)
    Write-Host ''
    Write-Host 'Connect this laptop to the UZH network. The EEG amplifier network is not enough.'
    Read-Host 'Press Enter to continue' | Out-Null

    if (-not (Test-Path -LiteralPath $localRoot -PathType Container)) {
        throw "Local EEG folder not found: $localRoot. Nothing was copied."
    }

    $serverRoot = Connect-MethlabData
    if ([string]::IsNullOrWhiteSpace($serverRoot)) {
        return 0
    }
    Write-Host ("Connected. Checking {0}" -f $serverRoot)
    Write-Host ''

    $found = @()
    foreach ($filter in @('*.cnt', '*.evt', '*.seg')) {
        $found += @(Get-ChildItem -LiteralPath $localRoot -Recurse -File -Filter $filter -ErrorAction SilentlyContinue)
    }

    $matched = @()
    $unrecognized = @()
    foreach ($file in $found) {
        if ($null -eq $file -or $file.Name -notmatch 'GCP') {
            continue
        }
        if (-not [string]::IsNullOrWhiteSpace((Get-SubjectId -Name $file.Name))) {
            $matched += $file
        } else {
            $unrecognized += $file.FullName
        }
    }

    $stable = @()
    $busy = @()
    $empty = @()
    if ($matched.Count -gt 0) {
        $sizesBefore = @{}
        foreach ($file in $matched) {
            $sizesBefore[$file.FullName] = [int64]$file.Length
        }
        Start-Sleep -Seconds 2
        foreach ($file in $matched) {
            $file.Refresh()
            if ([int64]$file.Length -eq 0) {
                $empty += $file.FullName
            } elseif ([int64]$file.Length -ne $sizesBefore[$file.FullName]) {
                $busy += $file.FullName
            } else {
                $stable += $file
            }
        }
    }

    $byName = @{}
    foreach ($file in $stable) {
        $key = $file.Name.ToLowerInvariant()
        if (-not $byName.ContainsKey($key)) {
            $byName[$key] = New-Object System.Collections.Generic.List[object]
        }
        $byName[$key].Add($file)
    }
    $duplicates = @()
    $uniqueFiles = @()
    foreach ($key in $byName.Keys) {
        if ($byName[$key].Count -gt 1) {
            $duplicates += $byName[$key][0].Name
        } else {
            $uniqueFiles += $byName[$key][0]
        }
    }

    $toCopy = @()
    $already = @()
    $sizeMismatch = @()
    $byId = @{}
    foreach ($file in $uniqueFiles) {
        $id = Get-SubjectId -Name $file.Name
        if (-not $byId.ContainsKey($id)) {
            $byId[$id] = New-Object System.Collections.Generic.List[object]
        }
        $byId[$id].Add($file)
    }
    foreach ($id in @($byId.Keys)) {
        $subjectDir = Join-Path $serverRoot $id
        if (-not (Test-UnderRoot -Root $serverRoot -Path $subjectDir)) {
            throw "Refusing to use a subject path outside OCC\GCP\data: $id"
        }
        $sizeMap = Get-ServerSizes -SubjectDir $subjectDir
        foreach ($file in @($byId[$id])) {
            if (Test-SameSizeOnServer -SizeMap $sizeMap -Name $file.Name -Length ([int64]$file.Length)) {
                $already += $file.Name
            } elseif (Test-AnyCopyOnServer -SizeMap $sizeMap -Name $file.Name) {
                $sizeMismatch += ("{0}\{1}" -f $id, $file.Name)
            } else {
                $toCopy += [pscustomobject]@{
                    Id     = [string]$id
                    Name   = [string]$file.Name
                    Source = [string]$file.FullName
                    Length = [int64]$file.Length
                }
            }
        }
    }

    Write-Host ("Already on methlab_data: {0}" -f $already.Count)
    if ($unrecognized.Count -gt 0) {
        Write-Host ''
        Write-Host 'These files contain GCP but are not named <ID>_GCP. They were left in the export folder.'
        Show-NameList -Lines $unrecognized
    }
    if ($empty.Count -gt 0) {
        Write-Host ''
        Write-Host 'These files are empty and were left in the export folder.'
        Show-NameList -Lines $empty
    }
    if ($busy.Count -gt 0) {
        Write-Host ''
        Write-Host 'These files were still changing, so the export may still be running. Run this script again when the export has finished.'
        Show-NameList -Lines $busy
    }
    if ($duplicates.Count -gt 0) {
        Write-Host ''
        Write-Host 'These filenames occur more than once in the export folder. Those copies were left there.'
        Show-NameList -Lines $duplicates
    }
    if ($sizeMismatch.Count -gt 0) {
        Write-Host ''
        Write-Host 'These filenames are already in the subject folder or in archive, but the size differs. They were left unchanged.'
        Write-Host 'If a previous copy was interrupted, remove that server file by hand and run this script again.'
        Show-NameList -Lines $sizeMismatch
    }

    if ($toCopy.Count -eq 0) {
        Write-Host ''
        Write-Host 'Nothing new to upload.'
        $problemCount = $sizeMismatch.Count + $busy.Count + $duplicates.Count + $empty.Count + $unrecognized.Count
        if ($problemCount -eq 0) {
            Write-Host 'Every GCP .cnt, .evt, and .seg file named <ID>_GCP is already in the subject folder or in archive.'
        }
        return 0
    }

    $subjectIds = @($toCopy | ForEach-Object { $_.Id } | Select-Object -Unique | Sort-Object { [int64]$_ })
    Write-Host ''
    Write-Host ("These {0} files from {1} subjects are not on methlab_data yet." -f $toCopy.Count, $subjectIds.Count)
    Write-Host 'A file already in the subject folder or in archive, with the same size, is left as it is.'
    Write-Host 'Nothing already stored as .cnt, .evt, or .seg is overwritten or deleted.'
    Write-Host ''
    foreach ($id in $subjectIds) {
        Write-Host $id
        $names = @($toCopy | Where-Object { $_.Id -eq $id } | ForEach-Object { $_.Name } | Sort-Object)
        foreach ($name in $names) {
            Write-Host ("    {0}" -f $name)
        }
    }
    Write-Host ''
    $answer = Read-Host 'Type yes to upload, or press Enter to cancel'
    if ($answer.Trim().ToLowerInvariant() -ne 'yes') {
        Write-Host 'Upload cancelled. Nothing was copied.'
        return 0
    }

    $before = @(Get-DirectoryNames -Root $serverRoot)
    $uploaded = @()
    $failed = @()

    foreach ($id in $subjectIds) {
        $subjectDir = Join-Path $serverRoot $id
        if (-not (Test-UnderRoot -Root $serverRoot -Path $subjectDir)) {
            $failed += ("{0}: the subject path is outside OCC\GCP\data." -f $id)
            continue
        }
        if (Test-Path -LiteralPath $subjectDir -PathType Leaf) {
            $failed += ("{0}: a file already occupies this subject path. It was left unchanged." -f $id)
            continue
        }
        if (-not (Test-Path -LiteralPath $subjectDir -PathType Container)) {
            try {
                New-Item -ItemType Directory -Path $subjectDir -ErrorAction Stop | Out-Null
            } catch {
                $failed += ("{0}: the subject folder could not be created. {1}" -f $id, $_.Exception.Message)
                continue
            }
        }

        $filesForId = @($toCopy | Where-Object { $_.Id -eq $id })
        foreach ($file in $filesForId) {
            $sizeMap = Get-ServerSizes -SubjectDir $subjectDir
            if (Test-SameSizeOnServer -SizeMap $sizeMap -Name $file.Name -Length $file.Length) {
                continue
            }
            if (Test-AnyCopyOnServer -SizeMap $sizeMap -Name $file.Name) {
                $failed += ("{0}: {1} is already on methlab_data with a different size. It was left unchanged." -f $id, $file.Name)
                continue
            }

            $dest = Join-Path $subjectDir $file.Name
            $temp = $dest + '.uploading'
            if (-not (Test-UnderRoot -Root $serverRoot -Path $dest)) {
                $failed += ("{0}: {1} would land outside OCC\GCP\data." -f $id, $file.Name)
                continue
            }
            if (-not (Test-Path -LiteralPath $file.Source -PathType Leaf)) {
                $failed += ("{0}: {1} was not found in the export folder." -f $id, $file.Name)
                continue
            }

            Write-Host ("Copying {0} to {1}" -f $file.Name, $subjectDir)
            try {
                Copy-ToTemp -Source $file.Source -TempPath $temp -DisplayName $file.Name
                $tempItem = Get-Item -LiteralPath $temp -ErrorAction Stop
                if ([int64]$tempItem.Length -ne $file.Length) {
                    $failed += ("{0}: {1} was not fully written. The unfinished file is {2}. The .cnt, .evt, or .seg name was not created." -f $id, $file.Name, $tempItem.FullName)
                    continue
                }
                $sizeMap = Get-ServerSizes -SubjectDir $subjectDir
                if (Test-AnyCopyOnServer -SizeMap $sizeMap -Name $file.Name) {
                    $failed += ("{0}: {1} appeared on methlab_data during the copy. The unfinished file is {2}. The existing recording was left unchanged." -f $id, $file.Name, $temp)
                    continue
                }
                [System.IO.File]::Move($temp, $dest)
                $finalItem = Get-Item -LiteralPath $dest -ErrorAction Stop
                if ([int64]$finalItem.Length -eq $file.Length) {
                    $uploaded += ("{0}\{1}" -f $id, $file.Name)
                } else {
                    $failed += ("{0}: {1} is on methlab_data but the size does not match. It was left in place." -f $id, $file.Name)
                }
            } catch {
                $detail = $_.Exception.Message
                if (Test-Path -LiteralPath $temp -PathType Leaf) {
                    $detail = "$detail An unfinished file was left at $temp."
                }
                $failed += ("{0}: {1}: {2}" -f $id, $file.Name, $detail)
            }
        }
    }

    $after = @()
    $lost = @()
    try {
        $after = @(Get-DirectoryNames -Root $serverRoot)
        $afterSet = @{}
        foreach ($name in $after) {
            $afterSet[$name.ToLowerInvariant()] = $true
        }
        foreach ($name in $before) {
            if (-not $afterSet.ContainsKey($name.ToLowerInvariant())) {
                $lost += $name
            }
        }
    } catch {
        $failed += ('The subject folder list could not be read after the copy. {0}' -f $_.Exception.Message)
    }

    Write-Host ''
    Write-Host ("Uploaded ({0})" -f $uploaded.Count)
    Show-NameList -Lines $uploaded
    if ($failed.Count -gt 0) {
        Write-Host ''
        Write-Host ("Not uploaded ({0})" -f $failed.Count)
        Show-NameList -Lines $failed
    }
    if ($lost.Count -gt 0) {
        Write-Host ''
        Write-Host 'WARNING: subject folders that were already on the server are no longer listed.'
        Write-Host 'This script has no delete step. Check the share before recording again.'
        Show-NameList -Lines $lost
    }
    Write-Host ''
    Write-Host ("Subject folders now in OCC\GCP\data: {0}" -f $after.Count)

    if ($failed.Count -gt 0 -or $lost.Count -gt 0) {
        return 1
    }
    return 0
}

$code = 0
try {
    $code = Invoke-EegUpload
} catch {
    Write-Host ''
    Write-Host $_.Exception.Message
    $code = 1
} finally {
    Disconnect-MethlabData
}

if ($script:stillMounted -and $code -eq 0) {
    $code = 1
}

exit $code
