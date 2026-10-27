# word_caption.ps1 - put a caption paragraph immediately after a table, in the open document.
#
#   powershell.exe -File word_caption.ps1 -Table 1 -Text "..." [-NoTrack] [-Verify]
#
# Split out of word_table.ps1 because `$tbl.Range.Next(4, 1)` does NOT work from PowerShell:
# COM returns a Range, PowerShell binds the arguments as if Next were an indexed property and
# tries to coerce the result to Int32, so it fails with "Cannot convert System.__ComObject to
# type System.Int32" and takes the next six statements down with it. The table itself had
# already been built by then, which is the worst shape of failure - it looks like the whole
# insertion worked.
#
# The reliable idiom is to walk the paragraph collection and take the first paragraph that
# starts at or after the table's end and is NOT itself inside a table. Table cells each hold
# their own paragraphs, so `Information(12)` (wdWithInTable) is the only safe test.
#
# -Verify prints the table's first rows instead of writing anything, so the content can be
# checked before a caption is committed to it.

param(
  [int]$Table = 1,
  [string]$Text = "",
  [string]$Doc = "Lenzo et al 2027.docx",
  [switch]$NoTrack,
  [switch]$Verify
)

try { $w = [Runtime.InteropServices.Marshal]::GetActiveObject('Word.Application') }
catch { Write-Output "FAILED: Word is not running."; exit 1 }
$d = $null
foreach ($x in $w.Documents) { if ($x.Name -eq $Doc) { $d = $x } }
if ($null -eq $d) { Write-Output "FAILED: '$Doc' is not open."; exit 1 }
if ($d.Tables.Count -lt $Table) {
  Write-Output ("FAILED: document has " + $d.Tables.Count + " table(s); asked for " + $Table)
  exit 1
}
$tbl = $d.Tables.Item($Table)

if ($Verify) {
  Write-Output ("table " + $Table + ": " + $tbl.Rows.Count + " rows x " + $tbl.Columns.Count + " cols")
  for ($r = 1; $r -le [Math]::Min(8, $tbl.Rows.Count); $r++) {
    $line = @()
    for ($c = 1; $c -le $tbl.Columns.Count; $c++) {
      try { $line += ($tbl.Cell($r, $c).Range.Text -replace "`r|`a", "") } catch { $line += "<merged>" }
    }
    Write-Output ("  " + ($line -join " | "))
  }
  exit 0
}
if ($Text.Length -eq 0) { Write-Output "FAILED: -Text is empty."; exit 1 }

$tblEnd = $tbl.Range.End
$target = $null
foreach ($p in $d.Paragraphs) {
  if ($p.Range.Start -ge $tblEnd -and $p.Range.Information(12) -eq $false) { $target = $p; break }
}
if ($null -eq $target) { Write-Output "FAILED: no paragraph after the table."; exit 1 }
Write-Output ("  paragraph after table: '" + ($target.Range.Text -replace "`r|`a", "") + "'")

try { $target.Range.Select() } catch { }
if (-not $NoTrack) {
  try { $d.TrackRevisions = $true } catch { }
  $isOn = $false; try { $isOn = $d.TrackRevisions } catch { }
  if (-not $isOn) {
    Write-Output "FAILED: could not enable track changes (comment card holds focus). Nothing written."
    exit 1
  }
}

# If the paragraph after the table is already empty, write into it rather than adding another;
# a stray empty paragraph between a table and its caption is the usual result otherwise.
$existing = ($target.Range.Text -replace "`r|`a", "").Trim()
if ($existing.Length -eq 0) {
  $cp = $target
} else {
  $target.Range.InsertParagraphBefore()
  $cp = $d.Paragraphs.Item($target.Range.Paragraphs.Item(1).Range.Paragraphs.Count)
  # Re-resolve by position rather than trusting an index arithmetic that the insert just moved.
  foreach ($p in $d.Paragraphs) {
    if ($p.Range.Start -ge $tblEnd -and $p.Range.Information(12) -eq $false) { $cp = $p; break }
  }
}
$cp.Range.InsertBefore($Text)
$cp.Style = "Normal"
$cp.Range.Font.Size = 9
$cp.Range.Font.Italic = $true
Write-Output ("  caption written, " + $Text.Length + " chars; revisions " + $d.Revisions.Count)
Write-Output "  NOT saved. Review in Word, then accept or reject."
