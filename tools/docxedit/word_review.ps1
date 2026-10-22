# word_review.ps1 - add comments / tracked edits to a document ALREADY OPEN in Word, via COM.
#
# WHY COM RATHER THAN EDITING THE .docx
# The manuscript lives in OneDrive and is kept open in Word so Leon can work on it while the
# analysis runs here. An open document is locked, so writing the .docx file directly fails with
# EPERM, and waiting for Word to close defeats the point. Driving Word itself sidesteps the lock
# entirely: Word owns the file, so Word does the writing, and the changes appear live.
#
# It is also MORE precise. Editing the XML, a comment can only be anchored to a whole paragraph
# without splitting runs; through Word, Range.Find anchors to the exact phrase.
#
# The Zotero field codes are safe here for the same reason: Word maintains its own fields, so
# adding a comment cannot orphan one the way a hand-edited run can.
#
#   powershell.exe -File word_review.ps1 -Tsv <path> [-Doc <name>] [-Track] [-List]
#
# TSV: anchor <TAB> comment     (lines starting # ignored)

param(
  [string]$Tsv,
  [string]$Doc = "Lenzo et al 2027.docx",
  [switch]$Track,
  [switch]$List,
  [string]$Author = "Claude"
)

try { $w = [Runtime.InteropServices.Marshal]::GetActiveObject('Word.Application') }
catch { Write-Output "FAILED: Word is not running. Open the document first."; exit 1 }

$d = $null
foreach ($x in $w.Documents) { if ($x.Name -eq $Doc) { $d = $x } }
if ($null -eq $d) {
  Write-Output "FAILED: '$Doc' is not open. Open documents:"
  foreach ($x in $w.Documents) { Write-Output ("  " + $x.Name) }
  exit 1
}

if ($List) {
  $i = 0
  foreach ($p in $d.Paragraphs) {
    $t = $p.Range.Text.Trim()
    if ($t.Length -gt 0) { Write-Output ("[{0,3}] {1}" -f $i, $t.Substring(0, [Math]::Min(110, $t.Length))) }
    $i++
  }
  exit 0
}

# Record and restore, so the user's own setting is not silently changed.
$prevTrack = $d.TrackRevisions
if ($Track) { $d.TrackRevisions = $true }

$before = $d.Comments.Count
$ok = 0; $miss = @()

foreach ($line in Get-Content -LiteralPath $Tsv) {
  if ($line.Trim().Length -eq 0 -or $line.StartsWith("#")) { continue }
  $parts = $line -split "`t", 2
  if ($parts.Count -lt 2) { continue }
  $anchor = $parts[0]; $text = $parts[1]

  # Fresh range over the whole document for every search, so one failed match cannot
  # silently shift where the next one lands.
  $r = $d.Content
  $f = $r.Find
  $f.ClearFormatting()
  $f.Text = $anchor
  $f.Forward = $true
  $f.Wrap = 0           # wdFindStop: do not wrap past the end
  $f.MatchCase = $false
  $f.MatchWildcards = $false

  if ($f.Execute()) {
    $null = $d.Comments.Add($r, $text)
    $ok++
    Write-Output ("  comment on: " + $anchor.Substring(0, [Math]::Min(60, $anchor.Length)))
  } else {
    $miss += $anchor
  }
}

foreach ($m in $miss) {
  Write-Output ("  NO MATCH: " + $m.Substring(0, [Math]::Min(70, $m.Length)))
}

if ($Track) { } else { $d.TrackRevisions = $prevTrack }

Write-Output ("  comments {0} -> {1}; revisions {2}; TrackRevisions={3}" -f `
  $before, $d.Comments.Count, $d.Revisions.Count, $d.TrackRevisions)
Write-Output "  NOT saved. Review in Word, then save yourself."
