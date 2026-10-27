# word_insert.ps1 - insert a whole new paragraph into the document ALREADY OPEN in Word, via COM.
#
# Whole-paragraph insertions are the only text edits that are safe to make from here: Word's
# Zotero field codes live inside existing runs, so an in-sentence edit can orphan one, while a
# new paragraph of its own cannot touch them. Citations in the inserted text are PLAIN TEXT and
# must be replaced with Zotero fields by hand, or they will vanish on a bibliography refresh.
#
#   powershell.exe -File word_insert.ps1 -TextFile <path> -After <paragraph index> [-NoTrack]
#
# -After takes the index printed by word_review.ps1 -List. The text file may hold SEVERAL
# paragraphs separated by BLANK LINES; a paragraph starting "# " becomes a Heading 1 and one
# starting "## " a Heading 2, with the marker stripped. Inserting is append-only, so unlike word_replace.ps1 there is no
# range to expand and no deletion to get wrong.

param(
  [string]$TextFile,
  [int]$After,
  [string]$Doc = "Lenzo et al 2027.docx",
  [switch]$NoTrack,
  [string]$Style = "Normal"
)

try { $w = [Runtime.InteropServices.Marshal]::GetActiveObject('Word.Application') }
catch { Write-Output "FAILED: Word is not running."; exit 1 }

$d = $null
foreach ($x in $w.Documents) { if ($x.Name -eq $Doc) { $d = $x } }
if ($null -eq $d) { Write-Output "FAILED: '$Doc' is not open."; exit 1 }

$raw = (Get-Content -LiteralPath $TextFile -Raw -Encoding UTF8).TrimEnd("`r", "`n")
if ($raw.Length -eq 0) { Write-Output "FAILED: text file is empty."; exit 1 }
$paras = [regex]::Split($raw, "(\r?\n){2,}") | Where-Object { $_.Trim().Length -gt 0 -and $_ -notmatch "^(\r?\n)+$" }

# Word indexes paragraphs from 1; the -List dump counts from 0.
$idx = $After + 1
if ($idx -lt 1 -or $idx -gt $d.Paragraphs.Count) {
  Write-Output ("FAILED: -After " + $After + " is out of range (document has " + $d.Paragraphs.Count + " paragraphs).")
  exit 1
}

$anchor = $d.Paragraphs.Item($idx)
Write-Output ("  inserting after: " + $anchor.Range.Text.Substring(0, [Math]::Min(70, $anchor.Range.Text.Length)))

# TrackRevisions is unreadable and unsettable while a comment card holds focus in the review
# pane, which is the normal state when Leon is reading comments. Selecting a body range moves
# focus back into the document. If tracking still cannot be turned on, ABORT rather than insert
# untracked text: an untracked insertion is indistinguishable from the author's own writing,
# which is the failure mode that looks like nothing happened.
try { $anchor.Range.Select() } catch { }
$prevTrack = $false
if (-not $NoTrack) {
  try { $d.TrackRevisions = $true; $prevTrack = $false } catch { }
  $isOn = $false
  try { $isOn = $d.TrackRevisions } catch { }
  if (-not $isOn) {
    Write-Output "FAILED: could not enable track changes (a comment card probably holds focus in the review pane)."
    Write-Output "        Click into the document body in Word, then run this again. Nothing was inserted."
    exit 1
  }
}

$n = 0
foreach ($para in $paras) {
  $text = $para.Trim()
  $st = $Style
  if ($text.StartsWith("## ")) { $st = "Heading 2"; $text = $text.Substring(3).Trim() }
  elseif ($text.StartsWith("# ")) { $st = "Heading 1"; $text = $text.Substring(2).Trim() }
  $d.Paragraphs.Item($idx + $n).Range.InsertParagraphAfter()
  $new = $d.Paragraphs.Item($idx + $n + 1)
  $new.Range.InsertBefore($text)
  $new.Style = $st
  $n++
  Write-Output ("  + " + $st + ", " + $text.Length + " chars: " + $text.Substring(0, [Math]::Min(55, $text.Length)))
}

try { $d.TrackRevisions = $prevTrack } catch { }

Write-Output ("  inserted " + $n + " paragraph(s) after " + $After + "; revisions now " + $d.Revisions.Count)
Write-Output "  NOT saved. Review in Word, then accept or reject."
