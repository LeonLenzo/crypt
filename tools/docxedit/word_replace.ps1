# word_replace.ps1 - replace the text of ONE paragraph in the open document, via COM.
#
#   powershell.exe -File word_replace.ps1 -Para <index> -Expect "<prefix>" -TextFile <path>
#                                         [-NoTrack] [-Keep]
#
# -Para is the 1-based Word paragraph index, which INCLUDES table-cell paragraphs, so take it
# from a loop over $d.Paragraphs rather than from word_review.ps1 -List (that dump skips empty
# paragraphs and counts from 0).
#
# -Expect IS REQUIRED and is the safety interlock: the paragraph's current text must start with
# it or nothing is written. Paragraph indices move whenever anything above them is inserted, and
# a replace that lands one paragraph off destroys real prose with no obvious trace. Every other
# script here fails loudly; this one has to as well, and an index alone cannot do that.
#
# ONLY FOR TEXT THIS TOOLCHAIN WROTE. The manuscript carries 22 Zotero citations as Word FIELD
# CODES inside existing runs. Replacing a paragraph's whole range deletes any field inside it,
# and a broken field does not announce itself: it looks fine until the bibliography refreshes
# and the entry vanishes. Generated captions have no fields, which is why this is safe for them
# and must not be pointed at authored prose. -Keep leaves the old paragraph in place as a
# tracked deletion target rather than replacing in place, if a reviewer wants to see both.

param(
  [int]$Para,
  [string]$Expect,
  [string]$TextFile,
  [string]$Doc = "Lenzo et al 2027.docx",
  [switch]$NoTrack,
  [switch]$Keep
)

if ($Para -lt 1) { Write-Output "FAILED: -Para must be 1 or greater."; exit 1 }
if ([string]::IsNullOrWhiteSpace($Expect)) {
  Write-Output "FAILED: -Expect is required. It is the interlock that stops a stale index"
  Write-Output "        from overwriting the wrong paragraph."
  exit 1
}

try { $w = [Runtime.InteropServices.Marshal]::GetActiveObject('Word.Application') }
catch { Write-Output "FAILED: Word is not running."; exit 1 }
$d = $null
foreach ($x in $w.Documents) { if ($x.Name -eq $Doc) { $d = $x } }
if ($null -eq $d) { Write-Output "FAILED: '$Doc' is not open."; exit 1 }

if ($Para -gt $d.Paragraphs.Count) {
  Write-Output ("FAILED: -Para " + $Para + " is out of range (" + $d.Paragraphs.Count + " paragraphs).")
  exit 1
}
$target = $d.Paragraphs.Item($Para)
$current = ($target.Range.Text -replace "`r|`a", "")

# Word stores curly quotes; normalise both sides or an apostrophe in -Expect never matches.
$norm = { param($s) $s -replace [char]0x2018, "'" -replace [char]0x2019, "'" `
                        -replace [char]0x201C, '"' -replace [char]0x201D, '"' }
if (-not ((& $norm $current).StartsWith((& $norm $Expect)))) {
  Write-Output "FAILED: paragraph $Para does not start with -Expect. Nothing written."
  Write-Output ("        expected: " + $Expect.Substring(0, [Math]::Min(60, $Expect.Length)))
  Write-Output ("        found   : " + $current.Substring(0, [Math]::Min(60, $current.Length)))
  exit 1
}

$new = (Get-Content -LiteralPath $TextFile -Raw -Encoding UTF8).Trim()
if ($new.Length -eq 0) { Write-Output "FAILED: new text is empty."; exit 1 }

try { $target.Range.Select() } catch { }
if (-not $NoTrack) {
  try { $d.TrackRevisions = $true } catch { }
  $isOn = $false; try { $isOn = $d.TrackRevisions } catch { }
  if (-not $isOn) {
    Write-Output "FAILED: could not enable track changes (comment card holds focus). Nothing written."
    exit 1
  }
}

# Replace the TEXT, not the paragraph: taking the whole range including its mark would delete
# the paragraph and shift every index below it, and the caption's style would go with it.
$r = $target.Range
$r.MoveEnd(1, -1) | Out-Null       # wdCharacter, drop the paragraph mark
$r.Text = $new
Write-Output ("  paragraph " + $Para + " replaced, " + $current.Length + " -> " + $new.Length + " chars")
Write-Output ("  revisions " + $d.Revisions.Count)
Write-Output "  NOT saved. Review in Word, then accept or reject."
