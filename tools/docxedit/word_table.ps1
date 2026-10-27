# word_table.ps1 - insert a TABLE into the document ALREADY OPEN in Word, via COM.
#
# word_insert.ps1 inserts paragraphs; this inserts a real Word table, which a supplementary
# table has to be if it is going to live in the manuscript rather than beside it as a .tsv.
#
#   powershell.exe -File word_table.ps1 -Tsv <path> -After <paragraph index> [-NoTrack]
#                                       [-Caption <text>] [-HeaderRows 1]
#
# -After takes the index printed by word_review.ps1 -List. The TSV's first row is the header.
# A data row whose FIRST cell is empty and whose remaining cells are empty is skipped; a row
# whose first cell starts "## " becomes a bold sub-heading spanning the row, with the marker
# stripped, which is how the exclusion block is separated from the funnel block without two
# tables.
#
# WHY A TRACKED INSERTION AND WHY IT CAN ABORT
# TrackRevisions is unreadable and unsettable while a comment card holds focus in Word's review
# pane, which is the normal state when the document is being read. Selecting a body range moves
# focus back into the document. If tracking still cannot be turned on this ABORTS rather than
# inserting untracked, because an untracked insertion is indistinguishable from the author's own
# writing - the same failure mode word_insert.ps1 was fixed for, where it looks like nothing
# happened.
#
# ZOTERO IS SAFE HERE for the same reason whole-paragraph insertions are: the 22 citations are
# Word FIELD CODES inside existing runs, and a new table touches none of them. Any citation in
# the caption text is PLAIN TEXT and must be replaced with a Zotero field by hand.
#
# Nothing is saved. The table appears live and the document is left dirty for a human to accept.

param(
  [string]$Tsv,
  [int]$After,
  [string]$Doc = "Lenzo et al 2027.docx",
  [switch]$NoTrack,
  [string]$Caption = "",
  [int]$HeaderRows = 1
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

$lines = @(Get-Content -LiteralPath $Tsv -Encoding UTF8 | Where-Object { $_.Trim().Length -gt 0 })
if ($lines.Count -lt 2) { Write-Output "FAILED: need a header row and at least one data row."; exit 1 }
$rows = @()
foreach ($l in $lines) { $rows += ,($l -split "`t") }
$nCol = $rows[0].Count
$nRow = $rows.Count

# Word indexes paragraphs from 1; the -List dump counts from 0.
$idx = $After + 1
if ($idx -lt 1 -or $idx -gt $d.Paragraphs.Count) {
  Write-Output ("FAILED: -After " + $After + " is out of range (document has " + $d.Paragraphs.Count + " paragraphs).")
  exit 1
}
$anchor = $d.Paragraphs.Item($idx)
Write-Output ("  inserting after: " + $anchor.Range.Text.Substring(0, [Math]::Min(70, $anchor.Range.Text.Length)))
Write-Output ("  table: " + $nRow + " rows x " + $nCol + " cols")

try { $anchor.Range.Select() } catch { }
if (-not $NoTrack) {
  try { $d.TrackRevisions = $true } catch { }
  $isOn = $false
  try { $isOn = $d.TrackRevisions } catch { }
  if (-not $isOn) {
    Write-Output "FAILED: could not enable track changes (a comment card probably holds focus in the review pane)."
    Write-Output "        Click into the document body in Word, then run this again. Nothing was inserted."
    exit 1
  }
}

# Build the anchor paragraph for the table, then the caption paragraph after it, so the caption
# ends up BELOW the table without having to re-find anything once the table exists.
$anchor.Range.InsertParagraphAfter()
$tblPara = $d.Paragraphs.Item($idx + 1)
$tblPara.Style = "Normal"
$rng = $tblPara.Range

$tbl = $d.Tables.Add($rng, $nRow, $nCol)
$tbl.Range.Style = "Normal"
$tbl.Range.Font.Size = 9
$tbl.Rows.AllowBreakAcrossPages = $false
$tbl.PreferredWidthType = 2       # wdPreferredWidthPercent
$tbl.PreferredWidth = 100

# Hairlines only, and only horizontal: a grid of boxes reads as a spreadsheet, not a table.
foreach ($b in 1, 2, 3, 4, 5, 6) { $tbl.Borders.Item($b).LineStyle = 0 }
$tbl.Borders.Item(-1).LineStyle = 1   # wdBorderTop
$tbl.Borders.Item(-1).LineWidth = 8   # 1pt
$tbl.Borders.Item(-3).LineStyle = 1   # wdBorderBottom
$tbl.Borders.Item(-3).LineWidth = 8

for ($r = 0; $r -lt $nRow; $r++) {
  $cells = $rows[$r]
  $isSub = $cells[0].StartsWith("## ")
  if ($isSub) {
    $txt = $cells[0].Substring(3).Trim()
    $c = $tbl.Cell($r + 1, 1)
    $c.Range.Text = $txt
    $c.Range.Font.Bold = $true
    $c.Range.Font.Italic = $true
    # Merge the row so the sub-heading reads as a band rather than a stray first cell.
    # GUARDED: the table is itself a tracked insertion, and merging cells inside one is the
    # operation most likely to be refused. A failed merge leaves a readable sub-heading in the
    # first cell, so it must not abort an otherwise complete table.
    try { $tbl.Cell($r + 1, 1).Merge($tbl.Cell($r + 1, $nCol)) }
    catch { Write-Output ("  note: could not merge the sub-heading row " + ($r + 1) + "; left in cell 1") }
    continue
  }
  for ($c = 0; $c -lt $nCol; $c++) {
    $v = if ($c -lt $cells.Count) { $cells[$c] } else { "" }
    $cell = $tbl.Cell($r + 1, $c + 1)
    $cell.Range.Text = $v
    if ($c -gt 0) { $cell.Range.ParagraphFormat.Alignment = 2 }   # right-align the numbers
  }
}

# Header: bold, with a rule under it. Repeat on page breaks so a split table stays readable.
for ($r = 1; $r -le $HeaderRows; $r++) {
  $tbl.Rows.Item($r).Range.Font.Bold = $true
  $tbl.Rows.Item($r).HeadingFormat = $true
  $tbl.Rows.Item($r).Borders.Item(-3).LineStyle = 1
  $tbl.Rows.Item($r).Borders.Item(-3).LineWidth = 8
}

if ($Caption.Length -gt 0) {
  # The table consumed the paragraph it was built on, so the caption goes after the table's
  # own trailing paragraph rather than at a remembered index.
  $after = $tbl.Range.Next(4, 1)   # wdParagraph
  if ($null -ne $after) {
    $after.InsertParagraphBefore()
    $cp = $after.Paragraphs.Item(1)
    $cp.Range.InsertBefore($Caption)
    $cp.Style = "Normal"
    $cp.Range.Font.Size = 9
    $cp.Range.Font.Italic = $true
    Write-Output ("  caption: " + $Caption.Length + " chars")
  } else {
    Write-Output "  WARNING: could not place the caption; add it by hand."
  }
}

Write-Output ("  tables now " + $d.Tables.Count + "; revisions " + $d.Revisions.Count)
Write-Output "  NOT saved. Review in Word, then accept or reject."
