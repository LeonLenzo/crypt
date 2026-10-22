# word_reply.ps1 - reply to existing comments in the open document.
#
# Sets Word's UserName for the duration so replies are attributed to Claude rather than to
# whoever is signed in. Word stamps the CURRENT user on every comment, which is why the first
# batch all came out as "Leon Lenzo"; the -Author argument on the add script was never applied
# and has been removed rather than left looking like it worked.
param([string]$Tsv, [string]$Doc = "Lenzo et al 2027.docx", [string]$As = "Claude")

$w = [Runtime.InteropServices.Marshal]::GetActiveObject('Word.Application')
$d = $null; foreach ($x in $w.Documents) { if ($x.Name -eq $Doc) { $d = $x } }
if ($null -eq $d) { Write-Output "FAILED: '$Doc' not open"; exit 1 }

$prevUser = $w.UserName
$prevInit = $w.UserInitials
$w.UserName = $As
$w.UserInitials = "CL"

$n = 0
foreach ($line in Get-Content -LiteralPath $Tsv -Encoding UTF8) {
  if ($line.Trim().Length -eq 0 -or $line.StartsWith("#")) { continue }
  $parts = $line -split "`t", 2
  if ($parts.Count -lt 2) { continue }
  $key = $parts[0]; $text = $parts[1]
  # Word stores curly quotes, so a key containing a straight ' or " never matches a comment
  # body typed in Word. Normalise both sides before comparing.
  $norm = { param($s) $s -replace [char]0x2018,"'" -replace [char]0x2019,"'" `
                        -replace [char]0x201C,'"' -replace [char]0x201D,'"' }
  $key = & $norm $key
  for ($i = 1; $i -le $d.Comments.Count; $i++) {
    $c = $d.Comments.Item($i)
    $body = & $norm (($c.Range.Text -replace "`r|`a"," "))
    if ($body -like "*$key*") {
      $null = $c.Replies.Add($c.Range, $text)
      Write-Output ("  replied to comment $i")
      $n++
      break
    }
  }
}
$w.UserName = $prevUser
$w.UserInitials = $prevInit
Write-Output ("  $n replies added; comments now " + $d.Comments.Count + "; NOT saved")
