# docxedit — editing the manuscript while it is open in Word

The manuscript is a `.docx` in OneDrive, kept OPEN in Word so it can be worked on while analysis
runs. That makes it a locked file, so these tools drive Word itself over COM rather than writing
the document.

    /mnt/c/Users/275080I/OneDrive - Curtin/Papers/Firsts/Review/Lenzo et al 2027.docx

| | |
|---|---|
| `word_review.ps1` | add comments (`-Tsv anchor<TAB>text`), `-List` dumps paragraphs |
| `word_reply.ps1` | reply to existing threads, attributed to Claude rather than the signed-in user |
| `docx_review.py` | file-level fallback for when Word is CLOSED. Works, but anchors only to whole paragraphs |

Run the PowerShell from WSL after copying to a Windows-visible path:

    cp tools/docxedit/word_review.ps1 /mnt/c/Users/275080I/
    powershell.exe -NoProfile -ExecutionPolicy Bypass \
      -File 'C:\Users\275080I\word_review.ps1' -Tsv 'C:\Users\275080I\comments.tsv'

Nothing is saved. Changes appear live in Word and the document is left dirty, so a human reviews
and saves.

## Comments, not in-sentence tracked edits

The document carries 22 Zotero citations as Word FIELD CODES. A broken field does not announce
itself: it looks like an ordinary citation until the bibliography is refreshed and the entry
vanishes. Comments cannot harm them, because Word maintains its own fields. Whole-paragraph
tracked insertions are also safe. In-sentence tracked edits are the operation that can split a
field, and are deliberately not implemented here.

## Five things that cost a cycle each

1. **`os.replace` fails on the OneDrive mount.** It is 9p/drvfs and rename-onto-existing raises
   EPERM. Copy the bytes in instead.
2. **Lock detection by `~$` owner file does not work.** Word truncates that name for long
   filenames and OneDrive may hide it, so the check passed while the document was open. Test
   the real condition: try to open the file `r+b`.
3. **Word stores CURLY quotes.** A TSV key containing a straight apostrophe never matches a
   comment body typed in Word. Both sides are normalised now.
4. **`Comments.Add` ignores any author argument.** Word stamps the signed-in user. Set
   `$w.UserName` before adding and restore it after.
5. **`printf '%s'` does not expand `\t`.** A TSV written that way has literal backslash-t and
   every row is silently skipped. Write it with Python.

## Counting tracked changes from the XML is wrong

`<w:ins` also matches `<w:instrText`, which is what Zotero field codes use. Counting that way
reported "29 tracked insertions" in a document that had none. Ask Word: `$d.Revisions.Count`.
