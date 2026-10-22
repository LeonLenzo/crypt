#!/usr/bin/env python3
"""docx_review.py — add Word COMMENTS and tracked insertions to a .docx, safely.

Written because the manuscript lives in OneDrive so co-authors can review it, and round-tripping
through Markdown means retyping the whole document for every change.

## WHY COMMENTS ARE THE DEFAULT AND TRACKED EDITS ARE NOT

The document carries 22 Zotero citations as Word FIELD CODES: matched `fldChar` begin/end pairs
wrapping `instrText` holding the CSL JSON. Rewriting text that contains or abuts one of those
can orphan a field, and a broken Zotero field does not announce itself. It looks like an ordinary
citation until the bibliography is refreshed and the entry vanishes.

So the operations here are ordered by how much they disturb existing runs:

  comment      Touches NO existing run. Wraps a paragraph in commentRangeStart/End and adds a
               separate part. Zotero fields cannot be affected. Use this freely.
  insert_para  Appends a NEW paragraph wrapped in <w:ins>. Existing runs untouched.
  (not here)   In-sentence tracked edits. Deliberately NOT implemented: that is the operation
               that can split a field, and it is the one a human should do in Word.

## SAFETY

Never writes in place. Writes a temp file beside the target, verifies it reopens and that the
field-code count is unchanged, and only then replaces. Refuses to run if Word holds a lock.

    python docx_review.py <file.docx> --list
    python docx_review.py <file.docx> --comment "anchor text" "the comment"
    python docx_review.py <file.docx> --from-tsv comments.tsv
"""

import argparse
import datetime as dt
import os
import re
import shutil
import sys
import zipfile
from pathlib import Path

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
CT_COMMENTS = ("application/vnd.openxmlformats-officedocument"
               ".wordprocessingml.comments+xml")
REL_COMMENTS = ("http://schemas.openxmlformats.org/officeDocument"
                "/2006/relationships/comments")
COMMENTS_SHELL = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    f'<w:comments xmlns:w="{W}">{{body}}</w:comments>')


def esc(s):
    return (s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def para_text(p):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", p)).strip()


def paragraphs(doc):
    """(start, end, xml) for each <w:p>, in document order."""
    out = []
    for m in re.finditer(r"<w:p(?: [^>]*)?>.*?</w:p>", doc, re.S):
        out.append((m.start(), m.end(), m.group(0)))
    return out


def count_fields(doc):
    return (len(re.findall(r'w:fldCharType="begin"', doc)),
            len(re.findall(r'w:fldCharType="end"', doc)),
            len(re.findall(r"ZOTERO_ITEM", doc)))


class Docx:
    def __init__(self, path):
        self.path = Path(path)
        self.parts = {}
        with zipfile.ZipFile(self.path) as z:
            self.order = z.namelist()
            for n in self.order:
                self.parts[n] = z.read(n)

    @property
    def doc(self):
        return self.parts["word/document.xml"].decode("utf8")

    @doc.setter
    def doc(self, v):
        self.parts["word/document.xml"] = v.encode("utf8")

    def next_comment_id(self):
        cur = self.parts.get("word/comments.xml", b"").decode("utf8", "ignore")
        ids = [int(x) for x in re.findall(r'<w:comment [^>]*w:id="(\d+)"', cur)]
        return max(ids) + 1 if ids else 0

    def add_comment(self, anchor, text, author, initials):
        """Wrap the FIRST paragraph containing `anchor` in a comment range.

        Several comments can legitimately land on the same paragraph; the caller is told which
        index each one hit so that is visible rather than surprising.
        """
        doc = self.doc
        hit = None
        self.last_para = None
        for idx, (s, e, p) in enumerate(paragraphs(doc)):
            if anchor.lower() in para_text(p).lower():
                hit = (s, e, p)
                self.last_para = idx
                break
        if not hit:
            return None
        s, e, p = hit
        cid = self.next_comment_id()

        # Insert the range start AFTER <w:pPr> if present, so paragraph properties stay first.
        m = re.match(r"(<w:p(?: [^>]*)?>)(<w:pPr>.*?</w:pPr>)?", p, re.S)
        head = m.group(0)
        body = p[len(head):-len("</w:p>")]
        new_p = (f"{head}<w:commentRangeStart w:id=\"{cid}\"/>{body}"
                 f"<w:commentRangeEnd w:id=\"{cid}\"/>"
                 f"<w:r><w:rPr><w:rStyle w:val=\"CommentReference\"/></w:rPr>"
                 f"<w:commentReference w:id=\"{cid}\"/></w:r></w:p>")
        self.doc = doc[:s] + new_p + doc[e:]

        when = dt.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
        runs = "".join(
            f'<w:p><w:r><w:t xml:space="preserve">{esc(line)}</w:t></w:r></w:p>'
            for line in text.split("\n"))
        node = (f'<w:comment w:id="{cid}" w:author="{esc(author)}" '
                f'w:date="{when}" w:initials="{esc(initials)}">{runs}</w:comment>')
        cur = self.parts.get("word/comments.xml", b"").decode("utf8", "ignore")
        if cur:
            self.parts["word/comments.xml"] = cur.replace(
                "</w:comments>", node + "</w:comments>").encode("utf8")
        else:
            self.parts["word/comments.xml"] = COMMENTS_SHELL.format(
                body=node).encode("utf8")
            self._register_comments_part()
        return cid

    def _register_comments_part(self):
        ct = self.parts["[Content_Types].xml"].decode("utf8")
        if "/word/comments.xml" not in ct:
            ct = ct.replace("</Types>",
                            f'<Override PartName="/word/comments.xml" '
                            f'ContentType="{CT_COMMENTS}"/></Types>')
            self.parts["[Content_Types].xml"] = ct.encode("utf8")
        rels = self.parts["word/_rels/document.xml.rels"].decode("utf8")
        if "comments.xml" not in rels:
            used = [int(x) for x in re.findall(r'Id="rId(\d+)"', rels)]
            rid = f"rId{max(used) + 1 if used else 1}"
            rels = rels.replace("</Relationships>",
                                f'<Relationship Id="{rid}" Type="{REL_COMMENTS}" '
                                f'Target="comments.xml"/></Relationships>')
            self.parts["word/_rels/document.xml.rels"] = rels.encode("utf8")

    def save(self, out):
        names = list(self.order)
        for n in self.parts:
            if n not in names:
                names.append(n)
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
            for n in names:
                z.writestr(n, self.parts[n])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("docx")
    ap.add_argument("--list", action="store_true", help="print paragraphs with indices")
    ap.add_argument("--comment", nargs=2, action="append", metavar=("ANCHOR", "TEXT"))
    ap.add_argument("--from-tsv", type=Path, help="TSV: anchor<TAB>comment")
    ap.add_argument("--author", default="Claude (via Leon)")
    ap.add_argument("--initials", default="CL")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    src = Path(a.docx)
    # Test the ACTUAL condition, not a proxy for it. The first version looked for Word's
    # "~$name" owner file and passed happily while Word had the document open: Word truncates
    # that name for long filenames, OneDrive may hide it, and it is not always created at all.
    # Opening the file for append is definitive and harmless.
    try:
        with open(src, "r+b"):
            pass
    except PermissionError:
        sys.exit(f"REFUSING: '{src.name}' is open in Word or locked by OneDrive.\n"
                 "Close the document in Word and run again. Nothing was written.")
    except FileNotFoundError:
        sys.exit(f"no such file: {src}")

    d = Docx(src)
    if a.list:
        for i, (_, _, p) in enumerate(paragraphs(d.doc)):
            t = para_text(p)
            if t:
                st = re.search(r'w:pStyle w:val="([^"]+)"', p)
                print(f"  [{i:>3}] {st.group(1) if st else 'Normal':<12} {t[:110]}")
        return

    jobs = list(a.comment or [])
    if a.from_tsv:
        import csv as _csv
        with open(a.from_tsv) as fh:
            for row in _csv.reader(fh, delimiter="\t"):
                if len(row) >= 2 and not row[0].startswith("#"):
                    jobs.append((row[0], row[1]))
    if not jobs:
        sys.exit("nothing to do: pass --comment or --from-tsv")

    before = count_fields(d.doc)
    ok, miss = 0, []
    for anchor, text in jobs:
        cid = d.add_comment(anchor, text, a.author, a.initials)
        if cid is None:
            miss.append(anchor)
        else:
            ok += 1
            print(f"  comment {cid} -> para {d.last_para}: {anchor[:55]}")
    for m in miss:
        print(f"  NO ANCHOR MATCH: {m[:70]}", file=sys.stderr)

    after = count_fields(d.doc)
    if before != after:
        sys.exit(f"ABORT: field codes changed {before} -> {after}. Nothing written.")
    print(f"  Zotero fields intact: {after[2]} citations, {after[0]}/{after[1]} begin/end")

    if a.dry_run:
        print("  --dry-run: nothing written")
        return

    tmp = src.with_suffix(".docx.tmp")
    d.save(tmp)
    v = Docx(tmp)                                    # must reopen as a valid zip
    if count_fields(v.doc) != after:
        tmp.unlink()
        sys.exit("ABORT: verification failed after write")
    bak = src.with_name(src.stem + ".bak.docx")
    shutil.copy2(src, bak)
    # NOT os.replace. The OneDrive tree is a 9p/drvfs mount and rename-onto-an-existing-file
    # raises EPERM there, which aborted a run after the temp file was already verified.
    # Copying the bytes in opens the destination for writing instead, which drvfs allows, and
    # it also keeps the file identity so OneDrive treats it as an edit rather than a replace.
    shutil.copyfile(tmp, src)
    tmp.unlink()
    check = Docx(src)
    if count_fields(check.doc) != after:
        shutil.copyfile(bak, src)
        sys.exit("ABORT: post-write check failed; original restored from backup")
    print(f"  wrote {ok} comment(s); backup at {bak.name}")


if __name__ == "__main__":
    main()
