#!/usr/bin/env python3
"""Remove hard line wrapping from prose, leaving code and structure intact.

Each prose paragraph becomes a single line containing no newline characters; blank lines continue to separate paragraphs. Code is untouched, and so is anything whose structure depends on line breaks: fenced code blocks, tables, lists, headings, block quotes, display mathematics and YAML front matter.

Comments and docstrings inside source files are prose and are unwrapped as well, since the line-length limit applies only to code.
"""

from __future__ import annotations

import re
import sys


# Markdown constructs whose meaning depends on the line break, and which are therefore copied through unchanged.
_STRUCTURAL = re.compile(
    r"^\s*("
    r"#{1,6}\s"          # heading
    r"|[-*+]\s"          # bullet list
    r"|\d+[.)]\s"        # ordered list
    r"|>"                # block quote
    r"|\|"               # table row
    r"|```"              # fence
    r"|:{3}"             # container
    r"|---\s*$"          # rule or front-matter delimiter
    r"|===\s*$"
    r"|\$\$"             # display mathematics
    r")"
)


def unwrap_markdown(text: str) -> str:
    """Join wrapped prose paragraphs in a Markdown document."""
    lines = text.split("\n")
    out: list[str] = []
    buffer: list[str] = []
    in_fence = False
    in_math = False

    def flush():
        if buffer:
            out.append(" ".join(s.strip() for s in buffer))
            buffer.clear()

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("```"):
            flush()
            in_fence = not in_fence
            out.append(line)
            continue
        if in_fence:
            out.append(line)
            continue

        # Display mathematics is copied verbatim: the line breaks inside it are meaningful to the renderer.
        if stripped == "$$":
            flush()
            in_math = not in_math
            out.append(line)
            continue
        if in_math:
            out.append(line)
            continue

        if not stripped:
            flush()
            out.append("")
            continue

        if _STRUCTURAL.match(line):
            flush()
            out.append(line)
            continue

        buffer.append(line)

    flush()
    return "\n".join(out)


def _unwrap_comment_block(block: list[str], marker: str, indent: str) -> list[str]:
    """Join a run of comment lines into paragraphs.

    A comment line carrying no text separates paragraphs, and a line that looks like a list item, a table row, a heading or code is left alone, since its break is meaningful.
    """
    result: list[str] = []
    buffer: list[str] = []
    prefix = ""

    def flush():
        if buffer:
            result.append(indent + marker + " " + prefix + " ".join(buffer))
            buffer.clear()

    for body in block:
        text = body.strip()
        if not text:
            flush()
            prefix = ""
            result.append(indent + marker.rstrip())
            continue
        # A table row, or an indented line, is code or layout and keeps its break.
        if text.startswith("|") or body.startswith("      "):
            flush()
            prefix = ""
            result.append(indent + marker + " " + body.rstrip())
            continue
        # A list marker opens a new item; its continuation lines join it, so that the whole item becomes one line.
        m = re.match(r"^([-*+]\s+|\d+[.)]\s+)(.*)$", text)
        if m:
            flush()
            prefix = m.group(1)
            buffer.append(m.group(2))
            continue
        buffer.append(text)

    flush()
    return result


def unwrap_fortran(text: str) -> str:
    """Unwrap runs of Fortran comment lines, leaving statements untouched."""
    lines = text.split("\n")
    out: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r"^(\s*)(!+>?!?|!!)\s?(.*)$", line)
        if m and line.strip().startswith("!"):
            indent = m.group(1)
            marker = re.match(r"^\s*(!+[>!]?)", line).group(1)
            # '!>' opens a documentation comment that continues with '!!'; the two are one block, written out with the opening marker.
            family = {"!>": "!!", "!!": "!!"}.get(marker, marker)
            block = []
            while i < len(lines):
                nxt = lines[i]
                m2 = re.match(r"^\s*(!+[>!]?)\s?(.*)$", nxt)
                if not (nxt.strip().startswith("!")):
                    break
                nxt_marker = re.match(r"^\s*(!+[>!]?)", nxt).group(1)
                if {"!>": "!!", "!!": "!!"}.get(nxt_marker, nxt_marker) != family:
                    break
                # A banner of repeated punctuation is a separator, not prose.
                body = m2.group(2)
                if re.fullmatch(r"[=\-*#]{4,}\s*", body):
                    break
                block.append(body)
                i += 1
            if block:
                out.extend(_unwrap_comment_block(block, marker, indent))
                continue
        out.append(line)
        i += 1
    return "\n".join(out)


def unwrap_python(text: str) -> str:
    """Unwrap docstrings and runs of hash comments in a Python file."""
    lines = text.split("\n")
    out: list[str] = []
    i = 0
    in_doc = False
    doc_quote = ""
    doc_indent = ""
    doc_body: list[str] = []

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if not in_doc:
            m = re.match(r'^(\s*)(?:[rRbBuUfF]*)("""|\'\'\')(.*)$', line)
            if m and not _closes_on_same_line(line, m.group(2)):
                in_doc = True
                doc_quote = m.group(2)
                doc_indent = m.group(1)
                out.append(line)
                doc_body = []
                i += 1
                continue

            if stripped.startswith("#"):
                indent = line[: len(line) - len(line.lstrip())]
                block = []
                while i < len(lines) and lines[i].strip().startswith("#"):
                    body = lines[i].strip()[1:]
                    if re.fullmatch(r"\s*[=\-*#]{4,}\s*", body):
                        break
                    block.append(body.strip())
                    i += 1
                if block:
                    out.extend(_unwrap_comment_block(block, "#", indent))
                    continue
            out.append(line)
            i += 1
            continue

        # Inside a docstring.
        if doc_quote in line:
            out.extend(_unwrap_doc(doc_body, doc_indent))
            out.append(line)
            in_doc = False
            doc_body = []
            i += 1
            continue
        doc_body.append(line)
        i += 1

    return "\n".join(out)


def _closes_on_same_line(line: str, quote: str) -> bool:
    return line.count(quote) >= 2


def _unwrap_doc(body: list[str], indent: str) -> list[str]:
    """Join paragraphs inside a docstring, preserving indented and listed text."""
    result: list[str] = []
    buffer: list[str] = []

    def flush():
        if buffer:
            result.append(indent + " ".join(buffer))
            buffer.clear()

    for raw in body:
        text = raw.strip()
        if not text:
            flush()
            result.append("")
            continue
        # Indented lines are examples or code; underlines mark section headings; list markers carry structure.
        if raw.startswith(indent + "    ") or re.fullmatch(r"[-=~]{3,}", text) \
                or re.match(r"^([-*+]\s|\d+[.)]\s|\|)", text):
            flush()
            result.append(raw)
            continue
        buffer.append(text)

    flush()
    return result


def process(path: str) -> bool:
    with open(path) as handle:
        original = handle.read()
    if path.endswith(".md"):
        updated = unwrap_markdown(original)
    elif path.endswith((".f90", ".F90")):
        updated = unwrap_fortran(original)
    elif path.endswith(".py"):
        updated = unwrap_python(original)
    else:
        return False
    if updated != original:
        with open(path, "w") as handle:
            handle.write(updated)
        return True
    return False


if __name__ == "__main__":
    changed = [p for p in sys.argv[1:] if process(p)]
    for p in changed:
        print("unwrapped", p)
    print(f"{len(changed)} file(s) changed")
