import os
import re

def _build_sphinx_paper():
    """
    Convert paper.md (JOSS/pandoc format) to paper_sphinx.md (MyST format)
    so citations render correctly on the docs site without touching paper.md.

    Conversions applied:
      [@key1; @key2]  ->  {cite:p}`key1,key2`   (parenthetical)
      @key            ->  {cite:t}`key`           (inline author-year)
    """
    src = os.path.join(os.path.dirname(__file__), "../../paper.md")
    dst = os.path.join(os.path.dirname(__file__), "paper_sphinx.md")

    with open(src) as f:
        text = f.read()

    def _multi_cite(m):
        keys = [k.strip().lstrip("@") for k in m.group(1).split(";")]
        return "{cite:p}`" + ",".join(keys) + "`"

    # bracketed citations first: [@key] or [@key1; @key2]
    text = re.sub(r"\[@([^\]]+)\]", _multi_cite, text)
    # remaining bare @key (inline author citations)
    text = re.sub(r"(?<!\[)@([A-Za-z]\w*)", r"{cite:t}`\1`", text)

    # Only write if content changed — prevents sphinx-autobuild reload loop
    if os.path.exists(dst):
        with open(dst) as f:
            if f.read() == text:
                return

    with open(dst, "w") as f:
        f.write(text)
