import fitz, re
doc = fitz.open(r"mdsp-i.pdf")

for pn in [2, 4, 7, 55, 170]:  # 0-indexed
    page = doc[pn]
    blocks = page.get_text("dict")["blocks"]
    print(f"=== Page {pn+1} ({len(blocks)} blocks) ===")
    for bi, b in enumerate(blocks):
        if "lines" not in b:
            print(f"  block{bi}: no lines, type={b.get('type','?')}")
            continue
        for li, l in enumerate(b["lines"]):
            sizes = [round(s["size"], 1) for s in l["spans"]]
            texts = [s["text"].strip() for s in l["spans"] if s["text"].strip()]
            if any(sz >= 15 for sz in sizes) or "单元" in "".join(texts) or "讲" in "".join(texts):
                print(f"  block{bi} line{li}: sizes={sizes}")
                print(f"    texts={texts}")
    print()

doc.close()
