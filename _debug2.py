import fitz, re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

pdf_path = r"mdsp-i.pdf"
doc = fitz.open(pdf_path)

SIZE_TO_LEVEL = {
    27.5: 1, 17.2: 2, 12.9: 3, 10.1: 4,
}

all_headings = []

for page_num in range(doc.page_count):
    page = doc[page_num]
    blocks = page.get_text("dict")["blocks"]
    page_spans = []

    for block in blocks:
        if "lines" not in block:
            continue
        for line in block["lines"]:
            for span in line["spans"]:
                size = round(span["size"], 1)
                text = span["text"].strip()
                if not text:
                    continue
                level = SIZE_TO_LEVEL.get(size)
                if level is not None:
                    page_spans.append((level, size, text, span["bbox"]))

    # Show first few pages with L1/L2 matches
    if page_num <= 10 and page_spans:
        l1 = [(l,s,t) for l,s,t,_ in page_spans if l<=2]
        if l1:
            print(f"p{page_num+1} raw spans: {l1}")

    merged = []
    i = 0
    while i < len(page_spans):
        level, size, text, bbox = page_spans[i]
        x0, y0, x1, y1 = bbox
        j = i + 1
        while j < len(page_spans) and page_spans[j][0] == level and page_spans[j][1] == size:
            bx0, by0, bx1, by1 = page_spans[j][3]
            if abs(by0 - y0) >= 6:
                break
            text += page_spans[j][2]
            x1 = max(x1, bx1)
            y1 = max(y1, by1)
            j += 1
        text = re.sub(r'\s+', ' ', text).strip()
        merged.append((page_num + 1, level, size, text, (x0, y0, x1, y1)))
        i = j

    for page, level, size, text, bbox in merged:
        if len(text) < 2:
            continue
        if re.match(r'^[\s\.\-·•/\$)（）]+$', text):
            continue
        if level in (1, 2, 3):
            all_headings.append((page, level, text, bbox))
        elif level == 4:
            if re.match(r'^\d+\.\d+', text):
                all_headings.append((page, level, text, bbox))

print(f"\nTotal headings found: {len(all_headings)}")
if all_headings:
    for lv in [1,2,3,4]:
        items = [(p,t) for p,l,t,_ in all_headings if l==lv]
        print(f"  Level {lv}: {len(items)} items")
        if items[:5]:
            for p,t in items[:5]:
                print(f"    p{p}: {t[:60]}")
else:
    print("  EMPTY!")
    
doc.close()
