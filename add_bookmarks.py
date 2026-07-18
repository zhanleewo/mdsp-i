"""从 mdsp-i.pdf 提取目录并添加书签（精确跳转至标题位置）"""
import fitz
import sys
import io
import re
from collections import defaultdict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

pdf_path = r"mdsp-i.pdf"
output_path = r"现代数字信号处理 I 笔记.pdf"

doc = fitz.open(pdf_path)

SIZE_TO_LEVEL = {
    33.2: 1,   # H1: 第x单元
    20.7: 2,   # H2: 第x讲
    15.5: 3,   # H3: x. xxx
    12.1: 4,   # H4: x.x xxx
}

# ---- 诊断模式：扫描全部字体大小 ----
SCAN_MODE = "--scan" in sys.argv
if SCAN_MODE:
    print(f"\n{'='*60}")
    print("PDF 字体大小扫描（诊断模式）")
    print(f"{'='*60}")
    size_samples = defaultdict(set)
    for page_num in range(min(doc.page_count, 30)):
        page = doc[page_num]
        blocks = page.get_text("dict")["blocks"]
        for block in blocks:
            if "lines" not in block:
                continue
            for line in block["lines"]:
                for span in line["spans"]:
                    size = round(span["size"], 1)
                    text = span["text"].strip()
                    if not text or len(text) < 2:
                        continue
                    if len(size_samples[size]) < 5:
                        size_samples[size].add(text[:50])
    print(f"\n{'字体大小':>8s} | {'示例文本'}")
    print("-" * 60)
    for sz in sorted(size_samples, reverse=True):
        samples = " | ".join(sorted(size_samples[sz])[:4])
        print(f"{sz:>8.1f}pt | {samples}")
    print(f"{'='*60}\n")
    print("请根据上述扫描结果更新 SIZE_TO_LEVEL 字典。")
    print("H1 = 单元标题 (第x单元)")
    print("H2 = 讲标题 (第x讲)")
    print("H3 = 小节标题 (x. xxx)")
    print("H4 = 子小节标题 (x.x xxx)\n")
    doc.close()
    sys.exit(0)

all_headings = []

# ---- 第一步：提取所有标题及其精确位置 ----
for page_num in range(doc.page_count):
    page = doc[page_num]
    page_rect = page.rect
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

    # 合并同一行内、同字号、同level的相邻span（只合并同行，不同行不合并）
    merged = []
    i = 0
    while i < len(page_spans):
        level, size, text, bbox = page_spans[i]
        x0, y0, x1, y1 = bbox
        j = i + 1
        while j < len(page_spans) and page_spans[j][0] == level and page_spans[j][1] == size:
            bx0, by0, bx1, by1 = page_spans[j][3]
            # 只在同一行内合并（y0 差距小于 6pt 视为同一行）
            if abs(by0 - y0) >= 6:
                break
            text += page_spans[j][2]
            x1 = max(x1, bx1)
            y1 = max(y1, by1)
            j += 1
        text = re.sub(r'\s+', ' ', text).strip()
        merged.append((page_num + 1, level, size, text, (x0, y0, x1, y1)))
        i = j

    # 过滤
    for page, level, size, text, bbox in merged:
        if len(text) < 2:
            continue
        # 纯标点/分隔符
        if re.match(r'^[\s\.\-·•/\\()（）]+$', text):
            continue
        # L1/L2/L3: 保留所有标题
        if level in (1, 2, 3):
            all_headings.append((page, level, text, bbox))
        # L4 (12.1pt): 只保留编号子小节（如 "1.1", "2.3"），过滤正文混入
        elif level == 4:
            if re.match(r'^\d+\.\d+', text):
                all_headings.append((page, level, text, bbox))

# ---- 第二步：构建 set_toc 书签列表 ----
toc_entries = []

if all_headings:
    # 插入虚拟 H1 作为引言等前置内容的父级
    first_page = all_headings[0][0]
    first_x0, first_y0 = all_headings[0][3][0], all_headings[0][3][1]
    toc_entries.append([1, "课程导论 / 前言", first_page, {
        "kind": fitz.LINK_GOTO,
        "page": first_page - 1,
        "to": fitz.Point(first_x0, first_y0),
    }])

    for page, level, text, (x0, y0, x1, y1) in all_headings:
        dest = {
            "kind": fitz.LINK_GOTO,
            "page": page - 1,
            "to": fitz.Point(x0, y0),
        }
        toc_entries.append([level, text, page, dest])

    doc.set_toc(toc_entries)
    doc.save(output_path, garbage=4, deflate=True)
else:
    print("警告: 未找到任何标题！请检查 SIZE_TO_LEVEL 字体大小配置。")

doc.close()

# ---- 第三步：输出结果 ----
print(f"原文件: {pdf_path}")
print(f"输出文件: {output_path}")
print(f"书签总数: {len(toc_entries)}")
if toc_entries:
    print(f"  H1 (单元): {sum(1 for e in toc_entries if e[0]==1)}")
    print(f"  H2 (讲):   {sum(1 for e in toc_entries if e[0]==2)}")
    print(f"  H3 (小节): {sum(1 for e in toc_entries if e[0]==3)}")
    print(f"  H4 (子小节): {sum(1 for e in toc_entries if e[0]==4)}")
print(f"\n每个书签均包含精确的页面坐标，点击后可跳转到标题所在位置。")
