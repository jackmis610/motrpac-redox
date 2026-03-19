#!/usr/bin/env python3
"""
16_sensitivity_table_pdf.py
Produces figures/table1_sensitivity.pdf — sensitivity analysis table
using reportlab.
"""

import os
from reportlab.lib import colors
from reportlab.lib.pagesizes import landscape, letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
)
from reportlab.lib.enums import TA_LEFT, TA_CENTER

os.makedirs("figures", exist_ok=True)

# ── Page setup ───────────────────────────────────────────────────────────────
PAGE = landscape(letter)   # 11 × 8.5 inches
MARGIN = 0.65 * inch

doc = SimpleDocTemplate(
    "figures/table1_sensitivity.pdf",
    pagesize=PAGE,
    leftMargin=MARGIN,
    rightMargin=MARGIN,
    topMargin=0.6 * inch,
    bottomMargin=0.6 * inch,
)

styles = getSampleStyleSheet()

title_style = ParagraphStyle(
    "Title",
    parent=styles["Normal"],
    fontSize=11,
    fontName="Helvetica-Bold",
    leading=15,
    spaceBefore=0,
    spaceAfter=10,
    alignment=TA_LEFT,
)

section_style = ParagraphStyle(
    "Section",
    parent=styles["Normal"],
    fontSize=10,
    fontName="Helvetica-Bold",
    leading=13,
    spaceBefore=6,
    spaceAfter=4,
    alignment=TA_LEFT,
)

footnote_style = ParagraphStyle(
    "Footnote",
    parent=styles["Normal"],
    fontSize=7.5,
    fontName="Helvetica",
    leading=11,
    spaceBefore=8,
    spaceAfter=0,
    alignment=TA_LEFT,
    textColor=colors.HexColor("#333333"),
)

# ── Colours ──────────────────────────────────────────────────────────────────
HEADER_BG    = colors.HexColor("#2166AC")   # blue header
HEADER_FG    = colors.white
CONSERV_BG   = colors.HexColor("#D6E4F5")   # light blue highlight
ALT_BG       = colors.HexColor("#F5F8FC")   # alternating row
WHITE        = colors.white
BLACK        = colors.black
RULE_COLOR   = colors.HexColor("#AAAAAA")

# ── Table data ───────────────────────────────────────────────────────────────
HEADER = ["Scenario", "n", "\u03b2", "95% CI", "R\u00b2", "p(\u03b2<1)"]

MALE_ROWS = [
    ["All tissues (full model)",              "18", "0.651", "[0.434, 0.868]", "0.716", "0.0018"],
    ["Minus BAT",                              "17", "1.220", "[0.943, 1.497]", "0.855", "0.945"],
    ["Minus BLOOD",                            "17", "0.446", "[0.364, 0.528]", "0.901", "1.6\u00d710\u207b\u00b9\u2070"],
    ["Minus BAT + BLOOD (conservative)",       "16", "0.551", "[0.279, 0.824]", "0.574", "0.0017"],
    ["Minus top-3 (BAT+BLOOD+SKM-GN)",         "15", "0.689", "[0.384, 0.995]", "0.646", "0.023"],
]

FEMALE_ROWS = [
    ["All tissues (full model)",              "18", "0.150", "[0.038, 0.261]", "0.333", "1.4\u00d710\u207b\u00b9\u00b9"],
    ["Minus BAT",                              "17", "0.115", "[0.033, 0.196]", "0.280", "4.7\u00d710\u207b\u00b9\u00b2"],
]

# ── Build a styled table ──────────────────────────────────────────────────────
def make_table(header, rows, conservative_row_idx=None):
    """
    conservative_row_idx: 0-based index in `rows` (not counting header)
                          that gets the highlight.
    """
    data = [header] + rows
    n_cols = len(header)
    n_rows = len(data)

    # Column widths
    col_widths = [3.8*inch, 0.45*inch, 0.55*inch, 1.35*inch, 0.55*inch, 1.05*inch]

    tbl = Table(data, colWidths=col_widths, repeatRows=1)

    style_cmds = [
        # Global font
        ("FONTNAME",    (0, 0), (-1, -1), "Helvetica"),
        ("FONTSIZE",    (0, 0), (-1, -1), 9),
        ("TOPPADDING",  (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("VALIGN",      (0, 0), (-1, -1), "MIDDLE"),

        # Header row
        ("BACKGROUND",  (0, 0), (-1, 0), HEADER_BG),
        ("TEXTCOLOR",   (0, 0), (-1, 0), HEADER_FG),
        ("FONTNAME",    (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",    (0, 0), (-1, 0), 9),
        ("ALIGN",       (0, 0), (-1, 0), "CENTER"),

        # Data rows: left-align scenario, center numerics
        ("ALIGN", (0, 1), (0, -1), "LEFT"),
        ("ALIGN", (1, 1), (-1, -1), "CENTER"),

        # Alternating rows
        *[
            ("BACKGROUND", (0, i), (-1, i), ALT_BG if i % 2 == 0 else WHITE)
            for i in range(1, n_rows)
        ],

        # Grid
        ("GRID",        (0, 0), (-1, -1), 0.4, RULE_COLOR),
        ("LINEBELOW",   (0, 0), (-1, 0),  1.0, colors.HexColor("#FFFFFF")),
    ]

    # Conservative row highlight
    if conservative_row_idx is not None:
        r = conservative_row_idx + 1   # +1 for header offset
        style_cmds += [
            ("BACKGROUND",  (0, r), (-1, r), CONSERV_BG),
            ("FONTNAME",    (0, r), (-1, r), "Helvetica-BoldOblique"),
        ]

    tbl.setStyle(TableStyle(style_cmds))
    return tbl

# ── Assemble document ─────────────────────────────────────────────────────────
title_text = (
    "Table 1. Sensitivity analysis: ETS\u2192buffering regression slope "
    "across tissue exclusion scenarios (8-week TRNSCRPT, males and females)"
)

footnote_text = (
    "\u03b2: OLS regression slope (buffering ~ ETS, cross-tissue). "
    "p(\u03b2<1): one-tailed t-test vs. slope=1 null. "
    "ETS index: mean log\u2082FC of 31 ETS+ATP synthase genes. "
    "Buffering index: mean log\u2082FC of 15 GSH/TRX/NNT/SOD genes. "
    "BAT excluded from conservative estimate due to multi-feature mapping artifact (Cook\u2019s D=16.7). "
    "BLOOD is the highest-leverage legitimate outlier (Cook\u2019s D=2.54)."
)

story = [
    Paragraph(title_text, title_style),
    Spacer(1, 0.08 * inch),
    Paragraph("Males", section_style),
    make_table(HEADER, MALE_ROWS, conservative_row_idx=3),
    Spacer(1, 0.18 * inch),
    Paragraph("Females", section_style),
    make_table(HEADER, FEMALE_ROWS, conservative_row_idx=None),
    Spacer(1, 0.12 * inch),
    Paragraph(footnote_text, footnote_style),
]

doc.build(story)
print("Saved figures/table1_sensitivity.pdf")
