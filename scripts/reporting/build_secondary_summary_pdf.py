#!/usr/bin/env python3
"""Build a Mac-friendly PDF summary of the secondary-sample fitting outputs."""

from __future__ import annotations

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    KeepTogether,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "output" / "pdf" / "Secondary_Sample_Peak_Fitting_Summary.pdf"

NAVY = colors.HexColor("#17365D")
BLUE = colors.HexColor("#1F4E78")
PALE_BLUE = colors.HexColor("#D9EAF7")
ORANGE = colors.HexColor("#C65911")
GREEN = colors.HexColor("#548235")
PALE_GOLD = colors.HexColor("#FFF2CC")
GREY = colors.HexColor("#44546A")
LIGHT_GREY = colors.HexColor("#E7E6E6")


def header_footer(canvas, doc) -> None:
    canvas.saveState()
    width, height = A4
    canvas.setFillColor(NAVY)
    canvas.rect(0, height - 14 * mm, width, 14 * mm, fill=1, stroke=0)
    canvas.setFillColor(colors.white)
    canvas.setFont("Helvetica-Bold", 9)
    canvas.drawString(17 * mm, height - 9 * mm, "GIWAXS secondary-sample peak fitting")
    canvas.setFillColor(colors.HexColor("#666666"))
    canvas.setFont("Helvetica", 8)
    canvas.drawString(17 * mm, 10 * mm, "Selected-frame screening report - generated 4 August 2026")
    canvas.drawRightString(width - 17 * mm, 10 * mm, f"Page {doc.page}")
    canvas.restoreState()


def make_styles() -> dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "Title",
            parent=base["Title"],
            fontName="Helvetica-Bold",
            fontSize=23,
            leading=27,
            textColor=NAVY,
            alignment=TA_CENTER,
            spaceAfter=8 * mm,
        ),
        "subtitle": ParagraphStyle(
            "Subtitle",
            parent=base["Normal"],
            fontName="Helvetica",
            fontSize=12,
            leading=16,
            textColor=GREY,
            alignment=TA_CENTER,
            spaceAfter=6 * mm,
        ),
        "h1": ParagraphStyle(
            "H1",
            parent=base["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=17,
            leading=21,
            textColor=NAVY,
            spaceBefore=3 * mm,
            spaceAfter=4 * mm,
        ),
        "h2": ParagraphStyle(
            "H2",
            parent=base["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=12,
            leading=15,
            textColor=BLUE,
            spaceBefore=2 * mm,
            spaceAfter=2 * mm,
        ),
        "body": ParagraphStyle(
            "Body",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=9.5,
            leading=13.2,
            textColor=colors.HexColor("#222222"),
            alignment=TA_LEFT,
            spaceAfter=2.5 * mm,
        ),
        "small": ParagraphStyle(
            "Small",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=8.2,
            leading=11,
            textColor=colors.HexColor("#333333"),
            spaceAfter=1.5 * mm,
        ),
        "path": ParagraphStyle(
            "Path",
            parent=base["BodyText"],
            fontName="Courier",
            fontSize=7.1,
            leading=9.5,
            textColor=colors.HexColor("#333333"),
            wordWrap="CJK",
        ),
        "callout": ParagraphStyle(
            "Callout",
            parent=base["BodyText"],
            fontName="Helvetica-Bold",
            fontSize=10.2,
            leading=14,
            textColor=NAVY,
            alignment=TA_LEFT,
        ),
    }


def bullet(text: str, styles: dict[str, ParagraphStyle]) -> Paragraph:
    return Paragraph(f"&#8226;&nbsp;&nbsp;{text}", styles["body"])


def callout(text: str, styles: dict[str, ParagraphStyle], fill=PALE_BLUE) -> Table:
    table = Table([[Paragraph(text, styles["callout"])]], colWidths=[174 * mm])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), fill),
                ("BOX", (0, 0), (-1, -1), 0.8, BLUE),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
            ]
        )
    )
    return table


def sample_section(
    scan: str,
    sample: str,
    frames: str,
    count: int,
    outcome: list[str],
    report_dir: str,
    styles: dict[str, ParagraphStyle],
) -> KeepTogether:
    rows = [
        [Paragraph("Selected frames", styles["small"]), Paragraph(frames, styles["small"])],
        [Paragraph("Final q-position rows", styles["small"]), Paragraph(str(count), styles["small"])],
        [Paragraph("Report folder", styles["small"]), Paragraph(report_dir, styles["path"])],
    ]
    table = Table(rows, colWidths=[42 * mm, 128 * mm])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (0, -1), LIGHT_GREY),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#BFBFBF")),
                ("INNERGRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#D9D9D9")),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
            ]
        )
    )
    content = [
        Paragraph(f"Scan {scan} - {sample}", styles["h1"]),
        table,
        Spacer(1, 3 * mm),
    ]
    content.extend(bullet(item, styles) for item in outcome)
    return KeepTogether(content)


def build() -> None:
    styles = make_styles()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc = BaseDocTemplate(
        str(OUTPUT),
        pagesize=A4,
        rightMargin=18 * mm,
        leftMargin=18 * mm,
        topMargin=21 * mm,
        bottomMargin=18 * mm,
        title="Secondary Sample Peak Fitting Summary",
        author="GIWAXS dissertation project",
        subject="Selected-frame peak fitting for scans 587217, 587225, 587241 and 587250",
    )
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="main")
    doc.addPageTemplates([PageTemplate(id="standard", frames=[frame], onPage=header_footer)])

    story = [
        Spacer(1, 14 * mm),
        Paragraph("Secondary-Sample Peak Fitting Summary", styles["title"]),
        Paragraph(
            "Mac-friendly consolidated version of the completion index and final QC notes",
            styles["subtitle"],
        ),
        callout(
            "586875 was excluded at the user's request. Scans 587217, 587225, "
            "587241 and 587250 are complete for selected-frame FR/IP/OOP peak screening.",
            styles,
        ),
        Spacer(1, 8 * mm),
        Paragraph("Scientific scope", styles["h1"]),
        bullet("Constrained Gaussian, Lorentzian and pseudo-Voigt alternatives were compared in sample-specific windows.", styles),
        bullet("Only the selected frames were fitted. No all-frame expansion or MCMC was used.", styles),
        bullet("Bootstrap was set to 0 for every cut, so reported confidence intervals are covariance approximations.", styles),
        bullet("Use the results for peak presence, fitted centre q0, d-spacing and cautious area trends.", styles),
        bullet("Do not use apparent FWHM from these supporting samples for CCL, crystallite size, strain or disorder.", styles),
        Spacer(1, 5 * mm),
        Paragraph("What to open", styles["h1"]),
        Paragraph(
            "The Excel files named <b>scan&lt;scan&gt;_simple_peaks_per_frame.xlsx</b> are the quickest summaries. "
            "For dissertation tables, use <b>accepted_peak_positions_after_visual_qc.csv</b>. Fit overlays and "
            "residuals remain beside each report in the key_frame_fits folder.",
            styles["body"],
        ),
        Paragraph("Base folder", styles["h2"]),
        Paragraph(str(ROOT), styles["path"]),
        PageBreak(),
        Paragraph("Completion overview", styles["h1"]),
    ]

    overview_data = [
        ["Scan", "Sample", "Frames", "Final q rows", "Visual QC"],
        ["586875", "nogas first run", "Excluded", "-", "Not processed"],
        ["587217", "bcnogas1", "10", "89", "1 row removed"],
        ["587225", "n2n2", "10", "135", "3 rows removed"],
        ["587241", "airn2", "13", "150", "No extra removals"],
        ["587250", "airn2ITO", "8", "159", "No extra removals"],
    ]
    overview = Table(overview_data, colWidths=[23 * mm, 42 * mm, 24 * mm, 28 * mm, 53 * mm], repeatRows=1)
    overview.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), NAVY),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
                ("FONTSIZE", (0, 0), (-1, -1), 8.5),
                ("ALIGN", (0, 0), (0, -1), "CENTER"),
                ("ALIGN", (2, 1), (3, -1), "CENTER"),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE_BLUE]),
                ("BOX", (0, 0), (-1, -1), 0.7, NAVY),
                ("INNERGRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#BFBFBF")),
                ("TOPPADDING", (0, 0), (-1, -1), 7),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 7),
            ]
        )
    )
    story.extend(
        [
            overview,
            Spacer(1, 7 * mm),
            callout(
                "All requested FR/IP/OOP outputs contain exactly the chosen frames, all optimizer calls "
                "succeeded, and every final run records bootstrap 0.",
                styles,
                fill=PALE_GOLD,
            ),
            Spacer(1, 7 * mm),
            Paragraph("Workbook symbols", styles["h1"]),
            bullet("Dagger: centre/presence only; do not use its fitted width for CCL.", styles),
            bullet("Double dagger: near the conservative Delta BIC or area-SNR reporting threshold.", styles),
            bullet("Broad: a fitted broad component whose centre is sensitive to the local background.", styles),
            bullet("A dash means no independent q0 passed the gate; it is not proof that a physical feature is absent.", styles),
            PageBreak(),
            sample_section(
                "587217",
                "bcnogas1",
                "0, 1, 2, 3, 4, 5, 6, 11, 25, 35",
                89,
                [
                    "Frame 35 is late cooling at 73.9 degrees C, not a final high-temperature frame.",
                    "One OOP frame-1 broad shoulder near q = 1.73 A<super>-1</super> was removed because its centre was background/model-sensitive.",
                    "The exploratory FR feature near q = 2.060 A<super>-1</super> did not provide an independently accepted centre.",
                    "Under-resolved sharp lines remain q/presence screens only; no width or CCL claim is allowed.",
                ],
                "results/dimitar_bcnogas1/scan_587217/key_frame_report_v2/",
                styles,
            ),
            Spacer(1, 6 * mm),
            sample_section(
                "587225",
                "n2n2",
                "0, 10, 19, 21, 24, 27, 30, 40, 70, 134",
                135,
                [
                    "Reference-like underlays were kept only to stabilize neighbouring peak fits and were excluded from assignment/result tables.",
                    "Three fitted positions were removed after visual QC: FR frame 27 near 0.941 A<super>-1</super>, OOP frame 24 secondary 1.75 A<super>-1</super>, and OOP frame 70 near 0.835 A<super>-1</super>.",
                    "Thirty explicitly position-only rows may support q positions but not fitted widths or areas.",
                ],
                "results/dimitar_n2n2/scan_587225/key_frame_report_v3/",
                styles,
            ),
            PageBreak(),
            sample_section(
                "587241",
                "airn2",
                "0, 9, 18, 20, 25, 30, 35, 40, 45, 50, 60, 80, 129",
                150,
                [
                    "The selected frames bracket the main transformation between 48.8 and 61.6 degrees C.",
                    "No additional positions were removed by visual QC; absent/noisy instances were already rejected by the numerical gates.",
                    "The FR/OOP line near 0.742 A<super>-1</super> is under-resolved and is centre/presence only.",
                    "Treat early FR features near 1.51 and 1.96 A<super>-1</super> and OOP features near 1.48 and 2.19 A<super>-1</super> as conservative q observations because symmetric profiles leave structured residuals.",
                ],
                "results/dimitar_airn2/scan_587241/key_frame_report_v1/",
                styles,
            ),
            Spacer(1, 6 * mm),
            sample_section(
                "587250",
                "airn2ITO",
                "0, 11, 22, 30, 40, 75, 100, 129",
                159,
                [
                    "No additional q-position rows were removed by the independent visual audit.",
                    "The strong high-q lines in FR, IP and OOP are visibly trackable but asymmetric. Use their centres only as descriptive q positions.",
                    "Do not interpret the widths, areas or CCL of those high-q lines; substrate/reference origin remains possible until assignment.",
                    "The broad IP component near 0.290 A<super>-1</super> at frame 75 is near threshold and must remain provisional rather than being described as an evolving peak.",
                ],
                "results/dimitar_airn2ito/scan_587250/key_frame_report_v3/",
                styles,
            ),
            PageBreak(),
            Paragraph("Recommended dissertation use", styles["h1"]),
            Paragraph(
                "These selected-frame results are suitable for drafting the Results and Discussion chapter after peak assignment. "
                "Use the Excel workbook to identify which fitted centres occur in each frame and sector, then copy only the "
                "visually reviewed values from the accepted-position CSV into dissertation tables.",
                styles["body"],
            ),
            bullet("Describe whether a feature is present, absent from the accepted table, emerging, disappearing or shifting beyond uncertainty.", styles),
            bullet("Compare FR, IP and OOP positions and intensities cautiously; sector-dependent intensity supports orientation, while peak assignment requires structural context.", styles),
            bullet("Do not call every fitted component a crystallographic reflection before assignment.", styles),
            bullet("Do not calculate CCL from these supporting-sample widths. The DROP40 calibration/width workflow is a separate analysis.", styles),
            Spacer(1, 7 * mm),
            callout(
                "Supervisor-ready question: Are these peak-presence and q-position trends acceptable for the dissertation, "
                "and do the proposed peak assignments match the expected material phases and orientation?",
                styles,
                fill=PALE_GOLD,
            ),
            Spacer(1, 8 * mm),
            Paragraph("Files represented in this PDF", styles["h1"]),
            Paragraph(
                "This document consolidates the Markdown completion index and the four final QC notes. The CSV tables, "
                "Excel summaries, fit overlays, residual panels and full audit tables remain the authoritative analysis files.",
                styles["body"],
            ),
        ]
    )
    doc.build(story)
    print(OUTPUT)


if __name__ == "__main__":
    build()

