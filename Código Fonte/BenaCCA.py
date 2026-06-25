"""BenaCCA: projeto de crossover passivo Butterworth.

O programa reúne em um único arquivo:
- cálculos dos filtros LPF e HPF;
- seleção dos componentes comerciais;
- geração dos gráficos de Bode;
- exportação de relatório em PDF;
- interface gráfica em Tkinter.
"""

from __future__ import annotations

import ctypes
import os
import sys
import tkinter as tk
from dataclasses import dataclass
from datetime import datetime
from math import isfinite, pi, sqrt
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Literal

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Image,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


class BodeNavigationToolbar(NavigationToolbar2Tk):
    """Barra do gráfico sem o ajuste manual de subplots."""

    toolitems = tuple(
        item
        for item in NavigationToolbar2Tk.toolitems
        if item[3] != "configure_subplots"
    )


# =============================================================================
# TABELAS DE COMPONENTES COMERCIAIS
# =============================================================================

FilterKind = Literal["LPF", "HPF"]
BodeTheme = Literal["light", "dark"]

# Tabela 1 do enunciado. Os valores escritos em mH são convertidos para H.
COMMERCIAL_INDUCTORS_H = np.array(
    [
        0.10,
        0.12,
        0.15,
        0.18,
        0.22,
        0.27,
        0.33,
        0.39,
        0.47,
        0.56,
        0.68,
        0.82,
        1.0,
        1.2,
        1.5,
        1.8,
        2.2,
        2.7,
        3.3,
        3.9,
        4.7,
        5.6,
        6.8,
        8.2,
        10.0,
        12.0,
        15.0,
    ],
    dtype=float,
) * 1e-3

# Tabela 2 do enunciado. Os valores escritos em uF são convertidos para F.
COMMERCIAL_CAPACITORS_F = np.array(
    [
        1.0,
        1.2,
        1.5,
        1.8,
        2.2,
        2.7,
        3.3,
        3.9,
        4.7,
        5.6,
        6.8,
        8.2,
        10.0,
        12.0,
        15.0,
        18.0,
        22.0,
        27.0,
        33.0,
        39.0,
        47.0,
        56.0,
        68.0,
        82.0,
        100.0,
    ],
    dtype=float,
) * 1e-6


# =============================================================================
# ESTRUTURAS DE DADOS
# =============================================================================


@dataclass(frozen=True)
class ComponentChoice:
    """Armazena o valor ideal e o valor comercial de um componente."""

    symbol: str
    ideal: float
    commercial: float

    @property
    def absolute_error(self) -> float:
        return abs(self.commercial - self.ideal)

    @property
    def percent_error(self) -> float:
        return 100.0 * self.absolute_error / self.ideal


@dataclass(frozen=True)
class FilterResult:
    """Resultados de um dos filtros do crossover."""

    kind: FilterKind
    name: str
    destination: str
    inductor: ComponentChoice
    capacitor: ComponentChoice
    ideal_cutoff_hz: float
    commercial_cutoff_hz: float
    commercial_natural_hz: float
    commercial_q: float
    max_magnitude_difference_db: float
    tolerance_cutoff_min_hz: float
    tolerance_cutoff_max_hz: float

    @property
    def cutoff_error_hz(self) -> float:
        return self.commercial_cutoff_hz - self.ideal_cutoff_hz

    @property
    def cutoff_error_percent(self) -> float:
        return 100.0 * self.cutoff_error_hz / self.ideal_cutoff_hz


@dataclass(frozen=True)
class ProjectResult:
    """Reúne todos os resultados calculados pelo programa."""

    cutoff_hz: float
    load_ohm: float
    tolerance_percent: float
    lpf: FilterResult
    hpf: FilterResult

    @property
    def audible_assessment(self) -> str:
        """Produz uma avaliação cuidadosa do possível impacto audível."""

        largest_shift = max(
            abs(self.lpf.cutoff_error_percent),
            abs(self.hpf.cutoff_error_percent),
        )
        largest_curve_difference = max(
            self.lpf.max_magnitude_difference_db,
            self.hpf.max_magnitude_difference_db,
        )
        if largest_shift < 3.0 and largest_curve_difference < 0.75:
            return (
                "As diferenças são pequenas e tendem a ser pouco perceptíveis em "
                "uma escuta comum, embora possam ser identificadas por medição."
            )
        if largest_shift < 8.0 and largest_curve_difference < 1.5:
            return (
                "A alteração tende a ser sutil, mas pode mudar levemente a região "
                "de transição entre o woofer e o tweeter."
            )
        return (
            "A alteração pode ser perceptível na região do crossover, modificando "
            "o equilíbrio entre o woofer e o tweeter. A confirmação exige medição "
            "dos alto-falantes reais e teste de escuta."
        )


# =============================================================================
# CÁLCULOS DO CROSSOVER
# =============================================================================


def parse_positive_number(text: str, field_name: str) -> float:
    """Converte números com ponto ou vírgula e exige um valor positivo."""

    normalized = text.strip().replace(" ", "").replace(",", ".")
    if not normalized:
        raise ValueError(f"Informe {field_name}.")
    try:
        value = float(normalized)
    except ValueError as exc:
        raise ValueError(f"{field_name.capitalize()} deve ser um número válido.") from exc
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{field_name.capitalize()} deve ser maior que zero.")
    return value


def ideal_components(cutoff_hz: float, load_ohm: float) -> tuple[float, float]:
    """Calcula L e C para um Butterworth passivo de segunda ordem.

    O denominador das topologias adotadas é igualado à forma normalizada:

        s²/wc² + sqrt(2)s/wc + 1

    Portanto:

        L = sqrt(2)R/wc
        C = 1/(sqrt(2)Rwc)
        wc = 2*pi*fc
    """

    if cutoff_hz <= 0 or load_ohm <= 0:
        raise ValueError("Frequência e impedância devem ser maiores que zero.")
    omega_c = 2.0 * pi * cutoff_hz
    inductance_h = sqrt(2.0) * load_ohm / omega_c
    capacitance_f = 1.0 / (sqrt(2.0) * load_ohm * omega_c)
    return inductance_h, capacitance_f


def nearest_commercial(ideal: float, available: np.ndarray) -> float:
    """Seleciona o valor comercial mais próximo do valor ideal."""

    distances = np.abs(available - ideal)
    minimum = float(np.min(distances))
    tied = available[np.isclose(distances, minimum, rtol=1e-12, atol=1e-18)]
    return float(np.min(tied))


def transfer_response(
    kind: FilterKind,
    frequencies_hz: np.ndarray,
    load_ohm: float,
    inductance_h: float,
    capacitance_f: float,
) -> np.ndarray:
    """Calcula a função de transferência complexa do filtro.

    LPF: indutor em série e capacitor em paralelo com a carga.
    HPF: capacitor em série e indutor em paralelo com a carga.

    Denominador comum:
        D(s) = R + Ls + RLCs²
    """

    s = 1j * 2.0 * pi * np.asarray(frequencies_hz, dtype=float)
    denominator = load_ohm + inductance_h * s + (
        load_ohm * inductance_h * capacitance_f * s**2
    )
    if kind == "LPF":
        numerator = load_ohm
    elif kind == "HPF":
        numerator = load_ohm * inductance_h * capacitance_f * s**2
    else:
        raise ValueError(f"Tipo de filtro desconhecido: {kind}")
    return numerator / denominator


def magnitude_db(response: np.ndarray) -> np.ndarray:
    """Converte a resposta complexa em magnitude expressa em dB."""

    return 20.0 * np.log10(np.maximum(np.abs(response), 1e-12))


def cutoff_from_response(
    kind: FilterKind,
    load_ohm: float,
    inductance_h: float,
    capacitance_f: float,
) -> float:
    """Localiza numericamente o ponto de -3,0103 dB do filtro real."""

    natural_hz = 1.0 / (2.0 * pi * sqrt(inductance_h * capacitance_f))
    frequencies = np.logspace(
        np.log10(natural_hz / 1000.0),
        np.log10(natural_hz * 1000.0),
        100_000,
    )
    db = magnitude_db(
        transfer_response(
            kind,
            frequencies,
            load_ohm,
            inductance_h,
            capacitance_f,
        )
    )
    target = -10.0 * np.log10(2.0)
    delta = db - target
    crossings = np.where(delta[:-1] * delta[1:] <= 0)[0]
    if crossings.size == 0:
        return natural_hz

    # Escolhe o cruzamento descendente no LPF e ascendente no HPF.
    candidates: list[int] = []
    for index in crossings:
        slope = db[index + 1] - db[index]
        if (kind == "LPF" and slope < 0) or (kind == "HPF" and slope > 0):
            candidates.append(int(index))
    index = min(
        candidates or [int(crossings[0])],
        key=lambda i: abs(np.log(frequencies[i] / natural_hz)),
    )

    # Interpolação no eixo logarítmico para melhorar a precisão.
    x1, x2 = np.log10(frequencies[index : index + 2])
    y1, y2 = db[index : index + 2]
    if y2 == y1:
        return float(frequencies[index])
    x_cutoff = x1 + (target - y1) * (x2 - x1) / (y2 - y1)
    return float(10.0**x_cutoff)


def max_curve_difference(
    kind: FilterKind,
    cutoff_hz: float,
    load_ohm: float,
    ideal_l: float,
    ideal_c: float,
    commercial_l: float,
    commercial_c: float,
) -> float:
    """Calcula a maior diferença entre as curvas perto do crossover."""

    frequencies = np.logspace(
        np.log10(cutoff_hz / 5.0),
        np.log10(cutoff_hz * 5.0),
        4000,
    )
    ideal_db = magnitude_db(
        transfer_response(kind, frequencies, load_ohm, ideal_l, ideal_c)
    )
    commercial_db = magnitude_db(
        transfer_response(
            kind,
            frequencies,
            load_ohm,
            commercial_l,
            commercial_c,
        )
    )
    return float(np.max(np.abs(ideal_db - commercial_db)))


def calculate_project(
    cutoff_hz: float,
    load_ohm: float,
    tolerance_percent: float = 5.0,
) -> ProjectResult:
    """Executa o dimensionamento completo do crossover."""

    if not isfinite(cutoff_hz) or not isfinite(load_ohm):
        raise ValueError("Frequência e impedância devem ser números finitos.")
    if cutoff_hz <= 0 or load_ohm <= 0:
        raise ValueError("Frequência e impedância devem ser maiores que zero.")
    if not isfinite(tolerance_percent) or not 0 <= tolerance_percent < 100:
        raise ValueError("A tolerância deve estar entre 0% e 100%.")

    ideal_l, ideal_c = ideal_components(cutoff_hz, load_ohm)
    commercial_l = nearest_commercial(ideal_l, COMMERCIAL_INDUCTORS_H)
    commercial_c = nearest_commercial(ideal_c, COMMERCIAL_CAPACITORS_F)

    inductor = ComponentChoice("L", ideal_l, commercial_l)
    capacitor = ComponentChoice("C", ideal_c, commercial_c)
    natural_hz = 1.0 / (2.0 * pi * sqrt(commercial_l * commercial_c))
    q_factor = load_ohm * sqrt(commercial_c / commercial_l)

    def build_filter(kind: FilterKind, name: str, destination: str) -> FilterResult:
        tolerance = tolerance_percent / 100.0
        tolerance_cutoffs = [
            cutoff_from_response(
                kind,
                load_ohm,
                commercial_l * l_factor,
                commercial_c * c_factor,
            )
            for l_factor in (1.0 - tolerance, 1.0 + tolerance)
            for c_factor in (1.0 - tolerance, 1.0 + tolerance)
        ]
        return FilterResult(
            kind=kind,
            name=name,
            destination=destination,
            inductor=inductor,
            capacitor=capacitor,
            ideal_cutoff_hz=cutoff_hz,
            commercial_cutoff_hz=cutoff_from_response(
                kind,
                load_ohm,
                commercial_l,
                commercial_c,
            ),
            commercial_natural_hz=natural_hz,
            commercial_q=q_factor,
            max_magnitude_difference_db=max_curve_difference(
                kind,
                cutoff_hz,
                load_ohm,
                ideal_l,
                ideal_c,
                commercial_l,
                commercial_c,
            ),
            tolerance_cutoff_min_hz=min(tolerance_cutoffs),
            tolerance_cutoff_max_hz=max(tolerance_cutoffs),
        )

    return ProjectResult(
        cutoff_hz=cutoff_hz,
        load_ohm=load_ohm,
        tolerance_percent=tolerance_percent,
        lpf=build_filter("LPF", "Passa-baixas", "Woofer"),
        hpf=build_filter("HPF", "Passa-altas", "Tweeter"),
    )


# =============================================================================
# GRÁFICOS DE BODE
# =============================================================================


def build_bode_figure(
    project: ProjectResult,
    figsize: tuple[float, float] = (12.5, 7.2),
    theme: BodeTheme = "light",
) -> Figure:
    """Cria o Bode completo: magnitude e fase dos dois filtros."""

    frequencies = np.logspace(
        np.log10(max(project.cutoff_hz / 20.0, 1.0)),
        np.log10(project.cutoff_hz * 20.0),
        1000,
    )
    palettes = {
        "light": {
            "figure": "white",
            "axes": "white",
            "text": "#172033",
            "muted": "#64748b",
            "grid": "#cbd5e1",
            "spine": "#cbd5e1",
            "legend": "white",
            "legend_edge": "#cbd5e1",
            "ideal": "#2563eb",
            "commercial": "#f97316",
        },
        "dark": {
            "figure": "#101114",
            "axes": "#181A20",
            "text": "#F3F4F6",
            "muted": "#A7B0C0",
            "grid": "#343946",
            "spine": "#3F4656",
            "legend": "#22252D",
            "legend_edge": "#343946",
            "ideal": "#60A5FA",
            "commercial": "#FDBA74",
        },
    }
    palette = palettes[theme]
    # Margens fixas mantêm o layout estável e deixam a interação mais leve.
    figure = Figure(figsize=figsize, dpi=100, facecolor=palette["figure"])
    figure.subplots_adjust(
        left=0.075,
        right=0.98,
        bottom=0.09,
        top=0.88,
        wspace=0.20,
        hspace=0.40,
    )
    axes = figure.subplots(2, 2, sharex="col")
    for axis in axes.flat:
        axis.set_facecolor(palette["axes"])
        axis.tick_params(axis="both", which="both", colors=palette["muted"])
        for spine in axis.spines.values():
            spine.set_color(palette["spine"])

    for row, result in enumerate((project.lpf, project.hpf)):
        ideal_response = transfer_response(
            result.kind,
            frequencies,
            project.load_ohm,
            result.inductor.ideal,
            result.capacitor.ideal,
        )
        commercial_response = transfer_response(
            result.kind,
            frequencies,
            project.load_ohm,
            result.inductor.commercial,
            result.capacitor.commercial,
        )
        ideal_magnitude = magnitude_db(ideal_response)
        commercial_magnitude = magnitude_db(commercial_response)
        ideal_phase = np.degrees(np.unwrap(np.angle(ideal_response)))
        commercial_phase = np.degrees(np.unwrap(np.angle(commercial_response)))

        magnitude_axis = axes[row, 0]
        phase_axis = axes[row, 1]

        ideal_magnitude_line = magnitude_axis.semilogx(
            frequencies,
            ideal_magnitude,
            color=palette["ideal"],
            linewidth=2.2,
            label="Componentes ideais",
        )[0]
        commercial_magnitude_line = magnitude_axis.semilogx(
            frequencies,
            commercial_magnitude,
            color=palette["commercial"],
            linewidth=2.0,
            linestyle="--",
            label="Componentes comerciais",
        )[0]
        ideal_magnitude_line.set_gid("bode:Ideal")
        commercial_magnitude_line.set_gid("bode:Comercial")
        magnitude_axis.axhline(
            -10.0 * np.log10(2.0),
            color=palette["muted"],
            linewidth=1.0,
            linestyle=":",
            label="-3,01 dB",
        )
        magnitude_axis.axvline(
            project.cutoff_hz, color=palette["ideal"], linewidth=1.0, alpha=0.65
        )
        magnitude_axis.axvline(
            result.commercial_cutoff_hz,
            color=palette["commercial"],
            linewidth=1.0,
            alpha=0.70,
        )
        magnitude_title = magnitude_axis.set_title(
            f"{result.name} — {result.destination} | Magnitude", loc="left"
        )
        magnitude_title.set_color(palette["text"])
        magnitude_axis.set_ylabel("Magnitude (dB)", color=palette["text"])
        magnitude_axis.set_ylim(-55, 5)
        magnitude_axis.grid(True, which="both", color=palette["grid"], alpha=0.55)
        magnitude_legend = magnitude_axis.legend(loc="best", fontsize=8)
        magnitude_legend.get_frame().set_facecolor(palette["legend"])
        magnitude_legend.get_frame().set_edgecolor(palette["legend_edge"])
        for label in magnitude_legend.get_texts():
            label.set_color(palette["text"])

        ideal_phase_line = phase_axis.semilogx(
            frequencies,
            ideal_phase,
            color=palette["ideal"],
            linewidth=2.2,
            label="Componentes ideais",
        )[0]
        commercial_phase_line = phase_axis.semilogx(
            frequencies,
            commercial_phase,
            color=palette["commercial"],
            linewidth=2.0,
            linestyle="--",
            label="Componentes comerciais",
        )[0]
        ideal_phase_line.set_gid("bode:Ideal")
        commercial_phase_line.set_gid("bode:Comercial")
        phase_axis.axvline(
            project.cutoff_hz, color=palette["ideal"], linewidth=1.0, alpha=0.65
        )
        phase_axis.axvline(
            result.commercial_cutoff_hz,
            color=palette["commercial"],
            linewidth=1.0,
            alpha=0.70,
        )
        phase_title = phase_axis.set_title(
            f"{result.name} — {result.destination} | Fase", loc="left"
        )
        phase_title.set_color(palette["text"])
        phase_axis.set_ylabel("Fase (graus)", color=palette["text"])
        phase_axis.grid(True, which="both", color=palette["grid"], alpha=0.55)
        phase_legend = phase_axis.legend(loc="best", fontsize=8)
        phase_legend.get_frame().set_facecolor(palette["legend"])
        phase_legend.get_frame().set_edgecolor(palette["legend_edge"])
        for label in phase_legend.get_texts():
            label.set_color(palette["text"])

    axes[1, 0].set_xlabel("Frequência (Hz)")
    axes[1, 1].set_xlabel("Frequência (Hz)")
    for axis in axes.flat:
        axis.xaxis.label.set_color(palette["text"])
        axis.yaxis.label.set_color(palette["text"])

    figure_title = figure.suptitle(
        f"BenaCCA — Bode completo ideal × comercial | {project.cutoff_hz:g} Hz, "
        f"{project.load_ohm:g} Ω",
        fontsize=14,
        fontweight="bold",
    )
    figure_title.set_color(palette["text"])
    return figure


def save_bode_figure(project: ProjectResult, output_path: str | Path) -> Path:
    """Salva o gráfico comparativo em PNG."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure = build_bode_figure(project)
    figure.savefig(
        output,
        dpi=180,
        bbox_inches="tight",
        facecolor=figure.get_facecolor(),
    )
    return output


# =============================================================================
# RELATÓRIO PDF
# =============================================================================


def format_pt(value: float, decimals: int = 3) -> str:
    """Formata um número usando vírgula decimal."""

    return f"{value:.{decimals}f}".replace(".", ",")


def filter_rows(result: FilterResult) -> list[list[str]]:
    """Monta as linhas da tabela de componentes de um filtro."""

    return [
        [
            result.name,
            "L",
            f"{format_pt(result.inductor.ideal * 1e3)} mH",
            f"{format_pt(result.inductor.commercial * 1e3)} mH",
            f"{format_pt(result.inductor.percent_error, 2)}%",
        ],
        [
            result.name,
            "C",
            f"{format_pt(result.capacitor.ideal * 1e6)} uF",
            f"{format_pt(result.capacitor.commercial * 1e6)} uF",
            f"{format_pt(result.capacitor.percent_error, 2)}%",
        ],
    ]


def analysis_text(project: ProjectResult) -> str:
    """Gera o texto usado na análise automática dos relatórios."""

    return (
        f"O indutor comercial difere "
        f"{format_pt(project.lpf.inductor.percent_error, 2)}% do valor ideal e "
        f"o capacitor difere "
        f"{format_pt(project.lpf.capacitor.percent_error, 2)}%. Com esses "
        f"componentes, o ponto de -3,01 dB foi estimado em "
        f"{format_pt(project.lpf.commercial_cutoff_hz, 1)} Hz para o LPF e "
        f"{format_pt(project.hpf.commercial_cutoff_hz, 1)} Hz para o HPF. "
        f"Considerando tolerância de ±{format_pt(project.tolerance_percent, 1)}% "
        f"em L e C, o corte pode variar de "
        f"{format_pt(project.lpf.tolerance_cutoff_min_hz, 1)} a "
        f"{format_pt(project.lpf.tolerance_cutoff_max_hz, 1)} Hz no LPF e de "
        f"{format_pt(project.hpf.tolerance_cutoff_min_hz, 1)} a "
        f"{format_pt(project.hpf.tolerance_cutoff_max_hz, 1)} Hz no HPF. "
        f"{project.audible_assessment}"
    )


def transfer_function_text(
    project: ProjectResult,
    commercial: bool = False,
) -> tuple[str, str]:
    """Retorna as funções de transferência com coeficientes numéricos."""

    result = project.lpf
    inductance = (
        result.inductor.commercial if commercial else result.inductor.ideal
    )
    capacitance = (
        result.capacitor.commercial if commercial else result.capacitor.ideal
    )
    second_order = inductance * capacitance
    first_order = inductance / project.load_ohm
    denominator = f"({second_order:.6e})s² + ({first_order:.6e})s + 1"
    return (
        f"H_LPF(s) = 1 / [{denominator}]",
        f"H_HPF(s) = ({second_order:.6e})s² / [{denominator}]",
    )


def generate_report_pdf(
    project: ProjectResult,
    pdf_path: str | Path,
) -> Path:
    """Gera somente o relatório PDF.

    O gráfico é criado temporariamente para ser incorporado ao documento e
    removido ao final, deixando apenas o arquivo PDF escolhido pelo usuário.
    """

    pdf = Path(pdf_path)
    if pdf.suffix.lower() != ".pdf":
        pdf = pdf.with_suffix(".pdf")
    pdf.parent.mkdir(parents=True, exist_ok=True)
    temporary_image = pdf.with_name(f".{pdf.stem}_bode_temporario.png")

    try:
        save_bode_figure(project, temporary_image)
        generate_pdf(project, pdf, temporary_image)
    finally:
        temporary_image.unlink(missing_ok=True)
    return pdf


def generate_pdf(project: ProjectResult, output: Path, image_path: Path) -> None:
    """Cria o relatório acadêmico em PDF."""

    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="BenaTitle",
            parent=styles["Title"],
            fontName="Helvetica-Bold",
            fontSize=24,
            textColor=colors.HexColor("#16324F"),
            alignment=TA_CENTER,
            spaceAfter=14,
        )
    )
    styles.add(
        ParagraphStyle(
            name="BenaHeading",
            parent=styles["Heading2"],
            fontName="Helvetica-Bold",
            textColor=colors.HexColor("#16324F"),
            spaceBefore=10,
            spaceAfter=7,
        )
    )
    styles.add(
        ParagraphStyle(
            name="BenaBody",
            parent=styles["BodyText"],
            alignment=TA_JUSTIFY,
            leading=15,
            spaceAfter=7,
        )
    )

    document = SimpleDocTemplate(
        str(output),
        pagesize=A4,
        rightMargin=1.8 * cm,
        leftMargin=1.8 * cm,
        topMargin=1.5 * cm,
        bottomMargin=1.5 * cm,
        title="BenaCCA - Relatório de crossover passivo",
    )

    story = [
        Spacer(1, 1.2 * cm),
        Paragraph("BenaCCA", styles["BenaTitle"]),
        Paragraph(
            "Projeto computacional de crossover passivo Butterworth de 2ª ordem",
            styles["Heading2"],
        ),
        Spacer(1, 0.5 * cm),
        Table(
            [
                ["Frequência de corte", f"{format_pt(project.cutoff_hz, 1)} Hz"],
                ["Impedância da carga", f"{format_pt(project.load_ohm, 2)} ohm"],
                [
                    "Tolerância considerada",
                    f"±{format_pt(project.tolerance_percent, 1)}%",
                ],
                ["Data de geração", datetime.now().strftime("%d/%m/%Y %H:%M")],
            ],
            colWidths=[5.0 * cm, 10.0 * cm],
            style=TableStyle(
                [
                    ("BACKGROUND", (0, 0), (0, -1), colors.HexColor("#E2E8F0")),
                    ("TEXTCOLOR", (0, 0), (0, -1), colors.HexColor("#16324F")),
                    ("FONTNAME", (0, 0), (0, -1), "Helvetica-Bold"),
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#94A3B8")),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("PADDING", (0, 0), (-1, -1), 7),
                ]
            ),
        ),
        Spacer(1, 0.6 * cm),
        Paragraph("1. Apresentação do problema", styles["BenaHeading"]),
        Paragraph(
            "O projeto dimensiona um crossover passivo de duas vias. O filtro "
            "passa-baixas envia as baixas frequências ao woofer, enquanto o filtro "
            "passa-altas envia as altas frequências ao tweeter.",
            styles["BenaBody"],
        ),
        Paragraph("2. Metodologia", styles["BenaHeading"]),
        Paragraph(
            "As duas seções utilizam a aproximação Butterworth de segunda ordem. "
            "No LPF, o indutor fica em série e o capacitor em paralelo com a carga. "
            "No HPF, o capacitor fica em série e o indutor em paralelo. As fórmulas "
            "são L = sqrt(2)R/wc e C = 1/(sqrt(2)Rwc), com wc = 2*pi*fc.",
            styles["BenaBody"],
        ),
        Paragraph("Funções de transferência ideais:", styles["BenaBody"]),
        Paragraph(
            "<br/>".join(transfer_function_text(project, commercial=False)),
            styles["BenaBody"],
        ),
        Paragraph("Funções com componentes comerciais:", styles["BenaBody"]),
        Paragraph(
            "<br/>".join(transfer_function_text(project, commercial=True)),
            styles["BenaBody"],
        ),
        Paragraph("3. Resultados dos componentes", styles["BenaHeading"]),
    ]

    component_data = [
        ["Filtro", "Elemento", "Ideal", "Comercial", "Erro"],
        *filter_rows(project.lpf),
        *filter_rows(project.hpf),
    ]
    story.append(
        Table(
            component_data,
            colWidths=[4.0 * cm, 2.0 * cm, 3.2 * cm, 3.2 * cm, 2.3 * cm],
            repeatRows=1,
            style=TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#16324F")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#F8FAFC")),
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#94A3B8")),
                    ("ALIGN", (1, 0), (-1, -1), "CENTER"),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("PADDING", (0, 0), (-1, -1), 6),
                ]
            ),
        )
    )
    story.extend(
        [
            Paragraph("4. Frequências e parâmetros reais", styles["BenaHeading"]),
            Table(
                [
                    ["Parâmetro", "LPF / Woofer", "HPF / Tweeter"],
                    [
                        "Corte real (-3,01 dB)",
                        f"{format_pt(project.lpf.commercial_cutoff_hz, 1)} Hz",
                        f"{format_pt(project.hpf.commercial_cutoff_hz, 1)} Hz",
                    ],
                    [
                        "Desvio do corte",
                        f"{format_pt(project.lpf.cutoff_error_percent, 2)}%",
                        f"{format_pt(project.hpf.cutoff_error_percent, 2)}%",
                    ],
                    [
                        "Diferença máx. da curva",
                        f"{format_pt(project.lpf.max_magnitude_difference_db, 2)} dB",
                        f"{format_pt(project.hpf.max_magnitude_difference_db, 2)} dB",
                    ],
                    [
                        "Q com valores comerciais",
                        format_pt(project.lpf.commercial_q, 4),
                        format_pt(project.hpf.commercial_q, 4),
                    ],
                    [
                        f"Faixa com tolerância ±{format_pt(project.tolerance_percent, 1)}%",
                        (
                            f"{format_pt(project.lpf.tolerance_cutoff_min_hz, 1)} a "
                            f"{format_pt(project.lpf.tolerance_cutoff_max_hz, 1)} Hz"
                        ),
                        (
                            f"{format_pt(project.hpf.tolerance_cutoff_min_hz, 1)} a "
                            f"{format_pt(project.hpf.tolerance_cutoff_max_hz, 1)} Hz"
                        ),
                    ],
                ],
                colWidths=[5.2 * cm, 5.0 * cm, 5.0 * cm],
                style=TableStyle(
                    [
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#16324F")),
                        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#94A3B8")),
                        ("ALIGN", (1, 1), (-1, -1), "CENTER"),
                        ("PADDING", (0, 0), (-1, -1), 6),
                    ]
                ),
            ),
            PageBreak(),
            Paragraph(
                "5. Gráfico de Bode comparativo — Magnitude e Fase",
                styles["BenaHeading"],
            ),
            Image(str(image_path), width=17.0 * cm, height=9.8 * cm),
            Paragraph("6. Análise crítica", styles["BenaHeading"]),
            Paragraph(analysis_text(project), styles["BenaBody"]),
            Paragraph("7. Conclusão", styles["BenaHeading"]),
            Paragraph(
                "O projeto atende ao objetivo de dimensionar os dois filtros, "
                "selecionar componentes disponíveis e quantificar o efeito da "
                "substituição. O principal desafio prático é que os valores "
                "comerciais não reproduzem exatamente o modelo teórico; por isso, "
                "a engenharia exige avaliar tolerâncias, resposta dos transdutores "
                "e medições do conjunto montado.",
                styles["BenaBody"],
            ),
        ]
    )
    document.build(story)


# =============================================================================
# INTERFACE GRÁFICA
# =============================================================================


class BenaCCAApp(tk.Tk):
    """Janela principal do software BenaCCA."""

    COLORS = {
        "header": "#0B0D12",
        "navy": "#DDE8FF",
        "blue": "#3B82F6",
        "orange": "#FDBA74",
        "bg": "#101114",
        "card": "#181A20",
        "card_soft": "#22252D",
        "input": "#111318",
        "text": "#F3F4F6",
        "muted": "#A7B0C0",
        "line": "#343946",
        "status": "#151821",
        "analysis_bg": "#2A1A10",
        "analysis_text": "#FFD7A6",
    }

    def __init__(self) -> None:
        super().__init__()
        self.title("BenaCCA — Projeto de Crossover Passivo")
        self.minsize(980, 720)
        # A interface inicia em tela cheia.
        self.fullscreen = True
        self.attributes("-fullscreen", True)
        self.configure(bg=self.COLORS["bg"])

        self.project: ProjectResult | None = None
        self.figure_canvas: FigureCanvasTkAgg | None = None
        self.graph_toolbar: BodeNavigationToolbar | None = None
        self.graph_motion_cid: int | None = None
        self.graph_leave_cid: int | None = None
        self.graph_draw_cid: int | None = None
        self.graph_resize_job: str | None = None
        self.graph_hover_artists: dict[object, tuple[object, list[object], object]] = {}
        self.graph_hover_backgrounds: dict[object, object] = {}
        self.graph_curve_data: dict[object, tuple[np.ndarray, list[object]]] = {}
        self.graph_active_axis: object | None = None
        self.graph_active_index: int | None = None

        self.frequency_var = tk.StringVar(value="2000")
        self.impedance_var = tk.StringVar(value="8")
        self.tolerance_var = tk.StringVar(value="5")
        self.status_var = tk.StringVar(
            value="Informe os parâmetros e clique em Calcular crossover."
        )

        self.configure_styles()
        self.build_layout()
        self.bind("<Return>", lambda _event: self.calculate())
        self.bind("<F11>", self.toggle_fullscreen)
        self.bind("<Escape>", self.leave_fullscreen)

    def toggle_fullscreen(self, _event: tk.Event | None = None) -> None:
        """Alterna entre tela cheia real e janela maximizada."""

        self.fullscreen = not self.fullscreen
        self.attributes("-fullscreen", self.fullscreen)
        if not self.fullscreen:
            self.state("zoomed")

    def leave_fullscreen(self, _event: tk.Event | None = None) -> None:
        """Sai da tela cheia e mantém o aplicativo maximizado."""

        self.fullscreen = False
        self.attributes("-fullscreen", False)
        self.state("zoomed")

    def configure_styles(self) -> None:
        """Define a aparência dos controles ttk."""

        style = ttk.Style(self)
        style.theme_use("clam")
        style.configure(
            ".",
            font=("Segoe UI", 10),
            background=self.COLORS["bg"],
            foreground=self.COLORS["text"],
        )
        style.configure("App.TFrame", background=self.COLORS["bg"])
        style.configure("Card.TFrame", background=self.COLORS["card"])
        style.configure("TSeparator", background=self.COLORS["line"])
        style.configure(
            "Title.TLabel",
            background=self.COLORS["header"],
            foreground="white",
            font=("Segoe UI Semibold", 22),
        )
        style.configure(
            "Subtitle.TLabel",
            background=self.COLORS["header"],
            foreground=self.COLORS["muted"],
            font=("Segoe UI", 10),
        )
        style.configure(
            "CardTitle.TLabel",
            background=self.COLORS["card"],
            foreground=self.COLORS["blue"],
            font=("Segoe UI Semibold", 13),
        )
        style.configure(
            "Body.TLabel",
            background=self.COLORS["card"],
            foreground=self.COLORS["text"],
        )
        style.configure(
            "Muted.TLabel",
            background=self.COLORS["card"],
            foreground=self.COLORS["muted"],
            font=("Segoe UI", 9),
        )
        style.configure(
            "Primary.TButton",
            background=self.COLORS["blue"],
            foreground="white",
            font=("Segoe UI Semibold", 10),
            padding=(16, 7),
            borderwidth=0,
        )
        style.map(
            "Primary.TButton",
            background=[("active", "#2563EB"), ("disabled", "#2E3542")],
            foreground=[
                ("disabled", self.COLORS["muted"]),
                ("active", "white"),
                ("!disabled", "white"),
            ],
        )
        style.configure(
            "Secondary.TButton",
            background=self.COLORS["card_soft"],
            foreground=self.COLORS["text"],
            padding=(13, 7),
            borderwidth=0,
        )
        style.map(
            "Secondary.TButton",
            background=[("active", "#2E3542"), ("disabled", self.COLORS["card"])],
            foreground=[
                ("disabled", self.COLORS["muted"]),
                ("!disabled", self.COLORS["text"]),
            ],
        )
        style.configure(
            "TEntry",
            fieldbackground=self.COLORS["input"],
            foreground=self.COLORS["text"],
            insertcolor=self.COLORS["text"],
            bordercolor=self.COLORS["line"],
            lightcolor=self.COLORS["line"],
            darkcolor=self.COLORS["line"],
            padding=(8, 5),
        )
        style.map(
            "TEntry",
            fieldbackground=[
                ("focus", "#151922"),
                ("disabled", self.COLORS["card_soft"]),
            ],
            foreground=[
                ("disabled", self.COLORS["muted"]),
                ("!disabled", self.COLORS["text"]),
            ],
            bordercolor=[("focus", self.COLORS["blue"])],
        )
        style.configure(
            "Vertical.TScrollbar",
            background=self.COLORS["card_soft"],
            troughcolor=self.COLORS["card"],
            bordercolor=self.COLORS["card"],
            arrowcolor=self.COLORS["muted"],
            lightcolor=self.COLORS["card"],
            darkcolor=self.COLORS["card"],
        )
        style.map(
            "Vertical.TScrollbar",
            background=[("active", self.COLORS["line"])],
            arrowcolor=[("active", self.COLORS["text"])],
        )
        # Remove o indicador de foco padrão e aplica o estilo das abas.
        style.layout(
            "Bena.TNotebook.Tab",
            [
                (
                    "Notebook.tab",
                    {
                        "sticky": "nswe",
                        "children": [
                            (
                                "Notebook.padding",
                                {
                                    "side": "top",
                                    "sticky": "nswe",
                                    "children": [
                                        (
                                            "Notebook.label",
                                            {"side": "top", "sticky": ""},
                                        )
                                    ],
                                },
                            )
                        ],
                    },
                )
            ],
        )
        style.configure(
            "Bena.TNotebook",
            background=self.COLORS["bg"],
            borderwidth=0,
            tabmargins=(0, 0, 0, 0),
        )
        style.configure(
            "Bena.TNotebook.Tab",
            background=self.COLORS["card_soft"],
            foreground=self.COLORS["muted"],
            borderwidth=0,
            padding=(24, 12),
            font=("Segoe UI Semibold", 10),
        )
        style.map(
            "Bena.TNotebook.Tab",
            background=[
                ("selected", self.COLORS["blue"]),
                ("active", "#2E3542"),
            ],
            foreground=[
                ("selected", "white"),
                ("active", self.COLORS["text"]),
            ],
            padding=[("selected", (26, 12))],
        )
        style.configure(
            "Treeview",
            rowheight=31,
            background=self.COLORS["input"],
            fieldbackground=self.COLORS["input"],
            foreground=self.COLORS["text"],
            bordercolor=self.COLORS["line"],
            lightcolor=self.COLORS["line"],
            darkcolor=self.COLORS["line"],
        )
        style.configure(
            "Treeview.Heading",
            background=self.COLORS["card_soft"],
            foreground="white",
            font=("Segoe UI Semibold", 9),
            padding=7,
            bordercolor=self.COLORS["line"],
        )
        style.map(
            "Treeview",
            background=[("selected", self.COLORS["blue"])],
            foreground=[("selected", "white")],
        )
        style.map("Treeview.Heading", background=[("active", "#2E3542")])

    def build_layout(self) -> None:
        """Monta o cabeçalho, o painel de entrada e a área de resultados."""

        header = tk.Frame(self, bg=self.COLORS["header"])
        header.pack(fill="x")
        header_inner = tk.Frame(header, bg=self.COLORS["header"])
        header_inner.pack(fill="x", padx=28, pady=(12, 11))
        title_box = tk.Frame(header_inner, bg=self.COLORS["header"])
        title_box.pack(side="left")
        ttk.Label(title_box, text="BenaCCA", style="Title.TLabel").pack(
            anchor="w"
        )
        ttk.Label(
            title_box,
            text="Dimensionamento e Análise de Crossover Passivo Butterworth",
            style="Subtitle.TLabel",
        ).pack(anchor="w")

        # Botão necessário enquanto a barra de título fica oculta.
        window_controls = tk.Frame(header_inner, bg=self.COLORS["header"])
        window_controls.pack(side="right", anchor="center")
        tk.Button(
            window_controls,
            text="Fechar",
            command=self.destroy,
            bg="#B42318",
            fg="white",
            activebackground="#D92D20",
            activeforeground="white",
            relief="flat",
            borderwidth=0,
            padx=18,
            pady=8,
            cursor="hand2",
            font=("Segoe UI Semibold", 9),
        ).pack(side="left")

        content = ttk.Frame(self, style="App.TFrame", padding=(18, 12, 18, 8))
        content.pack(fill="both", expand=True)
        content.columnconfigure(1, weight=1)
        content.rowconfigure(0, weight=1)

        self.build_input_panel(content)
        self.build_workspace(content)

        status = tk.Label(
            self,
            textvariable=self.status_var,
            bg=self.COLORS["status"],
            fg=self.COLORS["muted"],
            anchor="w",
            padx=22,
            pady=8,
            font=("Segoe UI", 9),
        )
        status.pack(fill="x", side="bottom")

    def build_input_panel(self, parent: ttk.Frame) -> None:
        """Cria os campos de entrada e os botões do programa."""

        panel = ttk.Frame(parent, style="Card.TFrame", padding=16)
        panel.grid(row=0, column=0, sticky="ns", padx=(0, 18))
        panel.configure(width=280)
        panel.grid_propagate(False)

        ttk.Label(panel, text="Parâmetros do projeto", style="CardTitle.TLabel").pack(
            anchor="w"
        )
        ttk.Label(
            panel,
            text="Use ponto ou vírgula nos valores decimais.",
            style="Muted.TLabel",
        ).pack(anchor="w", pady=(2, 12))

        self.labeled_entry(panel, "Frequência de corte (Hz)", self.frequency_var)
        self.labeled_entry(panel, "Impedância da carga (Ω)", self.impedance_var)
        self.labeled_entry(panel, "Tolerância dos componentes (%)", self.tolerance_var)

        ttk.Button(
            panel,
            text="Usar caso do enunciado: 2 kHz / 8 Ω",
            style="Secondary.TButton",
            command=self.load_standard_case,
        ).pack(fill="x", pady=(0, 10))

        ttk.Label(
            panel,
            text="Topologia",
            style="Body.TLabel",
            font=("Segoe UI Semibold", 9),
        ).pack(anchor="w", pady=(6, 3))
        topology = tk.Label(
            panel,
            text=(
                "LPF: L série + C paralelo\n"
                "HPF: C série + L paralelo\n"
                "Resposta: Butterworth, 2ª ordem"
            ),
            bg=self.COLORS["card_soft"],
            fg=self.COLORS["text"],
            justify="left",
            anchor="w",
            padx=8,
            pady=7,
            font=("Segoe UI", 9),
        )
        topology.pack(fill="x", pady=(0, 12))

        ttk.Button(
            panel,
            text="Calcular crossover",
            style="Primary.TButton",
            command=self.calculate,
        ).pack(fill="x", pady=(0, 8))
        ttk.Button(
            panel,
            text="Restaurar padrão",
            style="Secondary.TButton",
            command=self.reset,
        ).pack(fill="x")

        ttk.Separator(panel).pack(fill="x", pady=12)

        self.report_button = ttk.Button(
            panel,
            text="Gerar relatório PDF",
            style="Primary.TButton",
            command=self.generate_report,
            state="disabled",
        )
        self.report_button.pack(fill="x")

    def labeled_entry(
        self,
        parent: ttk.Frame,
        label: str,
        variable: tk.StringVar,
    ) -> None:
        """Adiciona um rótulo e um campo de entrada ao painel."""

        ttk.Label(
            parent,
            text=label,
            style="Body.TLabel",
            font=("Segoe UI Semibold", 9),
        ).pack(anchor="w", pady=(0, 4))
        entry = ttk.Entry(parent, textvariable=variable, font=("Segoe UI", 11))
        entry.pack(fill="x", ipady=4, pady=(0, 9))

    def build_workspace(self, parent: ttk.Frame) -> None:
        """Cria as abas de resultados e gráficos."""

        notebook = ttk.Notebook(parent, style="Bena.TNotebook", takefocus=False)
        notebook.grid(row=0, column=1, sticky="nsew")
        self.notebook = notebook

        self.results_tab = ttk.Frame(notebook, style="Card.TFrame", padding=18)
        self.graph_tab = ttk.Frame(notebook, style="Card.TFrame", padding=10)
        self.methodology_tab = ttk.Frame(notebook, style="Card.TFrame", padding=18)
        notebook.add(self.results_tab, text="Resultados")
        notebook.add(self.graph_tab, text="Bode: Magnitude e Fase")
        notebook.add(self.methodology_tab, text="Metodologia e Equações")
        notebook.bind("<<NotebookTabChanged>>", self.on_notebook_tab_changed)

        self.results_tab.columnconfigure(0, weight=1)
        # Mantém a tabela compacta mesmo quando há espaço vertical sobrando.
        self.results_tab.rowconfigure(4, weight=1)
        ttk.Label(
            self.results_tab,
            text="Resultados do dimensionamento",
            style="CardTitle.TLabel",
        ).grid(row=0, column=0, sticky="w", pady=(0, 12))

        columns = ("filter", "destination", "component", "ideal", "real", "error")
        self.results_tree = ttk.Treeview(
            self.results_tab,
            columns=columns,
            show="headings",
            height=8,
        )
        headings = {
            "filter": "Filtro",
            "destination": "Destino",
            "component": "Componente",
            "ideal": "Valor ideal",
            "real": "Valor comercial",
            "error": "Erro",
        }
        widths = {
            "filter": 120,
            "destination": 90,
            "component": 105,
            "ideal": 115,
            "real": 125,
            "error": 80,
        }
        for column in columns:
            self.results_tree.heading(column, text=headings[column])
            self.results_tree.column(
                column,
                width=widths[column],
                minwidth=70,
                anchor="center",
            )
        self.results_tree.grid(row=1, column=0, sticky="ew")

        summary_frame = ttk.Frame(
            self.results_tab,
            style="Card.TFrame",
            padding=(0, 16, 0, 0),
        )
        summary_frame.grid(row=2, column=0, sticky="ew")
        summary_frame.columnconfigure((0, 1, 2), weight=1)
        self.summary_labels: list[tk.Label] = []
        for column, title in enumerate(
            ("Corte comercial LPF", "Corte comercial HPF", "Fator Q comercial")
        ):
            card = tk.Frame(
                summary_frame,
                bg=self.COLORS["card_soft"],
                highlightbackground=self.COLORS["line"],
                highlightthickness=1,
            )
            card.grid(
                row=0,
                column=column,
                sticky="ew",
                padx=(0 if column == 0 else 6, 0 if column == 2 else 6),
            )
            tk.Label(
                card,
                text=title,
                bg=self.COLORS["card_soft"],
                fg=self.COLORS["muted"],
                font=("Segoe UI", 9),
            ).pack(anchor="w", padx=12, pady=(10, 2))
            value_label = tk.Label(
                card,
                text="—",
                bg=self.COLORS["card_soft"],
                fg=self.COLORS["orange"],
                font=("Segoe UI Semibold", 15),
            )
            value_label.pack(anchor="w", padx=12, pady=(0, 10))
            self.summary_labels.append(value_label)

        self.analysis_label = tk.Label(
            self.results_tab,
            text="Os resultados e a análise crítica aparecerão aqui.",
            bg=self.COLORS["analysis_bg"],
            fg=self.COLORS["analysis_text"],
            justify="left",
            anchor="nw",
            wraplength=760,
            padx=14,
            pady=12,
            font=("Segoe UI", 9),
        )
        self.analysis_label.grid(row=3, column=0, sticky="ew", pady=(14, 0))

        self.graph_tab.columnconfigure(0, weight=1)
        self.graph_tab.rowconfigure(2, weight=1)
        self.graph_hint = ttk.Label(
            self.graph_tab,
            text=(
                "Passe o mouse sobre uma curva para consultar os pontos. "
                "Use a barra abaixo para zoom, deslocamento e salvamento."
            ),
            style="Muted.TLabel",
        )
        self.graph_hint.grid(row=0, column=0, sticky="ew", pady=(0, 6))

        self.graph_toolbar_frame = ttk.Frame(
            self.graph_tab,
            style="Card.TFrame",
        )
        self.graph_toolbar_frame.grid(row=1, column=0, sticky="ew")

        self.graph_canvas_frame = ttk.Frame(
            self.graph_tab,
            style="Card.TFrame",
        )
        self.graph_canvas_frame.grid(row=2, column=0, sticky="nsew")
        self.graph_canvas_frame.columnconfigure(0, weight=1)
        self.graph_canvas_frame.rowconfigure(0, weight=1)
        self.graph_canvas_frame.bind("<Configure>", self.schedule_graph_resize)

        self.graph_placeholder = tk.Label(
            self.graph_canvas_frame,
            text="Calcule o crossover para gerar os gráficos.",
            bg=self.COLORS["card"],
            fg=self.COLORS["muted"],
            font=("Segoe UI", 12),
        )
        self.graph_placeholder.grid(row=0, column=0, sticky="nsew")
        self.build_methodology_tab()

    def build_methodology_tab(self) -> None:
        """Monta a documentação técnica exibida dentro do software."""

        container = ttk.Frame(self.methodology_tab, style="Card.TFrame")
        container.pack(fill="both", expand=True)
        container.columnconfigure(0, weight=1)
        container.rowconfigure(0, weight=1)

        scrollbar = ttk.Scrollbar(container, orient="vertical")
        scrollbar.grid(row=0, column=1, sticky="ns")

        self.methodology_text = tk.Text(
            container,
            wrap="word",
            yscrollcommand=scrollbar.set,
            bg=self.COLORS["input"],
            fg=self.COLORS["text"],
            insertbackground=self.COLORS["text"],
            selectbackground=self.COLORS["blue"],
            selectforeground="white",
            relief="flat",
            borderwidth=0,
            padx=18,
            pady=14,
            font=("Segoe UI", 10),
            spacing1=2,
            spacing3=7,
            cursor="arrow",
        )
        self.methodology_text.grid(row=0, column=0, sticky="nsew")
        scrollbar.configure(command=self.methodology_text.yview)

        self.methodology_text.tag_configure(
            "title",
            foreground=self.COLORS["text"],
            font=("Segoe UI Semibold", 16),
            spacing3=10,
        )
        self.methodology_text.tag_configure(
            "heading",
            foreground=self.COLORS["blue"],
            font=("Segoe UI Semibold", 12),
            spacing1=10,
            spacing3=5,
        )
        self.methodology_text.tag_configure(
            "equation",
            background=self.COLORS["card_soft"],
            foreground=self.COLORS["text"],
            font=("Consolas", 10),
            lmargin1=14,
            lmargin2=14,
            rmargin=14,
            spacing1=5,
            spacing3=5,
        )
        self.update_methodology()

    def update_methodology(self) -> None:
        """Atualiza fórmulas, funções numéricas e tabelas comerciais."""

        text = self.methodology_text
        text.configure(state="normal")
        text.delete("1.0", "end")

        text.insert("end", "Metodologia do projeto\n", "title")
        text.insert(
            "end",
            "O BenaCCA utiliza filtros passivos Butterworth de segunda ordem. "
            "No LPF, L fica em série e C em paralelo com o woofer. No HPF, "
            "C fica em série e L em paralelo com o tweeter.\n",
        )

        text.insert("end", "Equações de projeto\n", "heading")
        text.insert(
            "end",
            "ωc = 2πfc\n"
            "L = √2 · R / ωc\n"
            "C = 1 / (√2 · R · ωc)\n"
            "Q = R · √(C/L)\n",
            "equation",
        )

        text.insert("end", "Funções de transferência\n", "heading")
        text.insert(
            "end",
            "H_LPF(s) = R / (R + Ls + RLCs²)\n"
            "H_HPF(s) = RLCs² / (R + Ls + RLCs²)\n",
            "equation",
        )

        if self.project is None:
            text.insert(
                "end",
                "Calcule o crossover para visualizar aqui as funções de "
                "transferência com coeficientes numéricos.\n",
            )
        else:
            ideal_lpf, ideal_hpf = transfer_function_text(
                self.project, commercial=False
            )
            real_lpf, real_hpf = transfer_function_text(
                self.project, commercial=True
            )
            text.insert("end", "Coeficientes ideais\n", "heading")
            text.insert("end", f"{ideal_lpf}\n{ideal_hpf}\n", "equation")
            text.insert("end", "Coeficientes comerciais\n", "heading")
            text.insert("end", f"{real_lpf}\n{real_hpf}\n", "equation")

        text.insert("end", "Valores comerciais permitidos\n", "heading")
        inductors = ", ".join(
            f"{value * 1e3:g}" for value in COMMERCIAL_INDUCTORS_H
        )
        capacitors = ", ".join(
            f"{value * 1e6:g}" for value in COMMERCIAL_CAPACITORS_F
        )
        text.insert("end", f"Indutores (mH):\n{inductors}\n\n", "equation")
        text.insert("end", f"Capacitores (µF):\n{capacitors}\n", "equation")

        text.insert("end", "Análise de tolerância\n", "heading")
        text.insert(
            "end",
            "A faixa de corte é calculada combinando os quatro casos extremos: "
            "L mínimo/C mínimo, L mínimo/C máximo, L máximo/C mínimo e "
            "L máximo/C máximo. Essa análise aproxima o pior caso elétrico "
            "causado pela tolerância dos componentes.\n",
        )
        text.configure(state="disabled")

    def calculate(self) -> None:
        """Valida as entradas, calcula o projeto e atualiza a tela."""

        try:
            cutoff = parse_positive_number(
                self.frequency_var.get(),
                "a frequência de corte",
            )
            impedance = parse_positive_number(
                self.impedance_var.get(),
                "a impedância da carga",
            )
            tolerance = parse_positive_number(
                self.tolerance_var.get(),
                "a tolerância dos componentes",
            )
            if tolerance >= 100:
                raise ValueError("A tolerância deve ser menor que 100%.")
            self.project = calculate_project(
                cutoff,
                impedance,
                tolerance,
            )
        except ValueError as exc:
            messagebox.showerror("Dados inválidos", str(exc), parent=self)
            return

        self.update_results()
        self.update_graph()
        self.update_methodology()
        self.report_button.configure(state="normal")
        self.status_var.set(
            "Cálculo concluído. Gráficos e relatório PDF liberados."
        )

    def load_standard_case(self) -> None:
        """Carrega diretamente os parâmetros obrigatórios do enunciado."""

        self.frequency_var.set("2000")
        self.impedance_var.set("8")
        self.tolerance_var.set("5")
        self.status_var.set(
            "Caso do enunciado carregado: 2 kHz, 8 Ω e tolerância de 5%."
        )

    def update_results(self) -> None:
        """Preenche a tabela e os cartões de resumo."""

        if self.project is None:
            return
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)

        for result in (self.project.lpf, self.project.hpf):
            rows = (
                (
                    result.name,
                    result.destination,
                    "L (indutor)",
                    f"{result.inductor.ideal * 1e3:.3f} mH",
                    f"{result.inductor.commercial * 1e3:.3f} mH",
                    f"{result.inductor.percent_error:.2f}%",
                ),
                (
                    result.name,
                    result.destination,
                    "C (capacitor)",
                    f"{result.capacitor.ideal * 1e6:.3f} µF",
                    f"{result.capacitor.commercial * 1e6:.3f} µF",
                    f"{result.capacitor.percent_error:.2f}%",
                ),
            )
            for row in rows:
                self.results_tree.insert("", "end", values=row)

        self.summary_labels[0].configure(
            text=f"{self.project.lpf.commercial_cutoff_hz:.1f} Hz"
        )
        self.summary_labels[1].configure(
            text=f"{self.project.hpf.commercial_cutoff_hz:.1f} Hz"
        )
        self.summary_labels[2].configure(
            text=f"{self.project.lpf.commercial_q:.4f}"
        )

        analysis = (
            f"Análise automática: L varia "
            f"{self.project.lpf.inductor.percent_error:.2f}% e C varia "
            f"{self.project.lpf.capacitor.percent_error:.2f}%. "
            f"O corte muda {self.project.lpf.cutoff_error_percent:+.2f}% no LPF "
            f"e {self.project.hpf.cutoff_error_percent:+.2f}% no HPF. "
            f"Com tolerância de ±{self.project.tolerance_percent:.1f}%, a faixa "
            f"estimada é {self.project.lpf.tolerance_cutoff_min_hz:.1f}–"
            f"{self.project.lpf.tolerance_cutoff_max_hz:.1f} Hz no LPF e "
            f"{self.project.hpf.tolerance_cutoff_min_hz:.1f}–"
            f"{self.project.hpf.tolerance_cutoff_max_hz:.1f} Hz no HPF. "
            f"{self.project.audible_assessment}"
        )
        self.analysis_label.configure(text=analysis)

    def update_graph(self) -> None:
        """Exibe o gráfico de Bode dentro da segunda aba."""

        if self.project is None:
            return
        self.graph_placeholder.grid_remove()
        self.disconnect_graph_interaction()
        if self.figure_canvas is not None:
            self.figure_canvas.get_tk_widget().destroy()
        if self.graph_toolbar is not None:
            self.graph_toolbar.destroy()
            self.graph_toolbar = None

        # O canvas começa pequeno e depois acompanha o espaço disponível na aba.
        figure = build_bode_figure(self.project, figsize=(6.4, 4.8), theme="dark")
        self.figure_canvas = FigureCanvasTkAgg(
            figure,
            master=self.graph_canvas_frame,
        )
        canvas_widget = self.figure_canvas.get_tk_widget()
        canvas_widget.configure(
            width=1,
            height=1,
            highlightthickness=0,
            background=self.COLORS["card"],
        )
        canvas_widget.grid(row=0, column=0, sticky="nsew")

        self.graph_toolbar = BodeNavigationToolbar(
            self.figure_canvas,
            self.graph_toolbar_frame,
            pack_toolbar=False,
        )
        self.graph_toolbar.update()
        self.style_graph_toolbar()
        self.graph_toolbar.pack(fill="x")

        self.enable_graph_interaction()
        self.figure_canvas.draw()
        self.after_idle(self.fit_graph_to_container)

    def on_notebook_tab_changed(self, _event: tk.Event | None = None) -> None:
        """Reajusta o Bode quando sua aba se torna visível."""

        if self.notebook.select() == str(self.graph_tab):
            self.after_idle(self.fit_graph_to_container)

    def schedule_graph_resize(self, event: tk.Event) -> None:
        """Agrupa eventos de redimensionamento para evitar redesenhos excessivos."""

        if self.figure_canvas is None or event.width < 100 or event.height < 100:
            return
        if self.graph_resize_job is not None:
            self.after_cancel(self.graph_resize_job)
        self.graph_resize_job = self.after(80, self.fit_graph_to_container)

    def fit_graph_to_container(self) -> None:
        """Faz a figura ocupar apenas o espaço realmente disponível na aba."""

        self.graph_resize_job = None
        if self.figure_canvas is None:
            return
        width = self.graph_canvas_frame.winfo_width()
        height = self.graph_canvas_frame.winfo_height()
        if width < 100 or height < 100:
            return
        figure = self.figure_canvas.figure
        figure.set_size_inches(
            width / figure.dpi,
            height / figure.dpi,
            forward=False,
        )
        self.figure_canvas.draw_idle()

    def style_graph_toolbar(self) -> None:
        """Aplica o modo noturno aos controles nativos da barra do Matplotlib."""

        if self.graph_toolbar is None:
            return

        def apply_theme(widget: tk.Misc) -> None:
            widget_class = widget.winfo_class()
            options: dict[str, object] = {}
            if widget_class in {"Frame", "Labelframe"}:
                options = {"background": self.COLORS["card"]}
            elif widget_class in {"Button", "Checkbutton"}:
                options = {
                    "background": self.COLORS["card_soft"],
                    "foreground": self.COLORS["text"],
                    "activebackground": self.COLORS["line"],
                    "activeforeground": self.COLORS["text"],
                    "disabledforeground": self.COLORS["muted"],
                    "borderwidth": 0,
                    "highlightthickness": 0,
                    "relief": "flat",
                }
                if widget_class == "Checkbutton":
                    options["selectcolor"] = self.COLORS["card_soft"]
            elif widget_class == "Label":
                options = {
                    "background": self.COLORS["card"],
                    "foreground": self.COLORS["muted"],
                }
            elif widget_class == "Entry":
                options = {
                    "background": self.COLORS["input"],
                    "foreground": self.COLORS["text"],
                    "insertbackground": self.COLORS["text"],
                    "relief": "flat",
                }

            if options:
                try:
                    widget.configure(**options)
                except tk.TclError:
                    pass
            if (
                widget_class in {"Button", "Checkbutton"}
                and getattr(widget, "_image_file", None) is not None
            ):
                self.graph_toolbar._set_image_for_button(widget)
            for child in widget.winfo_children():
                apply_theme(child)

        apply_theme(self.graph_toolbar)

    @staticmethod
    def format_hover_frequency(frequency_hz: float) -> str:
        """Formata a frequência mostrada na leitura interativa."""

        if frequency_hz >= 1000:
            return f"{frequency_hz / 1000:.3f} kHz"
        return f"{frequency_hz:.2f} Hz"

    def enable_graph_interaction(self) -> None:
        """Adiciona leitura fluida dos pontos das curvas sob o cursor."""

        if self.figure_canvas is None:
            return

        figure = self.figure_canvas.figure
        self.graph_hover_artists.clear()
        self.graph_hover_backgrounds.clear()
        self.graph_curve_data.clear()
        self.graph_active_axis = None
        self.graph_active_index = None
        for axis in figure.axes:
            response_lines = [
                line
                for line in axis.get_lines()
                if (line.get_gid() or "").startswith("bode:")
            ]
            if response_lines:
                frequencies = np.asarray(response_lines[0].get_xdata(), dtype=float)
                self.graph_curve_data[axis] = (frequencies, response_lines)

            annotation = axis.annotate(
                "",
                xy=(0, 0),
                xytext=(14, 14),
                textcoords="offset points",
                color=self.COLORS["text"],
                bbox={
                    "boxstyle": "round,pad=0.45",
                    "facecolor": self.COLORS["card_soft"],
                    "edgecolor": self.COLORS["line"],
                    "alpha": 0.96,
                },
                arrowprops={
                    "arrowstyle": "->",
                    "color": self.COLORS["muted"],
                },
                fontsize=8,
                zorder=20,
            )
            annotation.set_visible(False)
            annotation.set_in_layout(False)
            annotation.set_animated(True)
            markers = [
                axis.plot(
                    [],
                    [],
                    marker="o",
                    markersize=5,
                    linestyle="none",
                    color=color,
                    markeredgecolor=self.COLORS["card"],
                    markeredgewidth=0.8,
                    zorder=19,
                )[0]
                for color in (self.COLORS["blue"], self.COLORS["orange"])
            ]
            for marker in markers:
                marker.set_animated(True)
            vertical_guide = axis.axvline(
                axis.get_xlim()[0],
                color=self.COLORS["muted"],
                linewidth=0.8,
                linestyle=":",
                alpha=0.55,
                visible=False,
                zorder=1,
            )
            vertical_guide.set_animated(True)
            self.graph_hover_artists[axis] = (
                annotation,
                markers,
                vertical_guide,
            )

        self.graph_motion_cid = self.figure_canvas.mpl_connect(
            "motion_notify_event",
            self.on_graph_mouse_move,
        )
        self.graph_leave_cid = self.figure_canvas.mpl_connect(
            "figure_leave_event",
            self.hide_graph_hover,
        )
        self.graph_draw_cid = self.figure_canvas.mpl_connect(
            "draw_event",
            self.capture_graph_backgrounds,
        )

    def capture_graph_backgrounds(self, _event: object | None = None) -> None:
        """Guarda o fundo pronto de cada eixo para atualizações por blitting."""

        if self.figure_canvas is None:
            return
        self.graph_hover_backgrounds = {
            axis: self.figure_canvas.copy_from_bbox(axis.bbox)
            for axis in self.graph_hover_artists
        }
        for annotation, markers, vertical_guide in self.graph_hover_artists.values():
            annotation.set_visible(False)
            for marker in markers:
                marker.set_visible(False)
            vertical_guide.set_visible(False)
        self.graph_active_axis = None
        self.graph_active_index = None

    def blit_graph_axis(self, axis: object) -> None:
        """Atualiza somente os indicadores móveis de um eixo."""

        if self.figure_canvas is None:
            return
        background = self.graph_hover_backgrounds.get(axis)
        if background is None:
            self.figure_canvas.draw_idle()
            return

        self.figure_canvas.restore_region(background)
        annotation, markers, vertical_guide = self.graph_hover_artists[axis]
        if vertical_guide.get_visible():
            axis.draw_artist(vertical_guide)
        for marker in markers:
            if marker.get_visible():
                axis.draw_artist(marker)
        if annotation.get_visible():
            axis.draw_artist(annotation)
        self.figure_canvas.blit(axis.bbox)

    def on_graph_mouse_move(self, event: object) -> None:
        """Atualiza o balão com os valores ideal e comercial mais próximos."""

        if (
            self.figure_canvas is None
            or getattr(event, "inaxes", None) not in self.graph_curve_data
            or getattr(event, "xdata", None) is None
        ):
            self.hide_graph_hover()
            return

        axis = event.inaxes
        frequencies, response_lines = self.graph_curve_data[axis]
        index = int(np.searchsorted(frequencies, event.xdata))
        index = max(0, min(index, frequencies.size - 1))
        if index > 0 and abs(frequencies[index - 1] - event.xdata) < abs(
            frequencies[index] - event.xdata
        ):
            index -= 1

        if axis is self.graph_active_axis and index == self.graph_active_index:
            return

        if self.graph_active_axis is not None and self.graph_active_axis is not axis:
            previous_axis = self.graph_active_axis
            previous_annotation, previous_markers, previous_guide = (
                self.graph_hover_artists[previous_axis]
            )
            previous_annotation.set_visible(False)
            for marker in previous_markers:
                marker.set_visible(False)
            previous_guide.set_visible(False)
            self.blit_graph_axis(previous_axis)

        frequency = float(frequencies[index])
        values = [
            float(np.asarray(line.get_ydata(), dtype=float)[index])
            for line in response_lines
        ]
        unit = "dB" if "Magnitude" in axis.get_ylabel() else "°"
        annotation, markers, vertical_guide = self.graph_hover_artists[axis]

        for other_axis, artists in self.graph_hover_artists.items():
            if other_axis is axis:
                continue
            artists[0].set_visible(False)
            for marker in artists[1]:
                marker.set_visible(False)
            artists[2].set_visible(False)

        for marker, value in zip(markers, values):
            marker.set_data([frequency], [value])
            marker.set_visible(True)
        vertical_guide.set_xdata([frequency, frequency])
        vertical_guide.set_visible(True)

        labels = [
            (line.get_gid() or "bode:Curva").split(":", 1)[1]
            for line in response_lines
        ]
        annotation.xy = (frequency, values[0])
        event_x = getattr(event, "x", axis.bbox.x0)
        event_y = getattr(event, "y", axis.bbox.y0)
        horizontal_offset = -14 if event_x > axis.bbox.x0 + axis.bbox.width / 2 else 14
        vertical_offset = -14 if event_y > axis.bbox.y0 + axis.bbox.height * 0.72 else 14
        annotation.xyann = (horizontal_offset, vertical_offset)
        annotation.set_horizontalalignment(
            "right" if horizontal_offset < 0 else "left"
        )
        annotation.set_verticalalignment(
            "top" if vertical_offset < 0 else "bottom"
        )
        annotation.set_text(
            "\n".join(
                [
                    f"f = {self.format_hover_frequency(frequency)}",
                    *[
                        f"{label}: {value:.2f} {unit}"
                        for label, value in zip(labels, values)
                    ],
                ]
            )
        )
        annotation.set_visible(True)
        self.graph_active_axis = axis
        self.graph_active_index = index
        self.blit_graph_axis(axis)

    def hide_graph_hover(self, _event: object | None = None) -> None:
        """Oculta os indicadores quando o cursor sai do gráfico."""

        if self.graph_active_axis is None:
            return
        axis = self.graph_active_axis
        annotation, markers, vertical_guide = self.graph_hover_artists[axis]
        annotation.set_visible(False)
        for marker in markers:
            marker.set_visible(False)
        vertical_guide.set_visible(False)
        self.graph_active_axis = None
        self.graph_active_index = None
        self.blit_graph_axis(axis)

    def disconnect_graph_interaction(self) -> None:
        """Desconecta os eventos associados à figura anterior."""

        if self.figure_canvas is not None:
            if self.graph_motion_cid is not None:
                self.figure_canvas.mpl_disconnect(self.graph_motion_cid)
            if self.graph_leave_cid is not None:
                self.figure_canvas.mpl_disconnect(self.graph_leave_cid)
            if self.graph_draw_cid is not None:
                self.figure_canvas.mpl_disconnect(self.graph_draw_cid)
        self.graph_motion_cid = None
        self.graph_leave_cid = None
        self.graph_draw_cid = None
        self.graph_hover_artists.clear()
        self.graph_hover_backgrounds.clear()
        self.graph_curve_data.clear()
        self.graph_active_axis = None
        self.graph_active_index = None

    def generate_report(self) -> None:
        """Solicita um destino e cria o relatório PDF."""

        if self.project is None:
            return
        selected = filedialog.asksaveasfilename(
            parent=self,
            title="Salvar relatório do BenaCCA",
            defaultextension=".pdf",
            filetypes=[("Documento PDF", "*.pdf")],
            initialfile="relatorio_benacca.pdf",
        )
        if not selected:
            return
        try:
            pdf = generate_report_pdf(self.project, selected)
        except (OSError, ValueError) as exc:
            messagebox.showerror("Erro ao gerar relatório", str(exc), parent=self)
            return
        self.status_var.set(f"Relatório gerado em: {pdf.parent}")
        messagebox.showinfo(
            "Relatório concluído",
            f"PDF gerado com sucesso:\n\n{pdf.name}",
            parent=self,
        )

    def reset(self) -> None:
        """Restaura os parâmetros padrão e limpa os resultados."""

        self.frequency_var.set("2000")
        self.impedance_var.set("8")
        self.tolerance_var.set("5")
        self.project = None
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        for label in self.summary_labels:
            label.configure(text="—")
        self.analysis_label.configure(
            text="Os resultados e a análise crítica aparecerão aqui."
        )
        self.disconnect_graph_interaction()
        if self.figure_canvas is not None:
            self.figure_canvas.get_tk_widget().destroy()
            self.figure_canvas = None
        if self.graph_toolbar is not None:
            self.graph_toolbar.destroy()
            self.graph_toolbar = None
        self.graph_placeholder.grid()
        self.update_methodology()
        self.report_button.configure(state="disabled")
        self.status_var.set(
            "Valores padrão restaurados: frequência de 2 kHz e carga de 8 Ω."
        )


# =============================================================================
# INICIALIZAÇÃO
# =============================================================================


def main() -> None:
    """Inicia o aplicativo."""

    # Mantém textos e gráficos nítidos quando o Windows usa escala de tela.
    if sys.platform == "win32":
        try:
            ctypes.windll.shcore.SetProcessDpiAwareness(2)
        except (AttributeError, OSError):
            try:
                ctypes.windll.user32.SetProcessDPIAware()
            except (AttributeError, OSError):
                pass

    app = BenaCCAApp()
    app.mainloop()


if __name__ == "__main__":
    main()
