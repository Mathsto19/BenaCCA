"""BenaCCA: projeto de crossover passivo Butterworth.

O programa reúne em um único arquivo:
- cálculos dos filtros LPF e HPF;
- seleção dos componentes comerciais;
- geração dos gráficos de Bode;
- exportação de relatório em PDF;
- interface gráfica em PyQt6.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime
from math import floor, isfinite, log10, pi, sqrt
from pathlib import Path
from typing import Literal

# Garante que o Matplotlib use o mesmo binding Qt da interface.
os.environ.setdefault("QT_API", "PyQt6")

import matplotlib

matplotlib.use("QtAgg")

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT,
)
from matplotlib.figure import Figure
from matplotlib.patches import Arc, Rectangle, Wedge
from PyQt6.QtCore import QEvent, Qt, QTimer
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)
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


class BodeNavigationToolbar(NavigationToolbar2QT):
    """Barra do gráfico sem o ajuste manual de subplots."""

    toolitems = tuple(
        item
        for item in NavigationToolbar2QT.toolitems
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

    # 1) Validação das entradas: tudo precisa ser finito e dentro dos limites.
    if not isfinite(cutoff_hz) or not isfinite(load_ohm):
        raise ValueError("Frequência e impedância devem ser números finitos.")
    if cutoff_hz <= 0 or load_ohm <= 0:
        raise ValueError("Frequência e impedância devem ser maiores que zero.")
    if not isfinite(tolerance_percent) or not 0 <= tolerance_percent < 100:
        raise ValueError("A tolerância deve estar entre 0% e 100%.")

    # 2) Valores ideais e seleção dos comerciais mais próximos.
    ideal_l, ideal_c = ideal_components(cutoff_hz, load_ohm)
    commercial_l = nearest_commercial(ideal_l, COMMERCIAL_INDUCTORS_H)
    commercial_c = nearest_commercial(ideal_c, COMMERCIAL_CAPACITORS_F)

    # 3) Parâmetros derivados dos componentes comerciais escolhidos.
    inductor = ComponentChoice("L", ideal_l, commercial_l)
    capacitor = ComponentChoice("C", ideal_c, commercial_c)
    natural_hz = 1.0 / (2.0 * pi * sqrt(commercial_l * commercial_c))
    q_factor = load_ohm * sqrt(commercial_c / commercial_l)

    # 4) Função interna que monta o resultado de um filtro (LPF ou HPF).
    def build_filter(kind: FilterKind, name: str, destination: str) -> FilterResult:
        # Faixa de corte no pior caso, variando L e C dentro da tolerância.
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

    # 5) Monta os dois filtros e devolve o resultado completo do projeto.
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

# Paletas de cores para os temas claro (relatório PDF) e escuro (interface).
# Ficam em nível de módulo porque são compartilhadas por todos os gráficos.
THEME_PALETTES: dict[str, dict[str, str]] = {
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


def style_axes(axis, palette: dict[str, str]) -> None:
    """Aplica cores de fundo, ticks e bordas a um eixo conforme o tema."""

    axis.set_facecolor(palette["axes"])
    axis.tick_params(axis="both", which="both", colors=palette["muted"])
    for spine in axis.spines.values():
        spine.set_color(palette["spine"])


def style_legend(legend, palette: dict[str, str]) -> None:
    """Pinta a moldura e os textos de uma legenda conforme o tema."""

    legend.get_frame().set_facecolor(palette["legend"])
    legend.get_frame().set_edgecolor(palette["legend_edge"])
    for label in legend.get_texts():
        label.set_color(palette["text"])


def build_bode_figure(
    project: ProjectResult,
    figsize: tuple[float, float] = (12.5, 7.2),
    theme: BodeTheme = "light",
) -> Figure:
    """Cria o Bode completo: magnitude e fase dos dois filtros."""

    # Eixo de frequência logarítmico, indo de 20x abaixo a 20x acima do corte.
    frequencies = np.logspace(
        np.log10(max(project.cutoff_hz / 20.0, 1.0)),
        np.log10(project.cutoff_hz * 20.0),
        1000,
    )

    palette = THEME_PALETTES[theme]

    # Cria a figura 2x2 (linha = filtro, coluna = magnitude/fase).
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
        style_axes(axis, palette)

    # Uma linha por filtro: calcula as respostas e desenha magnitude e fase.
    for row, result in enumerate((project.lpf, project.hpf)):
        # Respostas complexas com componentes ideais e comerciais.
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

        # Converte para magnitude (dB) e fase (graus) já desembrulhada.
        ideal_magnitude = magnitude_db(ideal_response)
        commercial_magnitude = magnitude_db(commercial_response)
        ideal_phase = np.degrees(np.unwrap(np.angle(ideal_response)))
        commercial_phase = np.degrees(np.unwrap(np.angle(commercial_response)))

        magnitude_axis = axes[row, 0]
        phase_axis = axes[row, 1]

        # --- Coluna da esquerda: gráfico de magnitude ---
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
        style_legend(magnitude_axis.legend(loc="best", fontsize=8), palette)

        # --- Coluna da direita: gráfico de fase ---
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
        style_legend(phase_axis.legend(loc="best", fontsize=8), palette)

    axes[1, 0].set_xlabel("Frequência (Hz)")
    axes[1, 1].set_xlabel("Frequência (Hz)")
    for axis in axes.flat:
        axis.xaxis.label.set_color(palette["text"])
        axis.yaxis.label.set_color(palette["text"])

    figure_title = figure.suptitle(
        f"Bode completo ideal × comercial | {project.cutoff_hz:g} Hz, "
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
# DESENHO DO ESQUEMÁTICO (símbolos elétricos)
# =============================================================================


def _draw_inductor_h(axis, x0, x1, y, color, bumps: int = 4) -> None:
    """Desenha um indutor (bobina) na horizontal, entre x0 e x1."""

    width = (x1 - x0) / bumps
    for index in range(bumps):
        center = x0 + width * (index + 0.5)
        axis.add_patch(
            Arc((center, y), width, width, theta1=0, theta2=180, color=color, lw=2.2)
        )


def _draw_inductor_v(axis, x, y0, y1, color, bumps: int = 4) -> None:
    """Desenha um indutor (bobina) na vertical, entre y0 e y1."""

    height = (y1 - y0) / bumps
    for index in range(bumps):
        center = y0 + height * (index + 0.5)
        axis.add_patch(
            Arc((x, center), height, height, theta1=-90, theta2=90, color=color, lw=2.2)
        )


def _draw_capacitor_h(axis, xc, y, color, gap: float = 0.34, plate: float = 0.9) -> None:
    """Desenha um capacitor (duas placas) na horizontal, centrado em xc."""

    for sign in (-1, 1):
        axis.plot(
            [xc + sign * gap / 2, xc + sign * gap / 2],
            [y - plate / 2, y + plate / 2],
            color=color,
            lw=2.4,
        )


def _draw_capacitor_v(axis, x, yc, color, gap: float = 0.34, plate: float = 0.9) -> None:
    """Desenha um capacitor (duas placas) na vertical, centrado em yc."""

    for sign in (-1, 1):
        axis.plot(
            [x - plate / 2, x + plate / 2],
            [yc + sign * gap / 2, yc + sign * gap / 2],
            color=color,
            lw=2.4,
        )


def _draw_load(axis, x, y0, y1, color, label: str) -> None:
    """Desenha a carga (alto-falante) como um retângulo rotulado."""

    axis.add_patch(
        Rectangle(
            (x - 0.45, y0 + 0.9),
            0.9,
            (y1 - y0) - 1.8,
            fill=False,
            edgecolor=color,
            lw=2.2,
        )
    )
    axis.text(
        x,
        (y0 + y1) / 2,
        label,
        ha="center",
        va="center",
        rotation=90,
        color=color,
        fontsize=8,
        fontweight="bold",
    )


def _draw_filter_schematic(axis, kind: FilterKind, result, palette) -> None:
    """Monta o esquemático completo de um filtro (LPF ou HPF)."""

    # 1) Cores e posições fixas dos trilhos e componentes.
    text = palette["text"]
    wire = palette["muted"]
    series_color = palette["ideal"]
    shunt_color = palette["commercial"]

    x_src, x_load = 1.0, 8.6
    y_bot, y_top = 1.0, 5.0
    x_shunt = 6.8

    # 2) Trilhos de entrada, terra e rótulos dos valores comerciais.
    axis.plot([x_src, x_src], [y_bot, y_top], color=wire, lw=1.8)
    axis.plot([x_src, x_load], [y_bot, y_bot], color=wire, lw=1.8)
    axis.text(x_src - 0.15, (y_bot + y_top) / 2, "Vin", rotation=90,
              ha="right", va="center", color=text, fontsize=8)

    inductance_label = f"L = {result.inductor.commercial * 1e3:.2f} mH"
    capacitance_label = f"C = {result.capacitor.commercial * 1e6:.2f} µF"

    # 3) Topologia: série e shunt mudam conforme LPF ou HPF.
    if kind == "LPF":
        # Indutor em série no trilho de cima.
        axis.plot([x_src, 3.0], [y_top, y_top], color=wire, lw=1.8)
        _draw_inductor_h(axis, 3.0, 5.4, y_top, series_color)
        axis.plot([5.4, x_load], [y_top, y_top], color=wire, lw=1.8)
        axis.text(4.2, y_top + 0.6, inductance_label, ha="center",
                  color=series_color, fontsize=8, fontweight="bold")
        # Capacitor em paralelo (shunt).
        axis.plot([x_shunt, x_shunt], [y_top, 3.45], color=wire, lw=1.8)
        _draw_capacitor_v(axis, x_shunt, 3.3, shunt_color)
        axis.plot([x_shunt, x_shunt], [3.15, y_bot], color=wire, lw=1.8)
        axis.text(x_shunt - 0.55, 3.3, capacitance_label, ha="right",
                  va="center", color=shunt_color, fontsize=8, fontweight="bold")
    else:
        # Capacitor em série no trilho de cima.
        axis.plot([x_src, 3.9], [y_top, y_top], color=wire, lw=1.8)
        _draw_capacitor_h(axis, 4.25, y_top, series_color)
        axis.plot([4.6, x_load], [y_top, y_top], color=wire, lw=1.8)
        axis.text(4.25, y_top + 0.6, capacitance_label, ha="center",
                  color=series_color, fontsize=8, fontweight="bold")
        # Indutor em paralelo (shunt).
        axis.plot([x_shunt, x_shunt], [y_top, 4.4], color=wire, lw=1.8)
        _draw_inductor_v(axis, x_shunt, 2.0, 4.4, shunt_color)
        axis.plot([x_shunt, x_shunt], [2.0, y_bot], color=wire, lw=1.8)
        axis.text(x_shunt - 0.55, 3.2, inductance_label, ha="right",
                  va="center", color=shunt_color, fontsize=8, fontweight="bold")

    # 4) Carga (alto-falante) e ajustes finais do eixo.
    _draw_load(axis, x_load, y_bot, y_top, text, f"{result.destination}\n{result.kind}")

    axis.set_xlim(0, 10)
    axis.set_ylim(0, 6.4)
    axis.set_aspect("equal")
    axis.axis("off")
    title = axis.set_title(f"{result.name} — {result.destination}", loc="left")
    title.set_color(text)


def _draw_selection_number_line(
    axis, ideal, commercial, available, scale, unit, palette
) -> None:
    """Desenha a reta com os valores comerciais e destaca o escolhido."""

    values = np.asarray(available) * scale
    axis.scatter(values, np.zeros_like(values), s=22,
                 color=palette["muted"], zorder=3)
    # Valor ideal (linha) e valor comercial escolhido (ponto destacado).
    axis.axvline(ideal * scale, color=palette["ideal"], lw=1.6, linestyle="--")
    axis.scatter([commercial * scale], [0], s=150,
                 color=palette["commercial"], zorder=5, edgecolor=palette["figure"])
    axis.annotate(
        f"ideal\n{ideal * scale:.3f}",
        xy=(ideal * scale, 0), xytext=(ideal * scale, 0.55),
        ha="center", color=palette["ideal"], fontsize=8,
    )
    axis.annotate(
        f"escolhido\n{commercial * scale:g}",
        xy=(commercial * scale, 0), xytext=(commercial * scale, -0.75),
        ha="center", color=palette["commercial"], fontsize=8, fontweight="bold",
    )
    axis.set_xscale("log")
    axis.set_ylim(-1.2, 1.2)
    axis.set_yticks([])
    axis.set_xlabel(unit, color=palette["text"])
    axis.tick_params(axis="x", colors=palette["muted"])
    for spine in ("left", "top", "right"):
        axis.spines[spine].set_visible(False)
    axis.spines["bottom"].set_color(palette["spine"])


def _draw_gauge(axis, cx, cy, radius, value, vmin, vmax, label, palette, target=None):
    """Desenha um medidor semicircular (estilo velocímetro)."""

    # 1) Arco de fundo e preenchimento proporcional ao valor.
    axis.add_patch(Wedge((cx, cy), radius, 0, 180, width=radius * 0.34,
                         facecolor=palette["axes"], edgecolor=palette["spine"]))
    fraction = min(max((value - vmin) / (vmax - vmin), 0.0), 1.0)
    angle = 180.0 * (1.0 - fraction)
    axis.add_patch(Wedge((cx, cy), radius, angle, 180, width=radius * 0.34,
                         facecolor=palette["ideal"]))

    # 2) Marcador opcional do valor de referência (ex.: Q ideal do Butterworth).
    if target is not None:
        t_fraction = min(max((target - vmin) / (vmax - vmin), 0.0), 1.0)
        t_angle = np.radians(180.0 * (1.0 - t_fraction))
        axis.plot(
            [cx + radius * 0.62 * np.cos(t_angle), cx + radius * 1.02 * np.cos(t_angle)],
            [cy + radius * 0.62 * np.sin(t_angle), cy + radius * 1.02 * np.sin(t_angle)],
            color=palette["muted"], lw=1.6, linestyle=":",
        )

    # 3) Ponteiro e rótulos numéricos.
    needle = np.radians(angle)
    axis.plot([cx, cx + radius * 0.86 * np.cos(needle)],
              [cy, cy + radius * 0.86 * np.sin(needle)],
              color=palette["commercial"], lw=2.4)
    axis.scatter([cx], [cy], s=30, color=palette["commercial"], zorder=5)
    axis.text(cx, cy - 0.22 * radius, f"{value:.3g}", ha="center", va="top",
              color=palette["text"], fontsize=12, fontweight="bold")
    axis.text(cx, cy - 0.5 * radius, label, ha="center", va="top",
              color=palette["muted"], fontsize=8)


def second_order_step(kind: FilterKind, natural_hz: float, q: float,
                      t: np.ndarray) -> np.ndarray:
    """Resposta ao degrau de um filtro de 2ª ordem (qualquer amortecimento).

    Usa a decomposição em polos para funcionar tanto no caso subamortecido
    (Butterworth, Q≈0,707) quanto em casos super/criticamente amortecidos.
    """

    # 1) Polos do sistema de 2ª ordem a partir de ωn e Q.
    omega_n = 2.0 * pi * natural_hz
    zeta = 1.0 / (2.0 * q)
    root = np.sqrt(complex(zeta**2 - 1.0))
    p1 = -zeta * omega_n + omega_n * root
    p2 = -zeta * omega_n - omega_n * root
    if abs(p1 - p2) < 1e-9 * omega_n:  # evita divisão por zero no caso crítico
        p2 = p1 * (1.0 + 1e-6)

    # 2) Resíduos e resposta no tempo conforme LPF ou HPF.
    if kind == "LPF":
        r1 = omega_n**2 / (p1 * (p1 - p2))
        r2 = omega_n**2 / (p2 * (p2 - p1))
        response = 1.0 + r1 * np.exp(p1 * t) + r2 * np.exp(p2 * t)
    else:
        r1 = p1 / (p1 - p2)
        r2 = p2 / (p2 - p1)
        response = r1 * np.exp(p1 * t) + r2 * np.exp(p2 * t)
    return np.real(response)


# =============================================================================
# GRÁFICOS EXTRAS (circuito, sistema e análise)
# =============================================================================


def build_schematic_figure(
    project: ProjectResult,
    figsize: tuple[float, float] = (11.0, 7.0),
    theme: BodeTheme = "light",
) -> Figure:
    """Esquemáticos dos dois filtros + reta de seleção dos componentes."""

    # 1) Figura 2×2 com o tema escolhido.
    palette = THEME_PALETTES[theme]
    figure = Figure(figsize=figsize, dpi=100, facecolor=palette["figure"])
    figure.subplots_adjust(left=0.06, right=0.96, top=0.9, bottom=0.1,
                           wspace=0.18, hspace=0.45)
    axes = figure.subplots(2, 2)
    for axis in axes.flat:
        axis.set_facecolor(palette["figure"])

    # 2) Linha de cima: esquemáticos do LPF e do HPF.
    _draw_filter_schematic(axes[0, 0], "LPF", project.lpf, palette)
    _draw_filter_schematic(axes[0, 1], "HPF", project.hpf, palette)

    # 3) Linha de baixo: reta de seleção do indutor e do capacitor.
    _draw_selection_number_line(
        axes[1, 0], project.lpf.inductor.ideal, project.lpf.inductor.commercial,
        COMMERCIAL_INDUCTORS_H, 1e3, "Indutores comerciais (mH)", palette,
    )
    axes[1, 0].set_title("Seleção do indutor", loc="left", color=palette["text"])
    _draw_selection_number_line(
        axes[1, 1], project.lpf.capacitor.ideal, project.lpf.capacitor.commercial,
        COMMERCIAL_CAPACITORS_F, 1e6, "Capacitores comerciais (µF)", palette,
    )
    axes[1, 1].set_title("Seleção do capacitor", loc="left", color=palette["text"])

    # 4) Título geral da figura.
    suptitle = figure.suptitle(
        "Topologia do crossover e escolha de componentes",
        fontsize=13, fontweight="bold",
    )
    suptitle.set_color(palette["text"])
    return figure


def build_system_figure(
    project: ProjectResult,
    figsize: tuple[float, float] = (11.0, 7.0),
    theme: BodeTheme = "light",
) -> Figure:
    """Resposta somada do sistema e divisão da faixa de áudio."""

    # 1) Figura com dois painéis: curvas em cima, faixa de áudio embaixo.
    palette = THEME_PALETTES[theme]
    figure = Figure(figsize=figsize, dpi=100, facecolor=palette["figure"])
    figure.subplots_adjust(left=0.08, right=0.97, top=0.9, bottom=0.1, hspace=0.45)
    top_axis, band_axis = figure.subplots(2, 1, height_ratios=[3, 1])

    # 2) Respostas individuais e soma das duas vias (com e sem inversão de fase).
    frequencies = np.logspace(
        np.log10(max(project.cutoff_hz / 30.0, 1.0)),
        np.log10(project.cutoff_hz * 30.0),
        1500,
    )
    lpf_response = transfer_response(
        "LPF", frequencies, project.load_ohm,
        project.lpf.inductor.commercial, project.lpf.capacitor.commercial,
    )
    hpf_response = transfer_response(
        "HPF", frequencies, project.load_ohm,
        project.hpf.inductor.commercial, project.hpf.capacitor.commercial,
    )
    # Soma elétrica das duas vias. Num crossover de 2ª ordem, as vias ficam
    # 180° defasadas no corte: em fase a soma cancela (notch); invertendo a
    # polaridade do tweeter o sistema soma corretamente.
    summed_in_phase = lpf_response + hpf_response
    summed_inverted = lpf_response - hpf_response

    # --- Gráfico de cima: vias individuais e soma do sistema ---
    style_axes(top_axis, palette)
    top_axis.semilogx(frequencies, magnitude_db(lpf_response),
                      color=palette["ideal"], lw=2.0, label="LPF (woofer)")
    top_axis.semilogx(frequencies, magnitude_db(hpf_response),
                      color=palette["commercial"], lw=2.0, label="HPF (tweeter)")
    top_axis.semilogx(frequencies, magnitude_db(summed_inverted),
                      color=palette["text"], lw=2.6,
                      label="Soma (tweeter invertido)")
    top_axis.semilogx(frequencies, magnitude_db(summed_in_phase),
                      color=palette["muted"], lw=1.6, linestyle="--",
                      label="Soma (mesma fase) → cancela")
    top_axis.axhline(0.0, color=palette["muted"], lw=1.0, linestyle=":")

    tolerance_cutoffs = [
        project.lpf.tolerance_cutoff_min_hz, project.lpf.tolerance_cutoff_max_hz,
        project.hpf.tolerance_cutoff_min_hz, project.hpf.tolerance_cutoff_max_hz,
    ]
    top_axis.axvspan(min(tolerance_cutoffs), max(tolerance_cutoffs),
                     color=palette["muted"], alpha=0.15,
                     label=f"faixa de corte (±{project.tolerance_percent:g}%)")
    top_axis.axvline(project.cutoff_hz, color=palette["muted"], lw=1.0, alpha=0.7)
    top_axis.set_ylim(-30, 6)
    top_axis.set_ylabel("Magnitude (dB)", color=palette["text"])
    top_axis.set_xlabel("Frequência (Hz)", color=palette["text"])
    top_axis.grid(True, which="both", color=palette["grid"], alpha=0.5)
    top_title = top_axis.set_title(
        "Resposta somada do sistema (woofer + tweeter)", loc="left")
    top_title.set_color(palette["text"])
    style_legend(top_axis.legend(loc="lower center", ncol=2, fontsize=8), palette)

    # --- Gráfico de baixo: qual faixa de frequência vai para cada alto-falante ---
    band_axis.set_xscale("log")
    band_axis.set_xlim(20, 20000)
    band_axis.set_ylim(0, 1)
    band_axis.axvspan(20, project.cutoff_hz, color=palette["ideal"], alpha=0.3)
    band_axis.axvspan(project.cutoff_hz, 20000, color=palette["commercial"], alpha=0.3)
    band_axis.axvspan(min(tolerance_cutoffs), max(tolerance_cutoffs),
                      color=palette["muted"], alpha=0.35)
    band_axis.axvline(project.cutoff_hz, color=palette["text"], lw=1.4)
    label_box = {"boxstyle": "round,pad=0.3", "facecolor": palette["axes"],
                 "edgecolor": palette["spine"], "alpha": 0.85}
    band_axis.text(np.sqrt(20 * project.cutoff_hz), 0.5, "WOOFER\n(graves)",
                   ha="center", va="center", color=palette["text"],
                   fontsize=9, fontweight="bold", bbox=label_box)
    band_axis.text(np.sqrt(project.cutoff_hz * 20000), 0.5, "TWEETER\n(agudos)",
                   ha="center", va="center", color=palette["text"],
                   fontsize=9, fontweight="bold", bbox=label_box)
    band_axis.text(project.cutoff_hz, 1.05, f"{project.cutoff_hz:g} Hz",
                   ha="center", va="bottom", color=palette["text"], fontsize=8)
    band_axis.set_yticks([])
    band_axis.set_xlabel("Frequência (Hz)", color=palette["text"])
    band_axis.tick_params(axis="x", colors=palette["muted"])
    for spine in band_axis.spines.values():
        spine.set_color(palette["spine"])
    band_title = band_axis.set_title("Divisão da faixa de áudio", loc="left")
    band_title.set_color(palette["text"])

    suptitle = figure.suptitle(
        f"Comportamento do sistema | {project.cutoff_hz:g} Hz, "
        f"{project.load_ohm:g} Ω", fontsize=13, fontweight="bold")
    suptitle.set_color(palette["text"])
    return figure


def build_analysis_figure(
    project: ProjectResult,
    figsize: tuple[float, float] = (11.0, 7.5),
    theme: BodeTheme = "light",
) -> Figure:
    """Polos e zeros, resposta ao degrau, barras de valores e medidores."""

    # 1) Grade 2×2 para os quatro painéis de análise.
    palette = THEME_PALETTES[theme]
    figure = Figure(figsize=figsize, dpi=100, facecolor=palette["figure"])
    figure.subplots_adjust(left=0.08, right=0.96, top=0.9, bottom=0.09,
                           wspace=0.28, hspace=0.4)
    axes = figure.subplots(2, 2)

    # --- Painel 1: polos e zeros no plano s ---
    pz_axis = axes[0, 0]
    style_axes(pz_axis, palette)
    omega_fc = 2.0 * pi * project.cutoff_hz
    for natural_hz, q, color, name in (
        (project.cutoff_hz, 1.0 / sqrt(2.0), palette["ideal"], "ideal"),
        (project.lpf.commercial_natural_hz, project.lpf.commercial_q,
         palette["commercial"], "comercial"),
    ):
        omega_n = 2.0 * pi * natural_hz
        zeta = 1.0 / (2.0 * q)
        root = np.sqrt(complex(zeta**2 - 1.0))
        poles = [-zeta * omega_n + omega_n * root, -zeta * omega_n - omega_n * root]
        pz_axis.scatter([p.real for p in poles], [p.imag for p in poles],
                        marker="x", s=90, color=color, label=f"polos ({name})")
    # Círculo de raio ωc: os polos do Butterworth ideal ficam sobre ele a 45°.
    angle = np.linspace(0, 2 * pi, 200)
    pz_axis.plot(omega_fc * np.cos(angle), omega_fc * np.sin(angle),
                 color=palette["muted"], lw=1.0, linestyle=":")
    pz_axis.scatter([0], [0], marker="o", s=90, facecolors="none",
                    edgecolors=palette["text"], label="zeros do HPF (s=0)")
    pz_axis.axhline(0, color=palette["spine"], lw=0.8)
    pz_axis.axvline(0, color=palette["spine"], lw=0.8)
    pz_axis.set_aspect("equal")
    pz_axis.set_xlabel("Re(s)", color=palette["text"])
    pz_axis.set_ylabel("Im(s)", color=palette["text"])
    pz_title = pz_axis.set_title("Polos e zeros (plano s)", loc="left")
    pz_title.set_color(palette["text"])
    style_legend(pz_axis.legend(loc="upper right", fontsize=7), palette)

    # --- Painel 2: resposta ao degrau (domínio do tempo) ---
    step_axis = axes[0, 1]
    style_axes(step_axis, palette)
    t = np.linspace(0, 5.0 / project.cutoff_hz, 700)
    for kind, color, name in (("LPF", palette["ideal"], "LPF (woofer)"),
                              ("HPF", palette["commercial"], "HPF (tweeter)")):
        natural = project.lpf.commercial_natural_hz
        q = project.lpf.commercial_q
        step_axis.plot(t * 1000.0, second_order_step(kind, natural, q, t),
                       color=color, lw=2.0, label=name)
    step_axis.axhline(0, color=palette["muted"], lw=0.8, linestyle=":")
    step_axis.set_xlabel("Tempo (ms)", color=palette["text"])
    step_axis.set_ylabel("Amplitude", color=palette["text"])
    step_axis.grid(True, color=palette["grid"], alpha=0.4)
    step_title = step_axis.set_title("Resposta ao degrau", loc="left")
    step_title.set_color(palette["text"])
    style_legend(step_axis.legend(loc="best", fontsize=7), palette)

    # --- Painel 3: barras comercial × ideal (ideal = 100%) ---
    bar_axis = axes[1, 0]
    style_axes(bar_axis, palette)
    l_ratio = 100.0 * project.lpf.inductor.commercial / project.lpf.inductor.ideal
    c_ratio = 100.0 * project.lpf.capacitor.commercial / project.lpf.capacitor.ideal
    positions = [0, 1]
    bars = bar_axis.bar(positions, [l_ratio, c_ratio], width=0.55,
                        color=[palette["ideal"], palette["commercial"]])
    bar_axis.axhline(100.0, color=palette["muted"], lw=1.2, linestyle="--")
    bar_axis.text(1.45, 100.0, "ideal", color=palette["muted"], fontsize=8,
                  va="center")
    labels = [
        f"{project.lpf.inductor.commercial * 1e3:.2f} mH",
        f"{project.lpf.capacitor.commercial * 1e6:.2f} µF",
    ]
    for bar, text in zip(bars, labels):
        bar_axis.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                      text, ha="center", color=palette["text"], fontsize=8)
    bar_axis.set_xticks(positions)
    bar_axis.set_xticklabels(["Indutor (L)", "Capacitor (C)"],
                             color=palette["text"])
    bar_axis.set_ylabel("% do valor ideal", color=palette["text"])
    bar_axis.set_ylim(0, 130)
    bar_title = bar_axis.set_title("Comercial × ideal", loc="left")
    bar_title.set_color(palette["text"])

    # --- Painel 4: medidores de Q e erro máximo de magnitude ---
    gauge_axis = axes[1, 1]
    gauge_axis.set_facecolor(palette["figure"])
    max_error = max(project.lpf.max_magnitude_difference_db,
                    project.hpf.max_magnitude_difference_db)
    _draw_gauge(gauge_axis, 0.27, 0.45, 0.24, project.lpf.commercial_q,
                0.3, 1.2, "Fator Q", palette, target=1.0 / sqrt(2.0))
    _draw_gauge(gauge_axis, 0.75, 0.45, 0.24, max_error, 0.0, 3.0,
                "Erro máx. (dB)", palette)
    gauge_axis.set_xlim(0, 1)
    gauge_axis.set_ylim(0, 0.85)
    gauge_axis.set_aspect("equal")
    gauge_axis.axis("off")
    gauge_title = gauge_axis.set_title("Indicadores", loc="left")
    gauge_title.set_color(palette["text"])

    suptitle = figure.suptitle(
        "Análise detalhada do projeto", fontsize=13, fontweight="bold")
    suptitle.set_color(palette["text"])
    return figure


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


def save_report_figure(
    project: ProjectResult,
    output_path: str | Path,
    builder,
    figsize: tuple[float, float],
) -> Path:
    """Salva uma das figuras do relatório em PNG temporário."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure = builder(project, figsize=figsize, theme="light")
    figure.savefig(
        output,
        dpi=180,
        bbox_inches="tight",
        facecolor=figure.get_facecolor(),
    )
    return output


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
    temporary_images = {
        "bode": pdf.with_name(f".{pdf.stem}_bode_temporario.png"),
        "circuit": pdf.with_name(f".{pdf.stem}_circuito_temporario.png"),
        "system": pdf.with_name(f".{pdf.stem}_sistema_temporario.png"),
        "analysis": pdf.with_name(f".{pdf.stem}_analise_temporario.png"),
    }

    try:
        save_bode_figure(project, temporary_images["bode"])
        save_report_figure(
            project,
            temporary_images["circuit"],
            build_schematic_figure,
            (11.0, 7.0),
        )
        save_report_figure(
            project,
            temporary_images["system"],
            build_system_figure,
            (11.0, 7.0),
        )
        save_report_figure(
            project,
            temporary_images["analysis"],
            build_analysis_figure,
            (11.0, 7.5),
        )
        generate_pdf(project, pdf, temporary_images)
    finally:
        for temporary_image in temporary_images.values():
            temporary_image.unlink(missing_ok=True)
    return pdf


def generate_pdf(
    project: ProjectResult,
    output: Path,
    image_paths: dict[str, Path],
) -> None:
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
    styles.add(
        ParagraphStyle(
            name="BenaSmall",
            parent=styles["BodyText"],
            fontSize=8.6,
            leading=11.3,
            spaceAfter=4,
        )
    )
    styles.add(
        ParagraphStyle(
            name="BenaCaption",
            parent=styles["Italic"],
            fontSize=8.5,
            leading=10.5,
            textColor=colors.HexColor("#475569"),
            alignment=TA_CENTER,
            spaceAfter=8,
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

    def body(text: str) -> Paragraph:
        return Paragraph(text, styles["BenaBody"])

    def small(text: str) -> Paragraph:
        return Paragraph(text, styles["BenaSmall"])

    def caption(text: str) -> Paragraph:
        return Paragraph(text, styles["BenaCaption"])

    def table_style(has_header: bool = True) -> TableStyle:
        commands = [
            ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#94A3B8")),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("PADDING", (0, 0), (-1, -1), 6),
        ]
        if has_header:
            commands.extend(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#16324F")),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ]
            )
        return TableStyle(commands)

    def trimmed_pt(value: float, decimals: int = 2) -> str:
        return format_pt(value, decimals).rstrip("0").rstrip(",")

    def commercial_values(values: np.ndarray, scale: float, unit: str) -> str:
        return ", ".join(
            f"{trimmed_pt(float(value) * scale)} {unit}" for value in values
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
        body(
            "O projeto dimensiona um crossover passivo de duas vias. O filtro "
            "passa-baixas envia as baixas frequências ao woofer, enquanto o filtro "
            "passa-altas envia as altas frequências ao tweeter. A transição foi "
            "modelada com filtros Butterworth de segunda ordem, mantendo a resposta "
            "suave ao redor da frequência de corte."
        ),
        Paragraph("2. Checklist do enunciado", styles["BenaHeading"]),
        Table(
            [
                ["Requisito", "Implementação no BenaCCA", "Status"],
                [
                    small("Receber frequência de corte e impedância da carga."),
                    small(
                        "Campos de entrada da interface e função "
                        "calculate_project(cutoff_hz, load_ohm, tolerance_percent)."
                    ),
                    small("Atendido"),
                ],
                [
                    small("Calcular valores ideais de L e C."),
                    small(
                        "Função ideal_components(), usando wc = 2*pi*fc, "
                        "L = sqrt(2)R/wc e C = 1/(sqrt(2)Rwc)."
                    ),
                    small("Atendido"),
                ],
                [
                    small("Sugerir componentes reais das tabelas comerciais."),
                    small(
                        "Listas COMMERCIAL_INDUCTORS_H e COMMERCIAL_CAPACITORS_F, "
                        "com escolha pelo menor erro absoluto."
                    ),
                    small("Atendido"),
                ],
                [
                    small("Gerar Bode comparativo ideal versus real."),
                    small(
                        "A função build_bode_figure() plota magnitude e fase para "
                        "componentes ideais e comerciais."
                    ),
                    small("Atendido"),
                ],
                [
                    small("Projetar LPF para woofer e HPF para tweeter."),
                    small(
                        "A função transfer_response() implementa as duas topologias "
                        "passivas de segunda ordem."
                    ),
                    small("Atendido"),
                ],
                [
                    small("Documentar lógica, resultados, análise crítica e conclusão."),
                    small(
                        "O repositório contém README, documentação acadêmica e este "
                        "PDF gerado com metodologia, resultados e figuras."
                    ),
                    small("Atendido"),
                ],
            ],
            colWidths=[4.6 * cm, 8.6 * cm, 2.1 * cm],
            repeatRows=1,
            style=table_style(),
        ),
        Paragraph("3. Bibliotecas e organização do código", styles["BenaHeading"]),
        body(
            "As bibliotecas ficam importadas no início do arquivo BenaCCA.py. "
            "Abaixo dos imports, o código é dividido por blocos comentados com "
            "linhas separadoras, como TABELAS DE COMPONENTES COMERCIAIS, CÁLCULOS "
            "DO CROSSOVER, GRÁFICOS DE BODE, GRÁFICOS EXTRAS, RELATÓRIO PDF e "
            "INTERFACE GRÁFICA."
        ),
        Table(
            [
                ["Biblioteca", "Uso no projeto"],
                [small("PyQt6"), small("Constrói a interface gráfica desktop.")],
                [
                    small("NumPy"),
                    small(
                        "Executa cálculos numéricos, vetores de frequência, "
                        "respostas complexas e busca do componente mais próximo."
                    ),
                ],
                [
                    small("Matplotlib"),
                    small(
                        "Gera os gráficos de Bode, circuito, resposta do sistema "
                        "e análise detalhada."
                    ),
                ],
                [
                    small("ReportLab"),
                    small("Monta e exporta este relatório técnico em PDF."),
                ],
                [
                    small("pathlib, tempfile, shutil, os"),
                    small(
                        "Organizam caminhos, arquivos temporários e tarefas de "
                        "apoio do sistema."
                    ),
                ],
            ],
            colWidths=[4.0 * cm, 11.2 * cm],
            repeatRows=1,
            style=table_style(),
        ),
        Paragraph("4. Metodologia", styles["BenaHeading"]),
        body(
            "As duas seções utilizam a aproximação Butterworth de segunda ordem. "
            "No LPF, o indutor fica em série e o capacitor em paralelo com a carga. "
            "No HPF, o capacitor fica em série e o indutor em paralelo. As fórmulas "
            "são L = sqrt(2)R/wc e C = 1/(sqrt(2)Rwc), com wc = 2*pi*fc."
        ),
        body("Funções de transferência ideais:"),
        body("<br/>".join(transfer_function_text(project, commercial=False))),
        body("Funções com componentes comerciais:"),
        body("<br/>".join(transfer_function_text(project, commercial=True))),
        Paragraph("5. Valores comerciais permitidos", styles["BenaHeading"]),
        body(
            "Indutores usados na seleção: "
            f"{commercial_values(COMMERCIAL_INDUCTORS_H, 1e3, 'mH')}."
        ),
        body(
            "Capacitores usados na seleção: "
            f"{commercial_values(COMMERCIAL_CAPACITORS_F, 1e6, 'uF')}."
        ),
        Paragraph("6. Resultados dos componentes", styles["BenaHeading"]),
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
            style=table_style(),
        )
    )
    story.extend(
        [
            Paragraph("7. Frequências e parâmetros reais", styles["BenaHeading"]),
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
                style=table_style(),
            ),
            PageBreak(),
            Paragraph(
                "8. Topologia e escolha dos componentes",
                styles["BenaHeading"],
            ),
            body(
                "A figura abaixo mostra os esquemáticos do LPF e do HPF, além da "
                "posição do valor ideal e do valor comercial escolhido nas tabelas "
                "permitidas pelo enunciado."
            ),
            Image(str(image_paths["circuit"]), width=17.0 * cm, height=10.8 * cm),
            caption("Figura 1 - Circuitos e seleção de L e C."),
            PageBreak(),
            Paragraph(
                "9. Gráfico de Bode comparativo - Magnitude e fase",
                styles["BenaHeading"],
            ),
            body(
                "As curvas azuis representam o projeto ideal e as curvas laranjas "
                "tracejadas representam o comportamento com os componentes "
                "comerciais escolhidos."
            ),
            Image(str(image_paths["bode"]), width=17.0 * cm, height=9.8 * cm),
            caption("Figura 2 - Bode comparativo ideal versus comercial."),
            Paragraph("10. Resposta do sistema", styles["BenaHeading"]),
            body(
                "A soma das duas vias mostra a divisão de energia entre woofer e "
                "tweeter. Como filtros de segunda ordem ficam defasados perto do "
                "corte, a figura também mostra a soma com o tweeter invertido."
            ),
            Image(str(image_paths["system"]), width=17.0 * cm, height=10.8 * cm),
            caption("Figura 3 - Resposta somada e faixa de atuação das vias."),
            PageBreak(),
            Paragraph(
                "11. Análise detalhada",
                styles["BenaHeading"],
            ),
            body(
                "A análise detalhada registra polos e zeros, resposta ao degrau, "
                "comparação percentual entre valores ideais e comerciais, fator Q "
                "e maior erro de magnitude."
            ),
            Image(str(image_paths["analysis"]), width=17.0 * cm, height=11.6 * cm),
            caption("Figura 4 - Polos, degrau, comparação de componentes e indicadores."),
            Paragraph("12. Análise crítica", styles["BenaHeading"]),
            body(analysis_text(project)),
            Paragraph("13. Conclusão", styles["BenaHeading"]),
            body(
                "O projeto atende ao objetivo de dimensionar os dois filtros, "
                "selecionar componentes disponíveis e quantificar o efeito da "
                "substituição. O principal desafio prático é que os valores "
                "comerciais não reproduzem exatamente o modelo teórico; por isso, "
                "a engenharia exige avaliar tolerâncias, resposta dos transdutores "
                "e medições do conjunto montado."
            ),
            body(
                "O modelo considera os alto-falantes como cargas resistivas de "
                f"{format_pt(project.load_ohm, 2)} ohm. Em uma caixa real, a "
                "impedância varia com a frequência; portanto, a validação final "
                "depende de medição acústica e teste de escuta."
            ),
        ]
    )
    document.build(story)


# =============================================================================
# INTERFACE GRÁFICA
# =============================================================================


class BenaCCAApp(QMainWindow):
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
    }

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("BenaCCA — Projeto de Crossover Passivo")
        self.setMinimumSize(980, 720)

        self._restore_fullscreen = True
        self._equation_dir = tempfile.mkdtemp(prefix="benacca_eq_")
        self.project: ProjectResult | None = None
        self.figure_canvas: FigureCanvasQTAgg | None = None
        self.graph_toolbar: BodeNavigationToolbar | None = None
        self.graph_motion_cid: int | None = None
        self.graph_leave_cid: int | None = None
        self.graph_draw_cid: int | None = None
        self.graph_hover_artists: dict[object, tuple[object, list[object], object]] = {}
        self.graph_hover_backgrounds: dict[object, object] = {}
        self.graph_curve_data: dict[object, tuple[np.ndarray, list[object]]] = {}
        self.graph_active_axis: object | None = None
        self.graph_active_index: int | None = None

        # Infraestrutura das abas visuais extras (circuito, sistema, análise).
        # Cada nome guarda o canvas atual, a barra e os layouts onde inseri-los.
        self.extra_builders = {
            "circuit": build_schematic_figure,
            "system": build_system_figure,
            "analysis": build_analysis_figure,
        }
        self.extra_canvases: dict[str, FigureCanvasQTAgg | None] = {}
        self.extra_toolbars: dict[str, NavigationToolbar2QT | None] = {}
        self.extra_canvas_layouts: dict[str, QVBoxLayout] = {}
        self.extra_toolbar_layouts: dict[str, QVBoxLayout] = {}
        self.extra_placeholders: dict[str, QLabel] = {}

        self.setStyleSheet(self.build_stylesheet())
        self.build_layout()

        # A interface inicia em tela cheia.
        self.showFullScreen()

    # ------------------------------------------------------------------
    # Aparência e tela cheia
    # ------------------------------------------------------------------
    def build_stylesheet(self) -> str:
        """Monta a folha de estilos (QSS) com o mesmo visual do tema escuro."""

        c = self.COLORS
        return f"""
        QWidget {{
            font-family: 'Segoe UI', 'DejaVu Sans', sans-serif;
            font-size: 10pt;
            color: {c['text']};
        }}
        #Central {{ background: {c['bg']}; }}

        #Header {{
            background: {c['header']};
            border-bottom: 2px solid {c['blue']};
        }}
        #TitleLabel {{ color: white; font-size: 22pt; font-weight: 600; }}
        #SubtitleLabel {{ color: {c['muted']}; font-size: 10pt; }}

        #MinimizeButton {{
            background: {c['card_soft']};
            color: {c['text']};
            border: 1px solid {c['line']};
            border-radius: 6px;
            min-width: 34px;
            max-width: 34px;
            min-height: 30px;
            max-height: 30px;
            font-size: 14pt;
            font-weight: 700;
        }}
        #MinimizeButton:hover {{ background: #2E3542; }}
        #RestartButton {{
            background: {c['card_soft']};
            color: {c['text']};
            border: 1px solid {c['line']};
            border-radius: 6px;
            min-width: 34px;
            max-width: 34px;
            min-height: 30px;
            max-height: 30px;
            font-size: 14pt;
            font-weight: 700;
        }}
        #RestartButton:hover {{ background: {c['blue']}; color: white; border-color: {c['blue']}; }}
        #CloseButton {{
            background: {c['card_soft']};
            color: {c['text']};
            border: 1px solid {c['line']};
            border-radius: 6px;
            min-width: 34px;
            max-width: 34px;
            min-height: 30px;
            max-height: 30px;
            font-size: 13pt;
            font-weight: 700;
        }}
        #CloseButton:hover {{ background: #D92D20; color: white; border-color: #D92D20; }}

        #Card {{ background: {c['card']}; }}
        #CardTitle {{ color: {c['blue']}; font-size: 13pt; font-weight: 600; }}
        #Muted {{ color: {c['muted']}; font-size: 9pt; }}
        #FieldLabel {{ color: {c['text']}; font-weight: 600; font-size: 9pt; }}
        #Topology {{
            background: {c['card_soft']};
            color: {c['text']};
            padding: 8px;
            font-size: 9pt;
        }}

        QLineEdit {{
            background: {c['input']};
            color: {c['text']};
            border: 1px solid {c['line']};
            border-radius: 2px;
            padding: 7px 8px;
            font-size: 11pt;
            selection-background-color: {c['blue']};
        }}
        QLineEdit:focus {{ border: 1px solid {c['blue']}; background: #151922; }}

        #PrimaryButton {{
            background: {c['blue']};
            color: white;
            border: none;
            padding: 9px 16px;
            font-weight: 600;
        }}
        #PrimaryButton:hover {{ background: #2563EB; }}
        #PrimaryButton:disabled {{ background: #2E3542; color: {c['muted']}; }}
        #SecondaryButton {{
            background: {c['card_soft']};
            color: {c['text']};
            border: none;
            padding: 8px 13px;
        }}
        #SecondaryButton:hover {{ background: #2E3542; }}
        #SecondaryButton:disabled {{ background: {c['card']}; color: {c['muted']}; }}

        #Separator {{ background: {c['line']}; max-height: 1px; min-height: 1px; }}

        QTabWidget::pane {{
            border: 1px solid {c['line']};
            border-top: 2px solid {c['blue']};
            background: {c['card']};
            top: -1px;
        }}
        QTabBar {{ background: {c['bg']}; }}
        QTabBar::tab {{
            background: {c['card_soft']};
            color: {c['muted']};
            padding: 12px 24px;
            margin-right: 3px;
            font-weight: 600;
            border: 1px solid {c['line']};
            border-bottom: none;
            border-top-left-radius: 6px;
            border-top-right-radius: 6px;
        }}
        QTabBar::tab:selected {{ background: {c['blue']}; color: white; border-color: {c['blue']}; }}
        QTabBar::tab:hover:!selected {{ background: #2E3542; color: {c['text']}; }}

        QTableWidget {{
            background: {c['input']};
            color: {c['text']};
            gridline-color: {c['line']};
            border: 1px solid {c['line']};
            selection-background-color: {c['blue']};
            selection-color: white;
        }}
        QTableWidget::item {{ padding: 6px; }}
        QHeaderView::section {{
            background: {c['card_soft']};
            color: white;
            font-weight: 600;
            font-size: 9pt;
            padding: 8px;
            border: none;
            border-right: 1px solid {c['line']};
        }}

        #SummaryCard {{ background: {c['card_soft']}; border: 1px solid {c['line']}; }}
        #SummaryTitle {{ color: {c['muted']}; font-size: 9pt; }}
        #SummaryValue {{ color: {c['orange']}; font-size: 15pt; font-weight: 600; }}

        QTextEdit {{
            background: {c['input']};
            color: {c['text']};
            border: none;
            selection-background-color: {c['blue']};
            selection-color: white;
        }}

        #Status {{
            background: {c['status']};
            color: {c['muted']};
            padding: 8px 22px;
            font-size: 9pt;
        }}

        QScrollBar:vertical {{
            background: {c['card']};
            width: 12px;
            margin: 0;
        }}
        QScrollBar::handle:vertical {{
            background: {c['card_soft']};
            min-height: 24px;
            border-radius: 4px;
        }}
        QScrollBar::handle:vertical:hover {{ background: {c['line']}; }}
        QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}
        QScrollBar::add-page:vertical, QScrollBar::sub-page:vertical {{
            background: {c['card']};
        }}
        """

    def keyPressEvent(self, event) -> None:  # noqa: N802 (assinatura do Qt)
        """Atalhos globais: F11 alterna tela cheia, Esc maximiza."""

        key = event.key()
        if key == Qt.Key.Key_F11:
            self.toggle_fullscreen()
        elif key == Qt.Key.Key_Escape:
            self.leave_fullscreen()
        elif key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            self.calculate()
        else:
            super().keyPressEvent(event)

    def toggle_fullscreen(self) -> None:
        """Alterna entre tela cheia real e janela maximizada."""

        if self.isFullScreen():
            self.showMaximized()
        else:
            self.showFullScreen()

    def leave_fullscreen(self) -> None:
        """Sai da tela cheia e mantém o aplicativo maximizado."""

        if self.isFullScreen():
            self.showMaximized()

    def minimize_window(self) -> None:
        """Minimiza a janela lembrando que ela deve voltar em tela cheia."""

        # Guarda o estado atual para restaurar exatamente como estava.
        self._restore_fullscreen = self.isFullScreen()
        self.showMinimized()

    def restart_program(self) -> None:
        """Reinicia o programa: limpa os valores e volta para a tela cheia."""

        self.reset()
        if self.notebook is not None:
            self.notebook.setCurrentIndex(0)
        self.showFullScreen()

    def changeEvent(self, event) -> None:  # noqa: N802 (assinatura do Qt)
        """Garante que a janela volte em tela cheia após ser minimizada."""

        if event.type() == QEvent.Type.WindowStateChange:
            restored = (
                not self.isMinimized()
                and not self.isFullScreen()
                and getattr(self, "_restore_fullscreen", False)
            )
            if restored:
                self._restore_fullscreen = False
                # O atraso evita conflito com a transição de estado do gerenciador.
                QTimer.singleShot(0, self.showFullScreen)
        super().changeEvent(event)

    def closeEvent(self, event) -> None:  # noqa: N802 (assinatura do Qt)
        """Remove as imagens temporárias das equações ao encerrar."""

        shutil.rmtree(self._equation_dir, ignore_errors=True)
        super().closeEvent(event)

    # ------------------------------------------------------------------
    # Construção do layout
    # ------------------------------------------------------------------
    def build_layout(self) -> None:
        """Monta o cabeçalho, o painel de entrada e a área de resultados."""

        central = QWidget()
        central.setObjectName("Central")
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        root.addWidget(self.build_header())

        content = QWidget()
        content.setObjectName("Central")
        content_layout = QHBoxLayout(content)
        content_layout.setContentsMargins(18, 12, 18, 8)
        content_layout.setSpacing(18)
        content_layout.addWidget(self.build_input_panel())
        content_layout.addWidget(self.build_workspace(), stretch=1)
        root.addWidget(content, stretch=1)

    def build_header(self) -> QWidget:
        """Cria o cabeçalho com título e os botões de janela."""

        header = QFrame()
        header.setObjectName("Header")
        layout = QHBoxLayout(header)
        layout.setContentsMargins(28, 12, 28, 11)

        title_box = QVBoxLayout()
        title_box.setSpacing(0)
        title = QLabel("BenaCCA")
        title.setObjectName("TitleLabel")
        subtitle = QLabel(
            "Dimensionamento e Análise de Crossover Passivo Butterworth"
        )
        subtitle.setObjectName("SubtitleLabel")
        title_box.addWidget(title)
        title_box.addWidget(subtitle)
        layout.addLayout(title_box)
        layout.addStretch(1)

        buttons = QHBoxLayout()
        buttons.setSpacing(8)

        minimize_button = QPushButton("\u2013")
        minimize_button.setObjectName("MinimizeButton")
        minimize_button.setToolTip("Minimizar")
        minimize_button.setCursor(Qt.CursorShape.PointingHandCursor)
        minimize_button.clicked.connect(self.minimize_window)
        buttons.addWidget(minimize_button)

        restart_button = QPushButton("\u21bb")
        restart_button.setObjectName("RestartButton")
        restart_button.setToolTip("Reiniciar o programa (limpa os valores)")
        restart_button.setCursor(Qt.CursorShape.PointingHandCursor)
        restart_button.clicked.connect(self.restart_program)
        buttons.addWidget(restart_button)

        close_button = QPushButton("\u2715")
        close_button.setObjectName("CloseButton")
        close_button.setToolTip("Fechar")
        close_button.setCursor(Qt.CursorShape.PointingHandCursor)
        close_button.clicked.connect(self.close)
        buttons.addWidget(close_button)

        layout.addLayout(buttons)
        return header

    def build_input_panel(self) -> QWidget:
        """Cria os campos de entrada e os botões do programa."""

        panel = QFrame()
        panel.setObjectName("Card")
        panel.setFixedWidth(280)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(0)

        title = QLabel("Parâmetros do projeto")
        title.setObjectName("CardTitle")
        layout.addWidget(title)

        hint = QLabel("Use ponto ou vírgula nos valores decimais.")
        hint.setObjectName("Muted")
        hint.setWordWrap(True)
        layout.addWidget(hint)
        layout.addSpacing(12)

        self.frequency_edit = self.labeled_entry(
            layout, "Frequência de corte (Hz)", "2000"
        )
        self.impedance_edit = self.labeled_entry(
            layout, "Impedância da carga (Ω)", "8"
        )
        self.tolerance_edit = self.labeled_entry(
            layout, "Tolerância dos componentes (%)", "5"
        )

        topology_title = QLabel("Topologia")
        topology_title.setObjectName("FieldLabel")
        layout.addWidget(topology_title)
        layout.addSpacing(3)
        topology = QLabel(
            "LPF: L série + C paralelo\n"
            "HPF: C série + L paralelo\n"
            "Resposta: Butterworth, 2ª ordem"
        )
        topology.setObjectName("Topology")
        layout.addWidget(topology)
        layout.addSpacing(12)

        calculate_button = QPushButton("Calcular crossover")
        calculate_button.setObjectName("PrimaryButton")
        calculate_button.setCursor(Qt.CursorShape.PointingHandCursor)
        calculate_button.clicked.connect(self.calculate)
        layout.addWidget(calculate_button)
        layout.addSpacing(8)

        reset_button = QPushButton("Restaurar padrão")
        reset_button.setObjectName("SecondaryButton")
        reset_button.setCursor(Qt.CursorShape.PointingHandCursor)
        reset_button.clicked.connect(self.reset)
        layout.addWidget(reset_button)

        separator = QFrame()
        separator.setObjectName("Separator")
        separator.setFrameShape(QFrame.Shape.HLine)
        layout.addSpacing(12)
        layout.addWidget(separator)
        layout.addSpacing(12)

        self.report_button = QPushButton("Gerar relatório PDF")
        self.report_button.setObjectName("PrimaryButton")
        self.report_button.setCursor(Qt.CursorShape.PointingHandCursor)
        self.report_button.setEnabled(False)
        self.report_button.clicked.connect(self.generate_report)
        layout.addWidget(self.report_button)

        layout.addStretch(1)
        return panel

    def labeled_entry(
        self,
        layout: QVBoxLayout,
        label: str,
        default: str,
    ) -> QLineEdit:
        """Adiciona um rótulo e um campo de entrada ao painel."""

        caption = QLabel(label)
        caption.setObjectName("FieldLabel")
        layout.addWidget(caption)
        layout.addSpacing(4)
        entry = QLineEdit(default)
        entry.returnPressed.connect(self.calculate)
        layout.addWidget(entry)
        layout.addSpacing(9)
        return entry

    def build_workspace(self) -> QWidget:
        """Cria as abas de resultados e gráficos."""

        notebook = QTabWidget()
        notebook.setObjectName("Bena")
        notebook.addTab(self.build_results_tab(), "Resultados")
        notebook.addTab(self.build_extra_tab("circuit"), "Circuito")
        notebook.addTab(self.build_graph_tab(), "Bode: Magnitude e Fase")
        notebook.addTab(self.build_extra_tab("system"), "Resposta do Sistema")
        notebook.addTab(self.build_extra_tab("analysis"), "Análise Detalhada")
        notebook.addTab(self.build_methodology_tab(), "Metodologia e Equações")
        self.notebook = notebook
        return notebook

    # ------------------------------------------------------------------
    # Abas visuais extras (circuito, sistema, análise)
    # ------------------------------------------------------------------
    def build_extra_tab(self, name: str) -> QWidget:
        """Cria uma aba de gráfico extra (circuito, sistema ou análise)."""

        tab = QWidget()
        tab.setObjectName("Card")
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)

        # Barra de ferramentas do Matplotlib (zoom, salvar, etc.).
        toolbar_holder = QWidget()
        toolbar_holder.setObjectName("Card")
        toolbar_layout = QVBoxLayout(toolbar_holder)
        toolbar_layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(toolbar_holder)

        # Área do gráfico; mostra um aviso enquanto não houver cálculo.
        canvas_holder = QWidget()
        canvas_holder.setObjectName("Card")
        canvas_layout = QVBoxLayout(canvas_holder)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(canvas_holder, stretch=1)

        placeholder = QLabel("Calcule o crossover para gerar esta visualização.")
        placeholder.setObjectName("Muted")
        placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        canvas_layout.addWidget(placeholder)

        # Guarda referências para inserir/remover o canvas posteriormente.
        self.extra_canvases[name] = None
        self.extra_toolbars[name] = None
        self.extra_canvas_layouts[name] = canvas_layout
        self.extra_toolbar_layouts[name] = toolbar_layout
        self.extra_placeholders[name] = placeholder
        return tab

    def build_results_tab(self) -> QWidget:
        """Monta a aba com a tabela e os cartões de resumo."""

        tab = QWidget()
        tab.setObjectName("Card")
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(18, 18, 18, 18)
        layout.setSpacing(12)

        title = QLabel("Resultados do dimensionamento")
        title.setObjectName("CardTitle")
        layout.addWidget(title)

        headings = [
            "Filtro",
            "Destino",
            "Componente",
            "Valor ideal",
            "Valor comercial",
            "Erro",
        ]
        self.results_table = QTableWidget(0, len(headings))
        self.results_table.setHorizontalHeaderLabels(headings)
        self.results_table.verticalHeader().setVisible(False)
        self.results_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.results_table.setSelectionMode(
            QAbstractItemView.SelectionMode.NoSelection
        )
        self.results_table.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.results_table.setShowGrid(True)
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.results_table.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.results_table.setMaximumHeight(330)
        layout.addWidget(self.results_table)

        summary = QWidget()
        summary.setObjectName("Card")
        summary_layout = QHBoxLayout(summary)
        summary_layout.setContentsMargins(0, 4, 0, 0)
        summary_layout.setSpacing(12)
        self.summary_value_labels: list[QLabel] = []
        for caption in (
            "Corte comercial LPF",
            "Corte comercial HPF",
            "Fator Q comercial",
        ):
            card = QFrame()
            card.setObjectName("SummaryCard")
            card_layout = QVBoxLayout(card)
            card_layout.setContentsMargins(12, 10, 12, 10)
            card_layout.setSpacing(2)
            card_title = QLabel(caption)
            card_title.setObjectName("SummaryTitle")
            value = QLabel("—")
            value.setObjectName("SummaryValue")
            card_layout.addWidget(card_title)
            card_layout.addWidget(value)
            summary_layout.addWidget(card)
            self.summary_value_labels.append(value)
        layout.addWidget(summary)
        layout.addStretch(1)
        return tab

    def build_graph_tab(self) -> QWidget:
        """Monta a aba do gráfico de Bode."""

        tab = QWidget()
        tab.setObjectName("Card")
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(6)

        self.toolbar_holder = QWidget()
        self.toolbar_holder.setObjectName("Card")
        self.toolbar_layout = QVBoxLayout(self.toolbar_holder)
        self.toolbar_layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toolbar_holder)

        self.canvas_holder = QWidget()
        self.canvas_holder.setObjectName("Card")
        self.canvas_layout = QVBoxLayout(self.canvas_holder)
        self.canvas_layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.canvas_holder, stretch=1)

        self.graph_placeholder = QLabel(
            "Calcule o crossover para gerar os gráficos."
        )
        self.graph_placeholder.setObjectName("Muted")
        self.graph_placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.canvas_layout.addWidget(self.graph_placeholder)
        return tab

    def build_methodology_tab(self) -> QWidget:
        """Monta a documentação técnica exibida dentro do software."""

        tab = QWidget()
        tab.setObjectName("Card")
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(18, 18, 18, 18)

        self.methodology_text = QTextEdit()
        self.methodology_text.setReadOnly(True)
        layout.addWidget(self.methodology_text)
        self.update_methodology()
        return tab

    # ------------------------------------------------------------------
    # Conteúdo dinâmico
    # ------------------------------------------------------------------
    def render_equation(
        self,
        key: str,
        latex: str,
        fontsize: int = 19,
        cache: bool = True,
    ) -> str:
        """Renderiza uma equação LaTeX em PNG e devolve a URL do arquivo.

        As equações fixas são guardadas em cache; as numéricas (que dependem do
        projeto) são sempre regeradas para refletir os valores atuais.
        """

        path = Path(self._equation_dir) / f"{key}.png"
        if not (cache and path.exists()):
            try:
                self._save_equation_png(path, latex, fontsize)
            except (ValueError, RuntimeError):
                # Caso o mathtext não consiga interpretar, mostra texto simples.
                fallback = latex.replace("$", "").replace("\\dfrac", "")
                self._save_equation_png(path, fallback, fontsize, math=False)
        return path.as_uri()

    def _save_equation_png(
        self,
        path: Path,
        text: str,
        fontsize: int,
        math: bool = True,
    ) -> None:
        """Desenha o texto (matemático ou simples) em um PNG transparente."""

        figure = Figure(figsize=(0.1, 0.1), dpi=150)
        figure.patch.set_alpha(0.0)
        FigureCanvasAgg(figure)
        figure.text(
            0.0,
            0.0,
            text,
            fontsize=fontsize,
            color=self.COLORS["text"],
            usetex=False,
            parse_math=math,
        )
        figure.savefig(
            path,
            transparent=True,
            bbox_inches="tight",
            pad_inches=0.06,
        )

    @staticmethod
    def latex_scientific(value: float) -> str:
        """Formata um número em notação científica para o LaTeX do mathtext."""

        if value == 0:
            return "0"
        exponent = floor(log10(abs(value)))
        mantissa = value / 10**exponent
        mantissa_text = f"{mantissa:.2f}".replace(".", "{,}")
        return rf"{mantissa_text}\times10^{{{exponent}}}"

    def coefficient_equations(self, commercial: bool) -> tuple[str, str]:
        """Monta o LaTeX das funções de transferência com números reais."""

        result = self.project.lpf
        inductance = (
            result.inductor.commercial if commercial else result.inductor.ideal
        )
        capacitance = (
            result.capacitor.commercial if commercial else result.capacitor.ideal
        )
        second_order = self.latex_scientific(inductance * capacitance)
        first_order = self.latex_scientific(inductance / self.project.load_ohm)
        denominator = rf"{second_order}\,s^{{2}} + {first_order}\,s + 1"
        return (
            rf"$H_{{LPF}}(s) = \dfrac{{1}}{{{denominator}}}$",
            rf"$H_{{HPF}}(s) = \dfrac{{{second_order}\,s^{{2}}}}{{{denominator}}}$",
        )

    def update_methodology(self) -> None:
        """Atualiza fórmulas, funções numéricas e tabelas comerciais."""

        c = self.COLORS

        def heading(text: str) -> str:
            return (
                f'<p style="color:{c["blue"]};font-size:13pt;font-weight:600;'
                f'margin:20px 0 6px 0;">{text}</p>'
            )

        def paragraph(text: str) -> str:
            return (
                f'<p style="color:{c["text"]};line-height:150%;'
                f'margin:0 0 8px 0;">{text}</p>'
            )

        def equation_card(
            keys_and_latex: list[tuple[str, str]],
            cache: bool = True,
        ) -> str:
            rows = ""
            for key, latex in keys_and_latex:
                uri = self.render_equation(key, latex, cache=cache)
                rows += (
                    f'<p align="center" style="margin:10px 0;">'
                    f'<img src="{uri}"></p>'
                )
            return (
                f'<table width="100%" cellspacing="0" cellpadding="0">'
                f'<tr><td style="background:{c["card_soft"]};'
                f'padding:8px 14px;">{rows}</td></tr></table>'
            )

        def value_grid(values, scale: float, unit: str, columns: int = 9) -> str:
            cells = [f"{value * scale:g}" for value in values]
            cell_style = (
                f"background:{c['input']};color:{c['text']};"
                f"border:1px solid {c['line']};padding:6px 4px;"
                f"font-size:9pt;"
            )
            rows = ""
            for start in range(0, len(cells), columns):
                chunk = cells[start : start + columns]
                chunk += [""] * (columns - len(chunk))
                line = "".join(
                    f'<td align="center" style="{cell_style}">{item}</td>'
                    for item in chunk
                )
                rows += f"<tr>{line}</tr>"
            caption = (
                f'<p style="color:{c["muted"]};font-size:9pt;'
                f'margin:4px 0;">{unit}</p>'
            )
            return (
                caption
                + f'<table width="100%" cellspacing="0" cellpadding="0">'
                + rows
                + "</table>"
                # Qt não cria espaço após tabelas, então um parágrafo vazio
                # garante a separação visual para a próxima seção.
                + '<p style="margin:0;">&nbsp;</p>'
            )

        parts = [
            f'<p style="color:{c["text"]};font-size:17pt;font-weight:600;'
            f'margin-bottom:10px;">Metodologia do projeto</p>',
            paragraph(
                "O BenaCCA utiliza filtros passivos Butterworth de segunda ordem. "
                "No LPF, L fica em série e C em paralelo com o woofer. No HPF, "
                "C fica em série e L em paralelo com o tweeter."
            ),
            heading("Equações de projeto"),
            equation_card(
                [
                    ("eq_omega", r"$\omega_c = 2\pi f_c$"),
                    ("eq_l", r"$L = \dfrac{\sqrt{2}\,R}{\omega_c}$"),
                    ("eq_c", r"$C = \dfrac{1}{\sqrt{2}\,R\,\omega_c}$"),
                    ("eq_q", r"$Q = R\,\sqrt{\dfrac{C}{L}}$"),
                ]
            ),
            heading("Funções de transferência"),
            equation_card(
                [
                    (
                        "eq_hlpf",
                        r"$H_{LPF}(s) = \dfrac{R}{R + Ls + RLCs^{2}}$",
                    ),
                    (
                        "eq_hhpf",
                        r"$H_{HPF}(s) = \dfrac{RLCs^{2}}{R + Ls + RLCs^{2}}$",
                    ),
                ]
            ),
        ]

        if self.project is None:
            parts.append(
                paragraph(
                    "Calcule o crossover para visualizar aqui as funções de "
                    "transferência com os coeficientes numéricos do seu projeto."
                )
            )
        else:
            ideal_lpf, ideal_hpf = self.coefficient_equations(commercial=False)
            real_lpf, real_hpf = self.coefficient_equations(commercial=True)
            parts.append(heading("Coeficientes ideais"))
            parts.append(
                equation_card(
                    [
                        ("eq_ideal_lpf", ideal_lpf),
                        ("eq_ideal_hpf", ideal_hpf),
                    ],
                    cache=False,
                )
            )
            parts.append(heading("Coeficientes comerciais"))
            parts.append(
                equation_card(
                    [
                        ("eq_real_lpf", real_lpf),
                        ("eq_real_hpf", real_hpf),
                    ],
                    cache=False,
                )
            )

        parts.append(heading("Valores comerciais permitidos"))
        parts.append(value_grid(COMMERCIAL_INDUCTORS_H, 1e3, "Indutores (mH)"))
        parts.append(value_grid(COMMERCIAL_CAPACITORS_F, 1e6, "Capacitores (µF)"))

        parts.append(heading("Análise de tolerância"))
        parts.append(
            paragraph(
                "A faixa de corte é calculada combinando os quatro casos extremos:"
            )
        )
        parts.append(
            paragraph(
                "&bull; L mínimo / C mínimo<br/>"
                "&bull; L mínimo / C máximo<br/>"
                "&bull; L máximo / C mínimo<br/>"
                "&bull; L máximo / C máximo"
            )
        )
        parts.append(
            paragraph(
                "Essa análise aproxima o pior caso elétrico causado pela "
                "tolerância dos componentes."
            )
        )
        self.methodology_text.setHtml("".join(parts))

    def calculate(self) -> None:
        """Valida as entradas, calcula o projeto e atualiza a tela."""

        # 1) Lê e valida os campos; em caso de erro, avisa e interrompe.
        try:
            cutoff = parse_positive_number(
                self.frequency_edit.text(),
                "a frequência de corte",
            )
            impedance = parse_positive_number(
                self.impedance_edit.text(),
                "a impedância da carga",
            )
            tolerance = parse_positive_number(
                self.tolerance_edit.text(),
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
            QMessageBox.critical(self, "Dados inválidos", str(exc))
            return

        # 2) Atualiza todas as áreas da interface com o novo projeto.
        self.update_results()
        self.update_graph()
        self.update_extra_views()  # circuito, sistema e análise detalhada
        self.update_methodology()
        self.report_button.setEnabled(True)

    def update_results(self) -> None:
        """Preenche a tabela e os cartões de resumo."""

        if self.project is None:
            return
        self.results_table.setRowCount(0)

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
                position = self.results_table.rowCount()
                self.results_table.insertRow(position)
                for column, value in enumerate(row):
                    item = QTableWidgetItem(value)
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.results_table.setItem(position, column, item)

        self.summary_value_labels[0].setText(
            f"{self.project.lpf.commercial_cutoff_hz:.1f} Hz"
        )
        self.summary_value_labels[1].setText(
            f"{self.project.hpf.commercial_cutoff_hz:.1f} Hz"
        )
        self.summary_value_labels[2].setText(
            f"{self.project.lpf.commercial_q:.4f}"
        )

    def update_graph(self) -> None:
        """Exibe o gráfico de Bode dentro da segunda aba."""

        if self.project is None:
            return

        # 1) Remove o gráfico e a barra de ferramentas anteriores, se existirem.
        self.graph_placeholder.hide()
        self.disconnect_graph_interaction()
        if self.figure_canvas is not None:
            self.canvas_layout.removeWidget(self.figure_canvas)
            self.figure_canvas.setParent(None)
            self.figure_canvas.deleteLater()
            self.figure_canvas = None
        if self.graph_toolbar is not None:
            self.toolbar_layout.removeWidget(self.graph_toolbar)
            self.graph_toolbar.setParent(None)
            self.graph_toolbar.deleteLater()
            self.graph_toolbar = None

        # 2) Cria a nova figura (tema escuro) e a coloca no canvas da aba.
        figure = build_bode_figure(self.project, figsize=(6.4, 4.8), theme="dark")
        self.figure_canvas = FigureCanvasQTAgg(figure)
        self.figure_canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self.figure_canvas.setStyleSheet(
            f"background-color: {self.COLORS['card']};"
        )
        self.canvas_layout.addWidget(self.figure_canvas)

        # 3) Adiciona a barra de navegação (zoom, salvar) já estilizada.
        self.graph_toolbar = BodeNavigationToolbar(
            self.figure_canvas, self.toolbar_holder
        )
        self.style_graph_toolbar()
        self.toolbar_layout.addWidget(self.graph_toolbar)

        # 4) Liga a leitura interativa dos pontos e desenha o gráfico.
        self.enable_graph_interaction()
        self.figure_canvas.draw()

    def style_graph_toolbar(self) -> None:
        """Aplica o modo noturno aos controles nativos da barra do Matplotlib."""

        self.apply_toolbar_style(self.graph_toolbar)

    def apply_toolbar_style(self, toolbar) -> None:
        """Estiliza qualquer barra de ferramentas do Matplotlib no tema escuro."""

        if toolbar is None:
            return
        c = self.COLORS
        toolbar.setStyleSheet(
            f"QToolBar {{ background: {c['card']}; border: none; spacing: 2px; }}"
            f"QToolButton {{ background: {c['card_soft']}; color: {c['text']};"
            f" border: none; padding: 4px; margin: 1px; }}"
            f"QToolButton:hover {{ background: {c['line']}; }}"
            f"QToolButton:checked {{ background: {c['blue']}; }}"
            f"QLabel {{ color: {c['muted']}; background: {c['card']}; }}"
        )

    def update_extra_views(self) -> None:
        """Recria as figuras das abas extras (circuito, sistema e análise)."""

        if self.project is None:
            return
        # Percorre circuito, sistema e análise, gerando uma figura por aba.
        for name, builder in self.extra_builders.items():
            figure = builder(self.project, theme="dark")
            self.swap_extra_canvas(name, figure)

    def swap_extra_canvas(self, name: str, figure: Figure) -> None:
        """Substitui o canvas e a barra de uma aba extra pela figura nova."""

        canvas_layout = self.extra_canvas_layouts[name]
        toolbar_layout = self.extra_toolbar_layouts[name]
        self.extra_placeholders[name].hide()

        # 1) Remove o canvas e a barra anteriores, se existirem.
        old_canvas = self.extra_canvases.get(name)
        if old_canvas is not None:
            canvas_layout.removeWidget(old_canvas)
            old_canvas.setParent(None)
            old_canvas.deleteLater()
        old_toolbar = self.extra_toolbars.get(name)
        if old_toolbar is not None:
            toolbar_layout.removeWidget(old_toolbar)
            old_toolbar.setParent(None)
            old_toolbar.deleteLater()

        # 2) Insere o novo canvas e a barra de ferramentas estilizada.
        canvas = FigureCanvasQTAgg(figure)
        canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        canvas.setStyleSheet(f"background-color: {self.COLORS['card']};")
        canvas_layout.addWidget(canvas)

        toolbar = BodeNavigationToolbar(canvas, self)
        self.apply_toolbar_style(toolbar)
        toolbar_layout.addWidget(toolbar)

        # 3) Guarda as referências e desenha a figura.
        self.extra_canvases[name] = canvas
        self.extra_toolbars[name] = toolbar
        canvas.draw()

    def clear_extra_views(self) -> None:
        """Remove as figuras das abas extras e mostra os avisos novamente."""

        # Para cada aba extra, descarta canvas/barra e restaura o placeholder.
        for name in self.extra_builders:
            canvas = self.extra_canvases.get(name)
            if canvas is not None:
                self.extra_canvas_layouts[name].removeWidget(canvas)
                canvas.setParent(None)
                canvas.deleteLater()
                self.extra_canvases[name] = None
            toolbar = self.extra_toolbars.get(name)
            if toolbar is not None:
                self.extra_toolbar_layouts[name].removeWidget(toolbar)
                toolbar.setParent(None)
                toolbar.deleteLater()
                self.extra_toolbars[name] = None
            self.extra_placeholders[name].show()

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
        selected, _ = QFileDialog.getSaveFileName(
            self,
            "Salvar relatório do BenaCCA",
            "Relatorio BenaCCA.pdf",
            "Documento PDF (*.pdf)",
        )
        if not selected:
            return
        try:
            pdf = generate_report_pdf(self.project, selected)
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, "Erro ao gerar relatório", str(exc))
            return
        QMessageBox.information(
            self,
            "Relatório concluído",
            f"PDF gerado com sucesso:\n\n{pdf.name}",
        )

    def reset(self) -> None:
        """Restaura os parâmetros padrão e limpa os resultados."""

        self.frequency_edit.setText("2000")
        self.impedance_edit.setText("8")
        self.tolerance_edit.setText("5")
        self.project = None
        self.results_table.setRowCount(0)
        for label in self.summary_value_labels:
            label.setText("—")
        self.disconnect_graph_interaction()
        if self.figure_canvas is not None:
            self.canvas_layout.removeWidget(self.figure_canvas)
            self.figure_canvas.setParent(None)
            self.figure_canvas.deleteLater()
            self.figure_canvas = None
        if self.graph_toolbar is not None:
            self.toolbar_layout.removeWidget(self.graph_toolbar)
            self.graph_toolbar.setParent(None)
            self.graph_toolbar.deleteLater()
            self.graph_toolbar = None
        self.graph_placeholder.show()
        self.update_methodology()
        self.report_button.setEnabled(False)


# =============================================================================
# INICIALIZAÇÃO
# =============================================================================


def main() -> None:
    """Inicia o aplicativo."""

    # No Qt6 o ajuste de High DPI já é automático, mantendo textos e gráficos
    # nítidos quando o sistema usa escala de tela.
    app = QApplication(sys.argv)
    app.setApplicationName("BenaCCA")
    default_font = QFont("Segoe UI", 10)
    app.setFont(default_font)

    window = BenaCCAApp()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
