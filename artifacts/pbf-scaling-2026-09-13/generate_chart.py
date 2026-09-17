import csv
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
rows = list(csv.DictReader((ROOT / "results.csv").open(newline="", encoding="utf-8")))

width, height = 1280, 760
left, top, right, bottom = 105, 70, 45, 105
plot_w = width - left - right
plot_h = height - top - bottom
image = Image.new("RGB", (width, height), "#f7f9fc")
draw = ImageDraw.Draw(image)

try:
    title_font = ImageFont.truetype("arialbd.ttf", 30)
    body_font = ImageFont.truetype("arial.ttf", 18)
    small_font = ImageFont.truetype("arial.ttf", 15)
except OSError:
    title_font = body_font = small_font = ImageFont.load_default()

x_min = math.log10(64)
x_max = math.log10(17576)
y_max = 1000.0


def x_pos(particles: float) -> float:
    return left + (math.log10(particles) - x_min) / (x_max - x_min) * plot_w


def y_pos(fps: float) -> float:
    return top + plot_h - fps / y_max * plot_h


draw.text((left, 20), "PBF frame-rate scaling with convex-hull rendering", fill="#172033", font=title_font)
for fps in range(0, 1001, 100):
    y = y_pos(fps)
    draw.line((left, y, width - right, y), fill="#dbe2ec", width=1)
    label = str(fps)
    box = draw.textbbox((0, 0), label, font=small_font)
    draw.text((left - 14 - (box[2] - box[0]), y - 8), label, fill="#596579", font=small_font)

particles = [int(row["particles"]) for row in rows]
for value in particles:
    x = x_pos(value)
    draw.line((x, top, x, top + plot_h), fill="#e8edf4", width=1)
    label = f"{value/1000:.1f}k" if value >= 10000 else f"{value:,}"
    box = draw.textbbox((0, 0), label, font=small_font)
    draw.text((x - (box[2] - box[0]) / 2, top + plot_h + 14), label, fill="#596579", font=small_font)

draw.line((left, top, left, top + plot_h), fill="#536075", width=2)
draw.line((left, top + plot_h, width - right, top + plot_h), fill="#536075", width=2)

sixty_y = y_pos(60)
draw.line((left, sixty_y, width - right, sixty_y), fill="#bd5b5b", width=2)
draw.text((width - right - 110, sixty_y - 24), "60 FPS", fill="#a44343", font=small_font)

series = [
    ("GPU GL4.3", "gpu_fps", "#2676d9"),
    ("CPU multi-thread", "cpu_mt_fps", "#e2862c"),
    ("CPU single-thread", "cpu_st_fps", "#2f9b67"),
]
for name, key, color in series:
    points = [(x_pos(float(row["particles"])), y_pos(float(row[key]))) for row in rows]
    draw.line(points, fill=color, width=4, joint="curve")
    for x, y in points:
        draw.ellipse((x - 6, y - 6, x + 6, y + 6), fill=color, outline="white", width=2)

legend_x = left + 20
legend_y = top + 20
for index, (name, _, color) in enumerate(series):
    y = legend_y + index * 30
    draw.line((legend_x, y + 9, legend_x + 34, y + 9), fill=color, width=4)
    draw.text((legend_x + 44, y), name, fill="#243047", font=body_font)

x_label = "Fluid particles (log scale)"
x_box = draw.textbbox((0, 0), x_label, font=body_font)
draw.text(((width - (x_box[2] - x_box[0])) / 2, height - 42), x_label, fill="#36445c", font=body_font)
draw.text((18, top + plot_h / 2 - 10), "FPS", fill="#36445c", font=body_font)
draw.text((left, height - 72), "120 fixed steps; first 30 discarded; 1280x720 offscreen hull pass", fill="#66738a", font=small_font)

image.save(ROOT / "pbf-fps-scaling.png")
