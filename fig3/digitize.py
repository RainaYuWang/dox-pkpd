"""Extract Fig. 3 from the archived PDF and regenerate observation CSV and overlay."""
import csv
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
PDF = HERE.parent / 'references/papers/41598_2017_Article_5902.pdf'
X_TICKS = [(177, 469), (729, 1021), (1282, 1573)]
Y_TICKS = [(327, 86), (663, 421), (998, 757)]
DOSES = [19.53125, 156.25, 1250]
DURATIONS = [6, 12, 24]


def main():
    with TemporaryDirectory() as tmp:
        prefix = Path(tmp) / 'fig3'
        subprocess.run(['pdfimages', '-f', '5', '-l', '5', '-png', PDF, prefix], check=True)
        source = Image.open(f'{prefix}-000.png').convert('RGB')
        assert source.size == (1650, 1101)
    with (HERE / "manual_points.csv").open(newline="") as file:
        manual = list(csv.DictReader(file))
    draw = ImageDraw.Draw(source)
    with (HERE / 'observations.csv').open('w', newline='') as file:
        writer = csv.writer(file, lineterminator='\n')
        writer.writerow(['panel', 'dose_nM', 'wash_h', 'time_h', 'cells', 'x_px', 'y_px', 'method'])
        for index, line in enumerate(line for line in (HERE / 'picks.txt').read_text().splitlines()
                                     if line and not line.startswith('#')):
            row, col = divmod(index, 3)
            x0, x600 = X_TICKS[col]
            y0, y40k = Y_TICKS[row]
            panel, *picks = line.split()
            rows = []
            for pick in picks:
                x, y = map(int, pick.rstrip('*').split(','))
                kind = 'ci_midpoint' if pick.endswith('*') else 'center'
                rows.append([panel, DOSES[row], DURATIONS[col], f'{600*(x-x0)/(x600-x0):.6f}',
                                 f'{40000*(y0-y)/(y0-y40k):.6f}', x, y, kind])
                color = '#1e68ff' if kind == 'ci_midpoint' else '#00996d'
                draw.ellipse((x-4, y-4, x+4, y+4), outline=color, width=1)
            for point in manual:
                if point['panel'] != panel:
                    continue
                t, cells = int(point['time_h']), int(point['cells'])
                rows.append([panel, DOSES[row], DURATIONS[col], t, cells, '', '', 'user_marked'])
                # The screenshots have their own pixel coordinates. Project only for the overlay.
                x = x0 + t * (x600-x0) / 600
                y = y0 - cells * (y0-y40k) / 40000
                draw.ellipse((x-3, y-3, x+3, y+3), outline='#9a168a', width=1)
            writer.writerows(sorted(rows, key=lambda r: float(r[3])))
    source.save(HERE / 'point_audit.png')


if __name__ == '__main__':
    main()
