# Poincaré-Halbebene Visualisierung

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MarkusHindinger/Poincare_Bestversion/main)

Dieses Repository enthält eine Python-Implementierung zur Visualisierung der hyperbolischen Geometrie in der Poincaré-Halbebene.

## 🚀 Direkt ausprobieren!

Klicken Sie auf den "launch binder" Button oben, um den Code direkt im Browser auszuführen - keine Installation erforderlich!

## Features

- **Geodätische Visualisierung**: Zeigt die kürzesten Verbindungen zwischen Punkten in der hyperbolischen Geometrie
- **Paralleltransport**: Visualisiert, wie Vektoren entlang geodätischer Linien transportiert werden
- **Interaktive Plots**: Nutzt Matplotlib für klare, interaktive Visualisierungen
- **Mathematische Korrektheit**: Implementiert die exakten mathematischen Formeln der hyperbolischen Geometrie

## Mathematischer Hintergrund

Die Poincaré-Halbebene ist ein Modell der hyperbolischen Geometrie. Sie besteht aus der oberen Halbebene der komplexen Zahlen (y > 0) mit einer speziellen Metrik.

- **Geodätische**: In diesem Modell sind Geodätische entweder vertikale Geraden oder Halbkreise, die senkrecht zur x-Achse stehen
- **Paralleltransport**: Vektoren werden entlang geodätischer Linien parallel verschoben, wobei sich ihre Richtung und Länge in Bezug auf die hyperbolische Metrik ändert

## Dateien im Repository

- `Poincare_Bestversion.py`: Die Hauptimplementierung der Klassen und Funktionen
- `Poincare_Tutorial.ipynb`: Ein interaktives Jupyter Notebook mit Erklärungen und Beispielen
- `requirements.txt`: Liste der benötigten Python-Pakete

## Lokale Installation

Wenn Sie den Code lokal ausführen möchten:

```bash
# Repository klonen
git clone https://github.com/MarkusHindinger/Poincare_Bestversion.git
cd Poincare_Bestversion

# Abhängigkeiten installieren
pip install -r requirements.txt

# Code ausführen
python Poincare_Bestversion.py
``` 