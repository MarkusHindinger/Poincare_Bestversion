import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')  # Setze das Backend explizit auf TkAgg
from matplotlib.patches import Arc
from matplotlib.path import Path
import matplotlib.patches as patches
import math

# Konfiguriere matplotlib für interaktive Anzeige
plt.ion()

class PoincareHalfPlane:
    """
    Klasse zur Berechnung von hyperbolischer Geometrie in der Poincaré-Halbebene.
    """
    
    def __init__(self):
        self.epsilon = 1e-10
    
    def check_valid_point(self, point):
        """Überprüft, ob ein Punkt in der Poincaré-Halbebene liegt (y > 0)."""
        x, y = point
        if y <= 0:
            raise ValueError(f"Punkt {point} liegt nicht in der Poincaré-Halbebene (y muss > 0 sein)")
        return True
    
    def hyperbolic_distance(self, p1, p2):
        """Berechnet den hyperbolischen Abstand zwischen zwei Punkten."""
        self.check_valid_point(p1)
        self.check_valid_point(p2)
        
        x1, y1 = p1
        x2, y2 = p2
        
        term = 1 + ((x2-x1)**2 + (y2-y1)**2) / (2*y1*y2)
        return np.arccosh(term)
    
    def geodesic(self, p1, p2, num_points=100):
        """Berechnet die Geodätische zwischen zwei Punkten."""
        self.check_valid_point(p1)
        self.check_valid_point(p2)
        
        x1, y1 = p1
        x2, y2 = p2
        
        if abs(x1 - x2) < self.epsilon:
            y_values = np.linspace(min(y1, y2), max(y1, y2), num_points)
            x_values = np.ones_like(y_values) * x1
            return np.column_stack((x_values, y_values))
        
        h = ((x1**2 + y1**2) - (x2**2 + y2**2)) / (2 * (x1 - x2))
        r = np.sqrt((x1 - h)**2 + y1**2)
        
        theta1 = np.arctan2(y1, x1 - h)
        theta2 = np.arctan2(y2, x2 - h)
        
        if theta1 > theta2:
            theta1, theta2 = theta2, theta1
            
        if abs(theta2 - theta1) > np.pi:
            if theta1 < 0:
                theta1 += 2 * np.pi
            else:
                theta2 += 2 * np.pi
        
        thetas = np.linspace(theta1, theta2, num_points)
        x_values = h + r * np.cos(thetas)
        y_values = r * np.sin(thetas)
        
        return np.column_stack((x_values, y_values))
    
    def parallel_transport(self, p1, v1, p2):
        """Führt den Paralleltransport eines Vektors durch."""
        self.check_valid_point(p1)
        self.check_valid_point(p2)
        
        x1, y1 = p1
        x2, y2 = p2
        vx1, vy1 = v1
        
        hyperbolic_norm_v1 = np.sqrt((vx1**2 + vy1**2) / y1**2)
        vx1, vy1 = vx1/hyperbolic_norm_v1, vy1/hyperbolic_norm_v1
        
        if abs(x1 - x2) < self.epsilon:
            tangent_geo_p1 = np.array([0, 1])
            tangent_geo_p2 = np.array([0, 1])
        else:
            h = ((x1**2 + y1**2) - (x2**2 + y2**2)) / (2 * (x1 - x2))
            r = np.sqrt((x1 - h)**2 + y1**2)
            
            r1_vec = np.array([x1 - h, y1])
            r1_len = np.sqrt(r1_vec[0]**2 + r1_vec[1]**2)
            r1_norm = r1_vec / r1_len
            tangent_geo_p1 = np.array([-r1_norm[1], r1_norm[0]])
            
            r2_vec = np.array([x2 - h, y2])
            r2_len = np.sqrt(r2_vec[0]**2 + r2_vec[1]**2)
            r2_norm = r2_vec / r2_len
            tangent_geo_p2 = np.array([-r2_norm[1], r2_norm[0]])
        
        dot_product = vx1 * tangent_geo_p1[0] + vy1 * tangent_geo_p1[1]
        norm_v1 = np.sqrt(vx1**2 + vy1**2)
        norm_t1 = np.sqrt(tangent_geo_p1[0]**2 + tangent_geo_p1[1]**2)
        cos_alpha = np.clip(dot_product / (norm_v1 * norm_t1), -1.0, 1.0)
        alpha = np.arccos(cos_alpha)
        
        cross_product = vx1 * tangent_geo_p1[1] - vy1 * tangent_geo_p1[0]
        
        if cross_product < 0:
            alpha = 2 * np.pi - alpha
        
        norm_t2 = np.sqrt(tangent_geo_p2[0]**2 + tangent_geo_p2[1]**2)
        tangent_geo_p2 = tangent_geo_p2 / norm_t2
        
        cos_alpha = np.cos(alpha)
        sin_alpha = np.sin(alpha)
        
        ortho_geo_p2 = np.array([tangent_geo_p2[1], -tangent_geo_p2[0]])
        v2 = tangent_geo_p2 * cos_alpha + ortho_geo_p2 * sin_alpha
        
        hyperbolic_norm_v2 = np.sqrt((v2[0]**2 + v2[1]**2) / y2**2)
        v2 = v2 / hyperbolic_norm_v2
        
        return v2

    def angle_between(self, v1, v2, y):
        """Berechnet den Winkel zwischen zwei Vektoren."""
        dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
        return np.arccos(dot)
    
    def vector_distance(self, p1, v1, p2, v2):
        """Berechnet die Distanz zwischen zwei Vektoren."""
        self.check_valid_point(p1)
        self.check_valid_point(p2)

        # Normiere die Vektoren bezüglich der hyperbolischen Metrik
        norm1 = np.sqrt((v1[0]**2 + v1[1]**2) / p1[1]**2)
        norm2 = np.sqrt((v2[0]**2 + v2[1]**2) / p2[1]**2)
        v1 = v1 / norm1
        v2 = v2 / norm2

        # Transportiere v1 nach p2 und berechne den Winkel
        v1_transported = self.parallel_transport(p1, v1, p2)
        angle = self.angle_between(v1_transported, v2, p2[1])
        
        # Berechne die Gesamtdistanz
        d = self.hyperbolic_distance(p1, p2)
        return np.sqrt(d**2 + angle**2)

class PoincareVisualizer:
    """
    Klasse für die Visualisierung in der Poincaré-Halbebene.
    """
    
    def __init__(self, poincare):
        self.poincare = poincare
    
    def draw_vector(self, ax, base, vec, color, label, scale=0.5):
        """Zeichnet einen Vektor."""
        ax.arrow(base[0], base[1], vec[0]*scale, vec[1]*scale,
                head_width=0.05, head_length=0.08,
                fc=color, ec=color, label=label)
    
    def draw_geodesic(self, ax, p1, p2, color, label):
        """Zeichnet eine Geodätische."""
        geo = self.poincare.geodesic(p1, p2)
        ax.plot(geo[:, 0], geo[:, 1], linestyle='--', color=color, label=label)
    
    def plot_basic_geodesic(self, p1, p2):
        """Erstellt eine einfache Geodätische-Visualisierung."""
        fig, ax = plt.subplots(figsize=(10, 6))
        geo_points = self.poincare.geodesic(p1, p2)
        
        ax.plot(geo_points[:, 0], geo_points[:, 1], 'b-', linewidth=2)
        ax.scatter([p1[0], p2[0]], [p1[1], p2[1]], color='red', s=50)
        ax.text(p1[0], p1[1], f'$p_1$', fontsize=12)
        ax.text(p2[0], p2[1], f'$p_2$', fontsize=12)
        
        self._setup_plot(ax, 'Geodätische in der Poincaré-Halbebene', [p1, p2])
        plt.show(block=True)

    def visualize_parallel_transport(self, p1, p2, v1, num_steps=10):
        """
        Visualisiert den Paralleltransport eines Vektors entlang einer Geodäte.
        
        Args:
            p1: Startpunkt (NumPy-Array)
            p2: Endpunkt (NumPy-Array)
            v1: Startvektor (NumPy-Array)
            num_steps: Anzahl der Zwischenschritte
        """
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Berechne die Geodätische und die Punkte darauf
        geo_points = self.poincare.geodesic(p1, p2)
        step_indices = np.linspace(0, len(geo_points)-1, num_steps, dtype=int)
        intermediate_points = [np.array(geo_points[i]) for i in step_indices]  # NumPy-Arrays
        intermediate_points = intermediate_points[::-1]  # Kehre die Reihenfolge um
        
        # Zeichne die Geodätische
        ax.plot(geo_points[:, 0], geo_points[:, 1], 'b-', linewidth=2, label='Geodäte')
        
        # Sammle alle Vektoren für die Plotgrenzen
        vectors_for_bounds = []
        display_length = 0.8
        
        # Normiere den Startvektor bezüglich der hyperbolischen Metrik
        v1 = np.array(v1)  # Stelle sicher, dass es ein NumPy-Array ist
        v1_norm = np.sqrt((v1[0]**2 + v1[1]**2) / p1[1]**2)
        v1 = v1 / v1_norm
        
        # Transportiere den Vektor schrittweise entlang der Geodäte
        current_point = intermediate_points[0]  # Starte bei p1
        current_vector = v1
        
        # Zeichne die Vektoren an jedem Zwischenpunkt
        for i, point in enumerate(intermediate_points):
            if i == 0:
                # Am Startpunkt verwenden wir den originalen Vektor
                color = 'red'
                label = 'Startvektor'
            else:
                # Für jeden weiteren Punkt transportieren wir vom vorherigen Punkt
                current_vector = self.poincare.parallel_transport(intermediate_points[i-1], current_vector, point)
                
                if i == len(intermediate_points) - 1:
                    color = 'green'
                    label = 'Endvektor'
                else:
                    color = 'orange'
                    label = 'Zwischenvektor' if i == 1 else None
            
            # Skaliere und zeichne den Vektor
            scaled_vec = self._scale_vector(current_vector, point[1], display_length)
            self.draw_vector(ax, point, scaled_vec, color, label)
            vectors_for_bounds.append((point, scaled_vec, 1.0))
            
            # Markiere den Punkt
            if i == 0:
                ax.plot(*point, 'ro', markersize=10, label='Startpunkt')
            elif i == len(intermediate_points) - 1:
                ax.plot(*point, 'go', markersize=10, label='Endpunkt')
            else:
                ax.plot(*point, 'ko', markersize=6)
        
        # Konfiguriere den Plot
        self._setup_plot(ax, 'Paralleltransport entlang einer Geodäte', 
                        intermediate_points, vectors_for_bounds)
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.show(block=True)

    def _scale_vector(self, vec, base_y, display_length):
        """
        Skaliert einen Vektor für die Anzeige.
        
        Args:
            vec: NumPy-Array des zu skalierenden Vektors
            base_y: y-Koordinate des Basispunkts
            display_length: Gewünschte Anzeigelänge
        """
        vec = np.array(vec)  # Stelle sicher, dass es ein NumPy-Array ist
        norm = np.linalg.norm(vec)
        if norm == 0:
            return vec
        scaling = display_length * base_y / np.sqrt(vec[0]**2 + vec[1]**2)
        return vec * scaling

    def _setup_plot(self, ax, title, points=None, vectors=None):
        """
        Konfiguriert die grundlegenden Plot-Einstellungen.
        
        Args:
            ax: Die Matplotlib-Axes
            title: Der Titel des Plots
            points: Liste von Punkten, die im Plot dargestellt werden
            vectors: Liste von Tupeln (Basispunkt, Vektor, Skalierung) für die Vektoren
        """
        ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        ax.grid(True, alpha=0.3)
        ax.set_title(title)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
        
        if points is not None:
            # Sammle alle relevanten Koordinaten
            x_coords = [p[0] for p in points]
            y_coords = [p[1] for p in points]
            
            # Füge Vektorendpunkte hinzu, wenn vorhanden
            if vectors is not None:
                for base, vec, scale in vectors:
                    end_x = base[0] + vec[0] * scale
                    end_y = base[1] + vec[1] * scale
                    x_coords.append(end_x)
                    y_coords.append(end_y)
            
            x_min, x_max = min(x_coords), max(x_coords)
            y_min, y_max = 0, max(y_coords)  # y_min bleibt 0 für die Halbebene
            
            # Füge großzügigeres Padding hinzu (30% des Bereichs)
            x_range = max(x_max - x_min, 1.0)  # Mindestens 1.0 als Range
            y_range = max(y_max - y_min, 1.0)
            
            padding_x = x_range * 0.3
            padding_y = y_range * 0.3
            
            ax.set_xlim(x_min - padding_x, x_max + padding_x)
            ax.set_ylim(y_min, y_max + padding_y)
        else:
            ax.set_ylim(0, 3)

class PoincareTests:
    """
    Testklasse für die Poincaré-Halbebene.
    """
    
    def __init__(self, poincare):
        self.poincare = poincare
        self.visualizer = PoincareVisualizer(poincare)
    
    def test_triangle_inequality(self, max_trials=None):
        """
        Testet die Dreiecksungleichung für Vektortransporte:
        d(v1,v3) ≤ d(v1,v2) + d(v2,v3)
        wobei v1 nach p2 und p3, und v2 nach p3 transportiert werden.
        """
        p1 = np.array([1.0, 1.0])
        v1 = np.array([1.0, 0.0])
        p3 = np.array([3.0, 2.0])
        v3 = np.array([0.0, 1.0])

        trial_count = 0
        while max_trials is None or trial_count < max_trials:
            trial_count += 1
            p2 = np.random.uniform(0.5, 5, size=2)
            p2[1] = abs(p2[1]) + 0.1
            v2 = np.random.randn(2)

            print(f"\n--- Versuch {trial_count} ---")
            print(f"p2 = {p2}, v2 = {v2}")

            # Berechne d(v1,v2): Transport von v1 nach p2 und Vergleich mit v2
            d12 = self.poincare.vector_distance(p1, v1, p2, v2)
            
            # Berechne d(v2,v3): Transport von v2 nach p3 und Vergleich mit v3
            d23 = self.poincare.vector_distance(p2, v2, p3, v3)
            
            # Berechne d(v1,v3): Transport von v1 nach p3 und Vergleich mit v3
            d13 = self.poincare.vector_distance(p1, v1, p3, v3)

            print(f"d(v1,v2): {d12:.4f}")
            print(f"d(v2,v3): {d23:.4f}")
            print(f"d(v1,v3): {d13:.4f}")
            print(f"d(v1,v2) + d(v2,v3) = {d12 + d23:.4f}")
            print(f"Verletzung = {d13 - (d12 + d23):.6f}")

            if d13 > d12 + d23 + 1e-6:
                print("❌ Verstoß gegen Dreiecksungleichung gefunden!")
                self._visualize_violation(p1, v1, p2, v2, p3, v3)
                return True
            
            if trial_count % 100 == 0:
                print(f"✅ Keine Verstöße nach {trial_count} Versuchen.")
        
        print(f"✅ Keine Verstöße gegen die Dreiecksungleichung in {trial_count} Versuchen gefunden.")
        return False

    def _visualize_violation(self, p1, v1, p2, v2, p3, v3):
        """Visualisiert eine Verletzung der Dreiecksungleichung."""
        fig, ax = plt.subplots(figsize=(12, 8))  # Größere Figur
        
        # Vektoren und ihre Transporte
        norm_v1 = v1 / np.sqrt((v1[0]**2 + v1[1]**2) / p1[1]**2)
        norm_v2 = v2 / np.sqrt((v2[0]**2 + v2[1]**2) / p2[1]**2)
        
        # Transporte berechnen
        v1_to_p2 = self.poincare.parallel_transport(p1, norm_v1, p2)  # v1 nach p2
        v1_to_p3 = self.poincare.parallel_transport(p1, norm_v1, p3)  # v1 nach p3
        v2_to_p3 = self.poincare.parallel_transport(p2, norm_v2, p3)  # v2 nach p3

        # Geodätische Linien
        self.visualizer.draw_geodesic(ax, p1, p2, 'gray', 'Geodäte p1→p2')
        self.visualizer.draw_geodesic(ax, p2, p3, 'black', 'Geodäte p2→p3')
        self.visualizer.draw_geodesic(ax, p1, p3, 'darkviolet', 'Geodäte p1→p3')

        # Vektoren zeichnen
        display_length = 0.8  # Erhöhte Vektorlänge
        vector_bases_and_data = [
            (p1, v1, 'red', 'v1 an p1'),
            (p2, v2, 'green', 'v2 an p2'),
            (p3, v3, 'blue', 'v3 an p3'),
            (p2, v1_to_p2, 'orange', 'v1→p2'),
            (p3, v1_to_p3, 'purple', 'v1→p3'),
            (p3, v2_to_p3, 'brown', 'v2→p3')
        ]
        
        # Sammle Vektordaten für die Plotgrenzen
        vectors_for_bounds = []
        for base, vec, color, label in vector_bases_and_data:
            scaled_vec = self._scale_vector(vec, base[1], display_length)
            self.visualizer.draw_vector(ax, base, scaled_vec, color, label)
            vectors_for_bounds.append((base, scaled_vec, 1.0))

        # Punkte markieren
        ax.plot(*p1, 'ro', label='p1', markersize=10)
        ax.plot(*p2, 'go', label='p2', markersize=10)
        ax.plot(*p3, 'bo', label='p3', markersize=10)

        # Konfiguriere Plot mit allen relevanten Punkten und Vektoren
        all_points = [p1, p2, p3]
        self.visualizer._setup_plot(ax, "Transportierte Vektoren & Dreiecksungleichung", 
                                  all_points, vectors_for_bounds)
        
        # Legende außerhalb des Plots platzieren
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()  # Automatische Anpassung des Layouts
        plt.show(block=True)  # Blockiere die Ausführung, bis das Fenster geschlossen wird
    
    def _scale_vector(self, vec, base_y, display_length):
        """Skaliert einen Vektor für die Anzeige."""
        norm = np.linalg.norm(vec)
        if norm == 0:
            return vec
        scaling = display_length * base_y / np.sqrt(vec[0]**2 + vec[1]**2)
        return vec * scaling

if __name__ == "__main__":
    # Erstelle Instanzen
    poincare = PoincareHalfPlane()
    visualizer = PoincareVisualizer(poincare)
    tests = PoincareTests(poincare)
    
    # Beispiel für Geodätische
    p1 = np.array([1, 1])
    p2 = np.array([3, 2])
    visualizer.plot_basic_geodesic(p1, p2)
    
    # Warte auf Benutzereingabe bevor der nächste Plot gezeigt wird
    input("Drücken Sie Enter für den nächsten Plot...")
    
    # Beispiel für Paralleltransport
    v1 = np.array([-1.0, -0.5])
    visualizer.visualize_parallel_transport(p1, p2, v1, num_steps=8)
    
    # Warte auf Benutzereingabe bevor der Test startet
    input("Drücken Sie Enter für den Dreiecksungleichungstest...")
    
    # Teste Dreiecksungleichung (läuft bis eine Verletzung gefunden wird)
    tests.test_triangle_inequality() 
    