"""
run_all.py - Ejecuta el pipeline completo en orden.

Uso (desde la raíz del repositorio):
    python scripts/run_all.py
"""
import importlib
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

STEPS = ["01_preprocesamiento", "02_analisis_univariado",
         "03_analisis_multivariado", "04_rutas_metabolicas"]

if __name__ == "__main__":
    t0 = time.time()
    for step in STEPS:
        importlib.import_module(step).main()
        print()
    print(f"Pipeline completo en {time.time() - t0:.1f} s. "
          "Resultados en results/tables y results/figures.")
