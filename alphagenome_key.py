import os
from alphagenome.models import dna_client

def get_dna_model():
    api_key = os.getenv("ALPHAGENOME_API_KEY")
    if not api_key: #Para entornos como Spyder, que pueden no encontrarla, leer desde .bashrc
        try:
            with open(os.path.expanduser("~/.bashrc")) as f:
                for line in f:
                    if "ALPHAGENOME_API_KEY" in line:
                        api_key = line.split("=")[1].strip().strip('"')
                        break
        except Exception:
            pass

    if not api_key:
        raise ValueError("No se encontró ALPHAGENOME_API_KEY ni en entorno ni en .bashrc")

    return dna_client.create(api_key, timeout=30)	
