import os
from alphagenome.models import dna_client

def get_dna_model():
    api_key = os.getenv("ALPHAGENOME_API_KEY")
    if not api_key:  # For environments such as Spyder, also read the key from .bashrc
        try:
            with open(os.path.expanduser("~/.bashrc")) as f:
                for line in f:
                    if "ALPHAGENOME_API_KEY" in line:
                        api_key = line.split("=")[1].strip().strip('"')
                        break
        except Exception:
            pass

    if not api_key:
        raise ValueError("ALPHAGENOME_API_KEY was not found in the environment or in .bashrc")

    return dna_client.create(api_key, timeout=30)	
