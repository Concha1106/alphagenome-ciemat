#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 14:47:54 2026

@author: e6260
"""

"""
Script: predict-contact-maps-sec23b.py
Description: Predice CONTACT_MAPS para la región SEC23B usando predict_interval().
@author: Concha
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from alphagenome.models import dna_client
from alphagenome.data import genome


# =========================
# 1. CONFIGURACIÓN
# =========================

OUTPUT_DIR = "results/contact-maps-sec23b"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Coordenadas de tu secuencia según Ensembl BLAT
chromosome = "chr20"
seq_start = 18510402
seq_end = 18511253

# Resolución aproximada de CONTACT_MAPS
BIN_SIZE = 2048


# =========================
# 2. CREAR CLIENTE
# =========================

api_key = os.environ.get("ALPHAGENOME_API_KEY")

if api_key is None:
    raise ValueError("No encuentro ALPHAGENOME_API_KEY.")

dna_model = dna_client.create(api_key)


# =========================
# 3. DEFINIR INTERVALO GENÓMICO
# =========================

interval = genome.Interval(
    chromosome=chromosome,
    start=seq_start,
    end=seq_end,
).resize(dna_client.SEQUENCE_LENGTH_1MB)

print("Intervalo enviado a AlphaGenome:")
print(interval)


# =========================
# 4. PREDICCIÓN CONTACT_MAPS
# =========================

output = dna_model.predict_interval(
    interval=interval,
    requested_outputs=[dna_client.OutputType.CONTACT_MAPS],
    ontology_terms=None,
)

contact_maps = output.contact_maps

print("Shape de CONTACT_MAPS:")
print(contact_maps.values.shape)

print("Metadata:")
print(contact_maps.metadata)

metadata_path = os.path.join(OUTPUT_DIR, "contact-maps-metadata.tsv")

contact_maps.metadata.to_csv(
    metadata_path,
    sep="\t",
    index=True
)

print(f"Metadata guardado en: {metadata_path}")

metadata_simple = contact_maps.metadata[[
    "ontology_curie",
    "biosample_name",
    "biosample_type",
    "Assay title",
    "name"
]]

print(metadata_simple)

metadata_simple.to_csv(
    os.path.join(OUTPUT_DIR, "contact-maps-metadata-simple.tsv"),
    sep="\t",
    index=True
)


# =========================
# 5. EXTRAER MATRIZ
# =========================

values = contact_maps.values

if values.ndim == 2:
    matrix = values

elif values.ndim == 3:
    track_index = 7
    matrix = values[:, :, track_index]
    print(f"Usando track_index = {track_index}")

else:
    raise ValueError(f"Dimensión inesperada: {values.shape}")

selected_metadata = contact_maps.metadata.iloc[track_index]

print("Track seleccionado:")
print(selected_metadata)
# =========================
# 6. CALCULAR BINS DE TU REGIÓN
# =========================

window_start = interval.start
window_end = interval.end

bin_start = int((seq_start - window_start) // BIN_SIZE)
bin_end = int(np.ceil((seq_end - window_start) / BIN_SIZE))

print("Tu secuencia cae aproximadamente en estos bins:")
print(f"bin_start = {bin_start}")
print(f"bin_end   = {bin_end}")

submatrix = matrix[bin_start:bin_end + 1, bin_start:bin_end + 1]

print("Shape de la submatriz de tu región:")
print(submatrix.shape)


# =========================
# 7. GUARDAR MATRIZ COMPLETA
# =========================

matrix_path = os.path.join(OUTPUT_DIR, "contact-maps-full-matrix.tsv")

pd.DataFrame(matrix).to_csv(
    matrix_path,
    sep="\t",
    index=False,
    header=False,
)

print(f"Matriz completa guardada en: {matrix_path}")


# =========================
# 8. GUARDAR SUBMATRIZ DE TU REGIÓN
# =========================

submatrix_path = os.path.join(OUTPUT_DIR, "contact-maps-sequence-region-submatrix.tsv")

pd.DataFrame(submatrix).to_csv(
    submatrix_path,
    sep="\t",
    index=False,
    header=False,
)

print(f"Submatriz de tu región guardada en: {submatrix_path}")




# =========================
# 9. VISUALIZAR MATRIZ COMPLETA
# =========================

fig_path = os.path.join(OUTPUT_DIR, "contact-maps-full-heatmap.png")

plt.figure(figsize=(7, 6))
plt.imshow(matrix, origin="lower", aspect="auto")
plt.colorbar(label="Predicted contact frequency")
plt.title("AlphaGenome CONTACT_MAPS - SEC23B 1 Mb window")
plt.xlabel("Genomic bin")
plt.ylabel("Genomic bin")
plt.tight_layout()
plt.savefig(fig_path, dpi=300)
plt.close()

print(f"Figura completa guardada en: {fig_path}")


# =========================
# 10. VISUALIZAR SUBMATRIZ
# =========================

subfig_path = os.path.join(OUTPUT_DIR, "contact-maps-sequence-region-heatmap.png")

plt.figure(figsize=(5, 4))
plt.imshow(submatrix, origin="lower", aspect="auto")
plt.colorbar(label="Predicted contact frequency")
plt.title("CONTACT_MAPS - sequence region")
plt.xlabel("Genomic bin")
plt.ylabel("Genomic bin")
plt.tight_layout()
plt.savefig(subfig_path, dpi=300)
plt.close()

print(f"Figura de submatriz guardada en: {subfig_path}")


rows = []

for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):

        bin_i_start = window_start + i * BIN_SIZE
        bin_i_end = bin_i_start + BIN_SIZE - 1

        bin_j_start = window_start + j * BIN_SIZE
        bin_j_end = bin_j_start + BIN_SIZE - 1

        rows.append({
            "track_index": track_index,
            "ontology_curie": selected_metadata["ontology_curie"],
            "biosample_name": selected_metadata["biosample_name"],
            "bin_i": i,
            "chr_i": chromosome,
            "start_i": bin_i_start,
            "end_i": bin_i_end,
            "bin_j": j,
            "chr_j": chromosome,
            "start_j": bin_j_start,
            "end_j": bin_j_end,
            "contact_value": matrix[i, j],
        })

contact_table = pd.DataFrame(rows)

contact_table_path = os.path.join(
    OUTPUT_DIR,
    "contact-maps-full-table-with-coordinates.tsv"
)

contact_table.to_csv(
    contact_table_path,
    sep="\t",
    index=False
)

print(f"Tabla completa con coordenadas guardada en: {contact_table_path}")


# =========================
# 12. CONTACTOS DE LA REGIÓN DE INTERÉS
# =========================

region_contacts = contact_table[
    (contact_table["bin_i"].between(bin_start, bin_end)) |
    (contact_table["bin_j"].between(bin_start, bin_end))
].copy()

region_contacts_path = os.path.join(
    OUTPUT_DIR,
    "contact-maps-region-contacts-with-coordinates.tsv"
)

region_contacts.to_csv(
    region_contacts_path,
    sep="\t",
    index=False
)

print(f"Contactos de la región guardados en: {region_contacts_path}")

# =========================
# 13. TOP CONTACTOS
# =========================

top_contacts = region_contacts.sort_values(
    "contact_value",
    ascending=False
).head(50)

top_contacts_path = os.path.join(
    OUTPUT_DIR,
    "contact-maps-region-top50-contacts.tsv"
)

top_contacts.to_csv(
    top_contacts_path,
    sep="\t",
    index=False
)

print(top_contacts)

print(f"Top contactos guardados en: {top_contacts_path}")
