# Evaluación de Long Branch Attraction (LBA) en *aeneus*

Este repositorio contiene los scripts y resultados del análisis de 
posibles artefactos de atracción de ramas largas (LBA) en filogenias 
de SNPs de *aeneus*, usando dos enfoques complementarios:

1. **Distancias raíz-punta** (R + adephylo): identificación de tips con 
   distancias anómalas respecto a la media + 2 SD.
2. **TreeShrink**: detección estadística de ramas anómalas.

## Estructura del repositorio

- `Rcodes/disttoRoot.R` — Script principal de análisis en R.
- `cladosInfo.tsv` — Asignación de muestras a clados monofiléticos.
- `asc/` — Resultados de IQ-TREE con corrección ASC.
- `noasc/` — Resultados de IQ-TREE sin corrección ASC.
- `treeshrink_asc_pergene/` — Resultados de TreeShrink (ASC).
- `treeshrink_noasc/` — Resultados de TreeShrink (NOASC).

## Requisitos

### Software
- R ≥ 4.0 con paquetes: `ape`, `adephylo`, `pheatmap`, `ggplot2`, 
  `ggrepel`, `dplyr`
- IQ-TREE multicore version 2.1.4-beta
- TreeShrink version 1.4.0

### Datos de entrada
Los alineamientos no están incluidos por su tamaño. Disponibles en:
- Zenodo: [DOI pendiente]

## Cómo reproducir el análisis

### 1. Filogenia con IQ-TREE

```bash
#Para reconstruccion sin ASC
iqtree \
    -s "$workdir/aeneusLD_filter-subsample.min200.phy" \
    --seqtype DNA \
    --prefix amex-RAD-NOASC \
    -m GTR+G \
    -alrt 1000 \
    -B 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --redo \
    --ninit 100 \
    --nbest 10 \
    --bnni \
    --ntop 50
#Para reconstruccion con ASC
iqtree \
    -s $workdir/amex-RAD-asc.varsites.phy \
    --seqtype DNA \
    --prefix amex-RAD-asc \
    -m GTR+G+ASC \
    -alrt 1000 \
    -B 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --redo \
    --ninit 100 \
    --nbest 10 \
    --bnni \
    --ntop 50
```

### 2. Análisis de distancias raíz-punta
```r
source("Rcodes/disttoRoot.R")
```

### 3. TreeShrink
```bash
run_treeshrink.py -t asc/amex-RAD-asc.treefile -o treeshrink_asc_pergene
run_treeshrink.py -t noasc/amex-RAD-NOASC.treefile -o treeshrink_noasc
```

## Resultados principales

- El Clado 1 muestra distancias raíz-punta sistemáticamente mayores 
  (fig_boxplot_clados.png).
- TreeShrink identifica a ElLimon_T0495, T0496 y T0500 como ramas 
  anómalas en el análisis NOASC.
- Las firmas de TreeShrink señalan un bloque coherente dentro del 
  Clado 1, sugiriendo un efecto colectivo de LBA.

## Comparación de las ramas utilizando la prueba AU
dentro del Clado 1 vamos a comparar la reconstrucción original con otra, donde forzamos
la monofilia de dos de los subclados en el clado 1, y posteriormente evaluaremos
cual de las topologías es mejor con iqtree

```bash
# Concatenar las dos topologías (una por línea)
cat asc/amex-RAD-asc.treefile > topologias.tre
cat asc/amex-RAD-asc-constrained.treefile >> topologias.tre

# Prueba AU + KH + SH
iqtree \
    -s "$workdir/amex-RAD-asc.varsites.phy" \
    -m GTR+G+ASC \
    -z topologias.tre \
    -zb 10000 -au -zw \
    -T "$SLURM_CPUS_PER_TASK" \
    --prefix asc/amex-RAD-AUtest
```


## Contacto
Marco G. — [mgpastos@gmail.com]
