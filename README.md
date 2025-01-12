# Material del Tabajo Fin de Máster

Los diferentes materiales se han distribuido en carpetas relacionadas con el manuscrito enviado en la plataforma de la universidad.

## Análisis de Supervivencia

El material se encuentra en la carpeta *OS_curves*, en el que se encuentran los archivos del código (.Rmd) y su pdf correspondiente. En el mismo script se han realizado todos los análisis de esta sección.

## Análisis de deconvolución

El material del análisis de deconvolución (carpeta *Deconvolution_analysis*) se ha distribuido en tres carpetas:

a) *Deconvolution_packages*: contiene los códigos de otros paquetes de R explorados para este trabajo.

b) *Bisque_analyses*: contiene archivos de códigos (.Rmd y su pdf), y los resultados de las fracciones celulares por cada estudio.

c) *quantiseqr_analyses*: contiene archivos de códigos (.Rmd y su pdf), y los resultados de las fracciones celulares por cada estudio. Además se han incluido los archivos correspondientes a las comparaciones de los resultados.

También se incluyen los archivos usados para presentar los resultados en heatmaps y los relacionados con los análisis estadísticos de la comparaciones de las fracciones celulares (csv, Rmd y pdf).

## Análisis estadístico

Se realizó un análisis estadístico (carpeta *Statistical_analysis*) de los resultados de la deconvolución con el código presentado en esta carpeta.

## Análisis de genes diferencialmente expresados (DEGs)

El material de este análisis (*DEG_analysis*) se ha distribuido en tres carpetas:

a) *Final_linear_models*: Contiene los archivos de los códigos utilizados (.Rmd) y los archivos Excel resultantes ("_stats.xlsx"), por cada cohorte seleccionada.

b) *GSE91061_vs_GSE50509*: Contiene los resultados de los solapamientos de los listados de DEGs resultantes de los modelos lineales de estas dos cohortes: los archivos .png muestran los Venn diagrams, y los archivos .csv contienen la distribución de los solapamientos.

c) *GSE91061_vs_GSE61992_vs_TCGASKCM*: Contiene los resultados de los solapamientos de los listados de DEGs resultantes de los modelos lineales de estas tres cohortes: los archivos .png muestran los Venn diagrams, y los archivos .csv contienen la distribución de los solapamientos.

## Análisis de pathways

Se presentan los códigos de los análisis de pathways en archivos Rmd y pdf por tipo celular.

La comparación de estos resultados se presentan dentro de dos carpetas:

a) *untreated_vs_treated*: consiste en la comparación de la cohorte GSE50509 con los estudios analizados, por cada tipo celular. Se presentan los resultados de los networks (pdf) obtenidas del análisis con EnrichmentMap en Cytoscape utilizando los archivos txt.

b) *treatments_comparison*: presenta los archivos txt de las comparaciones entre tratamientos. "Only_geneset" indica que el resultado corresponde al análisis de pathways usando los DEGs que nos se solapan entre tratamientos; "maxOverlap_geneset" indica que el resultado corresponde al análisis de pathways usando los DEGs que se solapan entre tratamientos.

## Anotaciones específicas

Los resultados de la exploración de anotaciones específicas se sitúan en esta carpeta (*Specific_annotations*).

El contenido consiste en heatmaps asociados a los pathways explorados (con los códigos Rmd y gráficos por pathway en pdf), y otras dos carpetas:

a) *Reference_annotations*: Presenta los archivos txt con los listados de genes de referencia descargados de la web [Gene Oncology con](https://geneontology.org/).

b) *Customized_GOBPs*: Contiene un archivo de Excel con el resumen de la comparación de los DEGs por cada estudio frente a los listados de referencia, la leyenda de los volcano plots, y carpetas de los estudios con imágenes de los volcano plots resultantes por cada tipo celular y pathway consultado.
