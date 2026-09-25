# Preparación del pipeline para Deniz

La versión adaptada concentra las entradas en un Excel con dos hojas y las opciones de ejecución en `config.R`. La guía de uso para Deniz está en `README.md`, en inglés.

## Qué cambia Deniz

1. Copia `input/study_template.xlsx` como `input/my_study.xlsx`.
2. En `Measurements` introduce sus participantes, grupo, número de visita y concentraciones. Añade una columna por aminoácido.
3. En `Analytes` añade una fila por aminoácido, con el mismo nombre de columna, unidad, LoQ, priors de CVI/CVG/CVA, pesos y fuentes.
4. Cambia `INPUT_FILE` en `config.R` y comienza con `MODE <- "validate"`.
5. Revisa los recuentos y usa `MODE <- "fit"` para estimar VB.

El CVA se introduce como porcentaje: 5 significa 5 %. El peso analítico es una fracción, por ejemplo 0,10. El peso del prior de CVI tiene una parametrización distinta, explicada en la guía.

## Conversión de tus archivos

| Entrada original | Entrada adaptada |
|---|---|
| `subject` de resultados y participantes | `Measurements$subject` |
| `Sex` del archivo de participantes | `Measurements$group`, enlazado por `subject` |
| `Visita` | `Measurements$sample_order` |
| Columnas de esteroides | Mismos nombres, eliminando espacios finales |
| Hoja `LOQ` del ejemplo | `Analytes$loq`, en ng/mL |
| `Magnitude` | `Analytes$analyte` |
| `prior_cv_within` | `prior_cvi_pct` |
| `prior_cv_between` | `prior_cvg_pct` |
| `prior_cv_method` | `prior_cva_pct` |
| `weight_cv_within` | `weight_cvi` |
| `weight_cv_inter` | `weight_cvg` |
| `weight_cva` | `weight_cva` |

El ejemplo adaptado conserva los cuatro resultados originales, sus grupos, los 11 LoQ y los 11 conjuntos de hiperparámetros. Se eliminaron filas vacías y columnas logísticas que no intervienen en el modelo. La edad tampoco era una covariable del modelo original.

Las filas del Excel de resultados están en orden de análisis aleatorizado. El ejemplo no tiene `fecha_muestra`, por lo que el orden de visitas debe tomarse de `Visita`. La nueva versión lo exige explícitamente.

## Qué requiere comprobación antes de usar resultados

La estructura del modelo bayesiano se conserva, pero hay cambios que pueden afectar resultados: orden temporal correcto, valor crítico de Cochran corregido, y datos idénticos entre análisis principal y sensibilidad tras excluir tendencias. Están detallados en `CHANGELOG.md`.

El ejemplo de cuatro muestras no permite reproducir los resultados del manuscrito. Se incluye un ejemplo simulado de 72 muestras para comprobar el recorrido de ejecución. Sus valores, priors y LoQ son ficticios.

Se han comprobado la conversión celda a celda y la presentación de las plantillas. En el entorno de preparación no está instalado R/Stan: quedan pendientes las pruebas R incluidas, la compilación y el muestreo, y una comparación con los datos completos de atletas. No se presenta esta adaptación como una nueva versión estadísticamente validada.

Para GitHub, subir el contenido de la carpeta extraída. `README.md` será la página inicial. Se incluye `.gitignore` para excluir resultados, ajustes y nuevos archivos de entrada privados del seguimiento habitual con Git. Esta entrega no se ha publicado en un repositorio.
