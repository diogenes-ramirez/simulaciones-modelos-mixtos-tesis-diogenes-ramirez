### Descripción de la Simulación

Este proyecto presenta una serie de simulaciones realizadas para evaluar el desempeño de un Criterio de Información Bayesiano (BIC) modificado, que incorpora el logaritmo del determinante de la matriz Hessiana de los parámetros (denominado BIC HES). El objetivo principal es analizar cómo esta modificación afecta la selección de modelos en diferentes contextos y compararla con otros criterios como BIC tradicional, AIC, AIC corregido (AICc), CAIF e ICOMP.

### 1. Simulaciones para el Modelo Lineal

#### 1.1. Procedimiento
- **Generación de Variables Regresoras**: Se generaron k variables regresoras independientes (k = 2, 4, 6, 8, 10), con una distribución normal (media = 10, desviación estándar = 1).
- **Generación de Parámetros**: Los parámetros correspondientes a las variables regresoras siguen una distribución normal con media = 10 y desviación estándar = 2.
- **Errores**: Errores distribuidos normalmente con media = 0 y desviación estándar = 1.
- **Generación de Variables de Ruido**: Se generaron variables regresoras adicionales (denominadas de ruido) con media = 3 y desviación estándar = 3, agregadas a los modelos para evaluar el impacto en la selección del modelo verdadero.
- **Cálculo de Criterios**: Se calcularon los criterios BIC, AIC, AICc y BIC HES para cada configuración.

#### 1.2. Resultados
- **Selección del Modelo Verdadero**: El criterio BIC HES mostró una mayor tasa de selección del modelo verdadero, especialmente en muestras pequeñas, comparado con BIC, AIC y AICc.
- **Comparaciones Directas**: Se realizaron comparaciones entre BIC HES y BIC, así como entre BIC HES y AIC, mostrando que BIC HES tiene una ventaja en la selección del modelo correcto.

### 2. Simulaciones para el Modelo Lineal Mixto

#### 2.1. Procedimiento
- **Generación de Datos**: Se generaron datos simulados para un modelo lineal mixto, incorporando efectos fijos y aleatorios, utilizando distribuciones normales para los parámetros y errores.
- **Tamaño de Muestra y Grupos**: Las simulaciones variaron el tamaño de muestra (n = 20, 50, 100, 500) y el número de grupos (N = 5, 10, 20).
- **Variables de Ruido**: Se incluyeron variables de ruido tanto para los efectos fijos como aleatorios.

#### 2.2. Resultados
- **Comparación de Criterios**: BIC HES demostró ser más efectivo en la selección del modelo verdadero en comparación con otros criterios, especialmente cuando se añadieron pocas variables de ruido.

### 3. Simulaciones en Modelos Binomial y Poisson

#### 3.1. Procedimiento
- **Generación de Datos**: Se generaron datos bajo diferentes configuraciones para modelos binomiales y Poisson, utilizando un enfoque de efectos aleatorios.
- **Evaluación de Criterios**: Se evaluaron múltiples criterios de selección de modelos, incluyendo BIC, AIC, AICc, BIC HES, CAIF e ICOMP.

#### 3.2. Resultados
- **Eficiencia de BIC HES**: BIC HES y su variante BIC HES SP fueron los más efectivos en la selección del modelo correcto, especialmente en presencia de efectos aleatorios significativos.

### 4. Prueba de Hipótesis

Se realizó una prueba de hipótesis para comparar la exactitud del BIC HES con otros criterios de información, utilizando una prueba exacta de Fisher. Los resultados sugieren que BIC HES es superior en muchos casos, especialmente en escenarios de muestras pequeñas y modelos mixtos.

### 5. Aplicaciones

#### 5.1. Modelo Lineal
Se aplicó el BIC HES en un modelo de regresión lineal utilizando datos reales del conjunto de datos BudgetUK, comparando su desempeño con otros criterios.

#### 5.2. Modelo Lineal Mixto
Dos aplicaciones de modelos lineales mixtos demostraron la robustez de BIC HES para la selección de modelos en contextos ecológicos y de comportamiento animal.

---

Este repositorio contiene los scripts y datos utilizados para realizar las simulaciones descritas, junto con las instrucciones necesarias para replicar los resultados y comparar diferentes criterios de selección de modelos en modelos lineales y mixtos.
