# Algoritmo de Peterson-Gorenstein-Zierler para códigos cíclicos sesgados

> Trabajo de Fin de Grado. Doble Grado en Ingeniería Informática y Matemáticas. Universidad de Granada

> Realizado por José María Martín Luque y tutorizado por Gabriel Navarro Garulo


El objetivo principal de este trabajo es estudiar, presentar e implementar el algoritmo de Peterson-Gorenstein-Zierler para códigos cíclicos sesgados.
Por tanto, nos encargamos de exponer todo el conocimiento necesario para abordar su estudio.
Para ello primero nos centraremos en los fundamentos necesarios de álgebra y teoría de códigos, especialmente en los códigos lineales y los códigos cíclicos.
Estudiaremos además la versión original del algoritmo PGZ para ayudarnos en la compresión del algoritmo para códigos cíclicos sesgados.
Posteriormente abordaremos los anillos de polinomios de Ore, que constituyen la base de los códigos cíclicos sesgados.
Una vez explicada esta clase de códigos, estaremos en disposición de estudiar e implementar el algoritmo que es el objeto de nuestro trabajo:  Peterson-Gorenstein-Zierler para códigos cíclicos sesgados.

Palabras clave: teoría de códigos, códigos cíclicos, polinomios torcidos, códigos cíclicos sesgados, Peterson-Gorenstein-Zierler

## Texto

Para compilar el texto es necesario utilizar `lualatex`. 
Se necesitan además las fuentes *EBGaramond*, *Garamond Math*, *Go Mono* y *Open Sans*.

## Código

El código se ha desarrollado con SageMath 9.0.
Para usar las clases desarrolladas simplemente es necesario cargar los archivos *.sage* de la carpeta *code* mediante la orden `load()` de Sage.