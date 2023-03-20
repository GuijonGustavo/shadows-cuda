/* 
# Made by:
# Sergio Mendoza <sergio@mendozza.org>
# Milton Santibañez <msantibanez@astro.unam.mx>
# Gustavo Magallanes-Guijon <gustavo.magallanes.guijon@ciencias.unam.mx>
# Instituto de Astronomia UNAM
# Ciudad Universitaria
# Ciudad de Mexico
# Mexico
# Fri 21 Oct 2020 05:40:36 PM UTC
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <cuda_runtime.h>
#include "Error.h"
#include "Functions.h"
#include <time.h>
#include "GpuTimer.h"
#include "CpuTimer.h"
#include <gsl/gsl_integration.h>
#include <omp.h>


