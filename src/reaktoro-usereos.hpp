#pragma once

#include <cstdint>

#define MF_API extern "C" __declspec(dllexport)

// Exported functions for USEREOS

MF_API void ReadConfigurationFile(char *filename, int *ierr, char *title);

MF_API void GetDimensions(int *nComponents, int *nPhaseMax, int *nAux);

MF_API void GetGlobalParameters(char *cmpNames, double *molWeights, char *phNames, char *auxNames, char *auxUnits,
                                int8_t *opt);

MF_API void PhaseEquilibrium(double *pres, double *temp, double *z, int *nPhase, double *props, int8_t *phaseId,
                             double *auxArray, int8_t *mode);
