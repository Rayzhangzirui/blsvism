#ifndef CFA_H
#define CFA_H

#include "Solute.h"
#include "util.h"
#include "globals.h"

//sequential initialize GB
void seqCFA(SoluteStruct &sol, double*** GB, GridData &grid);

double GB2outsidebox3D(const SoluteStruct& sol, GridData &grid);

#endif