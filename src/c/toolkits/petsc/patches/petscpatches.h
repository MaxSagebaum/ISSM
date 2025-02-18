/*\file petscpatches.h
 * \brief: our own patches for petsc use
 */

#ifndef _PETSC_PATCHES_H_
#define _PETSC_PATCHES_H_

#include <petscksp.h>

#include "./SolverEnum.h"
#include "../../toolkitsenums.h"
#include "../../../shared/io/Comm/IssmComm.h"

#if _HAVE_CODIPACK_
#include <adjoint_petsc/vec.h>
#include <adjoint_petsc/mat.h>
#include <adjoint_petsc/ksp.h>
#endif

class Parameters;

using PVec = adjoint_petsc::ADVec;
using PMat = adjoint_petsc::ADMat;
using PKSP = adjoint_petsc::ADKSP;
using PArray = adjoint_petsc::WrapperArray;

PVec NewVec(int size,ISSM_MPI_Comm comm,bool fromlocalsize=false);
PMat NewMat(int M,int N,ISSM_MPI_Comm comm);
PMat NewMat(int M,int N,IssmPDouble sparsity,ISSM_MPI_Comm comm);
PMat NewMat(int M,int N,int connectivity,int numberofdofspernode, ISSM_MPI_Comm comm);

int VecToMPISerial(IssmDouble** pgathered_vector, PVec vector,ISSM_MPI_Comm comm,bool broadcast=true);
void MatFree(PMat* pmat);
void ISFree(IS* pis);
void VecFree(PVec* pvec);
void KSPFree(PKSP* pksp);
int MatPartition(PMat* poutmatrix,PMat matrixA,IssmDouble* row_partition_vector,int row_partition_vector_size, IssmDouble* col_partition_vector,int col_partition_vector_size);
void PetscOptionsDetermineSolverType(int* psolver_type);
void MatToMPISerial(IssmDouble** poutmatrix,PMat matrix,ISSM_MPI_Comm comm,bool broadcast=true);
PVec  SerialToVec(IssmDouble* vector,int vector_size);
InsertMode ISSMToPetscInsertMode(InsMode mode);
NormType ISSMToPetscNormMode(NormMode mode);
MatType ISSMToPetscMatrixType(MatrixType type);

#endif
