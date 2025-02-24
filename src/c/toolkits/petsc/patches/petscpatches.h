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
#include </home/msagebaum/Kaiserslautern/Programms/adjoint_petsc/include/adjoint_petsc/vec.h>
#endif

class Parameters;

using PVec = adjoint_petsc::ADVec;

PVec NewVec(int size,ISSM_MPI_Comm comm,bool fromlocalsize=false);
Mat NewMat(int M,int N,ISSM_MPI_Comm comm);
Mat NewMat(int M,int N,double sparsity,ISSM_MPI_Comm comm);
Mat NewMat(int M,int N,int connectivity,int numberofdofspernode, ISSM_MPI_Comm comm);

int VecToMPISerial(IssmDouble** pgathered_vector, PVec vector,ISSM_MPI_Comm comm,bool broadcast=true);
void MatFree(Mat* pmat);
void ISFree(IS* pis);
void VecFree(PVec* pvec);
void KSPFree(KSP* pksp);
int MatPartition(Mat* poutmatrix,Mat matrixA,double* row_partition_vector,int row_partition_vector_size, double* col_partition_vector,int col_partition_vector_size);
void PetscOptionsDetermineSolverType(int* psolver_type);
void MatToMPISerial(double** poutmatrix,Mat matrix,ISSM_MPI_Comm comm,bool broadcast=true);
Vec  SerialToVec(double* vector,int vector_size);
InsertMode ISSMToPetscInsertMode(InsMode mode);
NormType ISSMToPetscNormMode(NormMode mode);
MatType ISSMToPetscMatrixType(MatrixType type);

#endif
