/*!\file:  Vector.h
 * \brief wrapper to vector objects. The goal is to control which API (PETSc,Scalpack, Plapack?)
 * implements our underlying vector format.
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include <cstring>
#include "../../shared/Enum/Enum.h"
#include "../petsc/petscincludes.h"
#include "../issm/issmtoolkit.h"

#ifdef _HAVE_CODIPACK_
#include "../codipack/CoDiPackDebug.h"
#endif
/*}}}*/

enum vectortype { PetscVecType, IssmVecType };

template <class doubletype>
class Vector{

	public:

		int  type;
		#ifdef _HAVE_PETSC_
		PetscVec<doubletype>* pvector;
		#endif
		IssmVec<doubletype>* ivector;

		/*Vector constructors, destructors */
		Vector(){ /*{{{*/

			InitCheckAndSetType();
		}
		/*}}}*/
		Vector(int M,bool fromlocalsize=false){ /*{{{*/

			InitCheckAndSetType();

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector=new PetscVec<doubletype>(M,fromlocalsize);
				#endif
			}
			else this->ivector=new IssmVec<doubletype>(M,fromlocalsize);

		}
		/*}}}*/
		Vector(int m,int M){ /*{{{*/

			InitCheckAndSetType();

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
					this->pvector=new PetscVec<doubletype>(m,M);
				 #endif
			}
			else this->ivector=new IssmVec<doubletype>(m,M);
		}
		/*}}}*/
		Vector(doubletype* serial_vec,int M){ /*{{{*/

			ArrayDebugOutput("constructor from serial", M, serial_vec);
			InitCheckAndSetType();

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector=new PetscVec<doubletype>(serial_vec,M);
				#endif
			}
			else this->ivector=new IssmVec<doubletype>(serial_vec,M);
		}
		/*}}}*/
		~Vector(){ /*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				delete this->pvector;
				#endif
			}
			else delete this->ivector;
		}
		/*}}}*/
		#ifdef _HAVE_PETSC_
		Vector(PVec petsc_vector){ /*{{{*/

			this->type=PetscVecType;
			this->ivector=NULL;
			this->pvector=new PetscVec<doubletype>(petsc_vector);

		}
		/*}}}*/
		#endif
		void InitCheckAndSetType(void){ /*{{{*/

			#ifdef _HAVE_PETSC_
			pvector=NULL;
			#endif
			ivector=NULL;

			/*retrieve toolkittype: */
			char* toolkittype=ToolkitOptions::GetToolkitType();
			_assert_(toolkittype);

			/*set vector type: */
			if(strcmp(toolkittype,"petsc")==0){
				#ifdef _HAVE_PETSC_
				type=PetscVecType;
				#else
				_error_("cannot create petsc vector without PETSC compiled!");
				#endif
			}
			else if(strcmp(toolkittype,"issm")==0){
				/*let this choice stand:*/
				type=IssmVecType;
			}
			else{
				_error_("unknow toolkit type ");
			}

			/*Free resources: */
			xDelete<char>(toolkittype);
		}
		/*}}}*/

		/*Vector specific routines*/
		void Echo(void){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Echo();
				#endif
			}
			else this->ivector->Echo();

		}
		/*}}}*/
		void EchoDebug(std::string message){_assert_(this);/*{{{*/

			if(type==PetscVecType){
#ifdef _HAVE_PETSC_
				this->pvector->EchoDebug(message);
#endif
			}
			else this->ivector->EchoDebug(message);
		}
		/*}}}*/
		void Assemble(void){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Assemble();
				#endif
			}
			else this->ivector->Assemble();
			EchoDebug("assemble");
		}
		/*}}}*/
		void SetValues(int ssize, int* list, doubletype* values, InsMode mode){ _assert_(this);/*{{{*/
			ArrayDebugOutput("setValues in", ssize, values);

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->SetValues(ssize,list,values,mode);
				#endif
			}
			else this->ivector->SetValues(ssize,list,values,mode);

		}
		/*}}}*/
		void SetValue(int dof, doubletype value, InsMode mode){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->SetValue(dof,value,mode);
				#endif
			}
			else this->ivector->SetValue(dof,value,mode);

		}
		/*}}}*/
		void GetValue(doubletype* pvalue,int dof){_assert_(this);/*{{{*/

			EchoDebug("GetValue");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->GetValue(pvalue,dof);
				#endif
			}
			else this->ivector->GetValue(pvalue,dof);

			ArrayDebugOutput("getValue out", 1, pvalue);

		}
		/*}}}*/
		void GetSize(int* pM){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->GetSize(pM);
				#endif
			}
			else this->ivector->GetSize(pM);

		}
		/*}}}*/
		bool IsEmpty(void){/*{{{*/
			int M;

			_assert_(this);
			this->GetSize(&M);

			if(M==0)
				return true;
			else
				return false;
		}
		/*}}}*/
		void GetLocalSize(int* pM){_assert_(this);/*{{{*/

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->GetLocalSize(pM);
				#endif
			}
			else this->ivector->GetLocalSize(pM);

		}
		/*}}}*/
		void GetLocalVector(doubletype** pvector,int** pindices){_assert_(this);/*{{{*/

			EchoDebug("GetLocalVector");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->GetLocalVector(pvector,pindices);
				#endif
			}
			else this->ivector->GetLocalVector(pvector,pindices);

			int size;
			this->GetLocalSize(&size);
			ArrayDebugOutput("getLocalVector out", size, *pvector);

		}
		/*}}}*/
		Vector<doubletype>* Duplicate(void){_assert_(this);/*{{{*/

			Vector<doubletype>* output=NULL;

			output=new Vector<doubletype>();

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				output->pvector=this->pvector->Duplicate();
				#endif
			}
			else output->ivector=this->ivector->Duplicate();

			return output;
		} /*}}}*/
		void Set(doubletype value){_assert_(this);/*{{{*/

			ArrayDebugOutput("set in", 1, &value);

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Set(value);
				#endif
			}
			else this->ivector->Set(value);
			EchoDebug("Set");

		}
		/*}}}*/
		void AXPY(Vector* X, doubletype a){_assert_(this);/*{{{*/

			ArrayDebugOutput("AXPY a", 1, &a);
			X->EchoDebug("AXPY x");
			EchoDebug("AXPY y in");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->AXPY(X->pvector,a);
				#endif
			}
			else this->ivector->AXPY(X->ivector,a);
			EchoDebug("AXPY y out");

		}
		/*}}}*/
		void AYPX(Vector* X, doubletype a){_assert_(this);/*{{{*/

			ArrayDebugOutput("AyPX a", 1, &a);
			X->EchoDebug("AYPX x");
			EchoDebug("AYPX y in");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->AYPX(X->pvector,a);
				#endif
			}
			else this->ivector->AYPX(X->ivector,a);
			EchoDebug("AYPX y out");
		}
		/*}}}*/
		doubletype* ToMPISerial(void){/*{{{*/

			EchoDebug("ToMpiSerial");
			doubletype* vec_serial=NULL;

			_assert_(this);
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				vec_serial=this->pvector->ToMPISerial();
				#endif
			}
			else vec_serial=this->ivector->ToMPISerial();

			int size;
			this->GetLocalSize(&size);
			ArrayDebugOutput("ToMpiSerial out", size, vec_serial);

			return vec_serial;

		}
		/*}}}*/
		doubletype* ToMPISerial0(void){/*{{{*/

			EchoDebug("ToMpiSerial0");
			doubletype* vec_serial=NULL;

			_assert_(this);
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				vec_serial=this->pvector->ToMPISerial0();
				#else
				_error_("Cannot serialize PETSc Vec without PETSc");
				#endif
			}
			else vec_serial=this->ivector->ToMPISerial0();

			int size;
			this->GetLocalSize(&size);
			ArrayDebugOutput("ToMpiSerial0 out", size, vec_serial);

			return vec_serial;

		}
		/*}}}*/
		void Shift(doubletype shift){_assert_(this);/*{{{*/

			ArrayDebugOutput("shift in", 1, &shift);

			EchoDebug("Shift in");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Shift(shift);
				#endif
			}
			else this->ivector->Shift(shift);
			EchoDebug("Shift out");
		}
		/*}}}*/
		void Copy(Vector* to){_assert_(this);/*{{{*/

			EchoDebug("Copy in");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Copy(to->pvector);
				#endif
			}
			else this->ivector->Copy(to->ivector);

			to->EchoDebug("Copy out");

		}
		/*}}}*/
		doubletype Max(void){_assert_(this);/*{{{*/

			EchoDebug("Max");
			doubletype max=0;

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				max=this->pvector->Max();
				#endif
			}
			else _error_("operation not supported yet");

			ArrayDebugOutput("max out", 1, &max);
			return max;
		}
		/*}}}*/
		doubletype Norm(NormMode norm_type){_assert_(this);/*{{{*/

			doubletype norm=0;

			EchoDebug("Norm");

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				norm=this->pvector->Norm(norm_type);
				#endif
			}
			else norm=this->ivector->Norm(norm_type);

			ArrayDebugOutput("norm out", 1, &norm);
			return norm;
		}
		/*}}}*/
		void Scale(doubletype scale_factor){_assert_(this);/*{{{*/

			ArrayDebugOutput("scale s in", 1, &scale_factor);
			EchoDebug("Scale in");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Scale(scale_factor);
				#endif
			}
			else this->ivector->Scale(scale_factor);

			EchoDebug("Scale out");
		}
		/*}}}*/
		doubletype Dot(Vector* vector){_assert_(this);/*{{{*/

			EchoDebug("Dot 1");
			vector->EchoDebug("Dot 2");
			doubletype dot;

			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				dot=this->pvector->Dot(vector->pvector);
				#endif
			}
			else dot=this->ivector->Dot(vector->ivector);
			ArrayDebugOutput("dot out", 1, &dot);
			return dot;
		}
		/*}}}*/
		void PointwiseDivide(Vector* x,Vector* y){_assert_(this);/*{{{*/

			x->EchoDebug("PointwiseDevide x");
			y->EchoDebug("PointwiseDevide y");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->PointwiseDivide(x->pvector,y->pvector);
				#endif
			}
			else this->ivector->PointwiseDivide(x->ivector,y->ivector);
			EchoDebug("PointwiseDevide r");
		}
		/*}}}*/
		void PointwiseMult(Vector* x,Vector* y){_assert_(this);/*{{{*/

			x->EchoDebug("PointwiseMult x");
			y->EchoDebug("PointwiseMult y");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->PointwiseMult(x->pvector,y->pvector);
				#endif
			}
			else this->ivector->PointwiseMult(x->ivector,y->ivector);
			EchoDebug("PointwiseMult r");
		}
		/*}}}*/
		void Pow(doubletype scale_factor){_assert_(this);/*{{{*/

			ArrayDebugOutput("pow s in", 1, &scale_factor);
			EchoDebug("pow in");
			if(type==PetscVecType){
				#ifdef _HAVE_PETSC_
				this->pvector->Pow(scale_factor);
				#endif
			}
			else this->ivector->Pow(scale_factor);
			EchoDebug("pow out");
		}
		/*}}}*/
void Sum(doubletype* pvalue){ /*{{{*/
	_assert_(this);/*{{{*/

	EchoDebug("sum");
	if(type==PetscVecType){
		#ifdef _HAVE_PETSC_
		this->pvector->Sum(pvalue);
		#endif
	}
	else this->ivector->Sum(pvalue);

	ArrayDebugOutput("sum out", 1, pvalue);
}
/*}}}*/
}; /*}}}*/
#endif //#ifndef _VECTOR_H_
