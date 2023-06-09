#include "bddSystem.h"

/**Function*************************************************************

  Synopsis    [Get the BDD characteriza the sparsity of the tensor.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
DdNode* BDDSystem::sparsityDD(Tensor *tensor) const
{
	DdNode *dd, *tmp;
	dd = Cudd_Not(Cudd_ReadOne(_ddManager));
	Cudd_Ref(dd);

	for(int i = 0; i < _w; ++i)
	{
		for(int j = 0, end_j = tensor->_r; j < end_j; ++j)
		{
			tmp = Cudd_bddOr(_ddManager, dd, tensor->_allBDD[i][j]);
			Cudd_Ref(tmp);
			Cudd_RecursiveDeref(_ddManager, dd);
			dd = tmp;
		}
	}

	return dd;
}

/**Function*************************************************************

  Synopsis    [Get the sparsity value of the tensor.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
double BDDSystem::sparsity(Tensor *tensor) const
{
	DdNode* dd = sparsityDD(tensor);
	double sp = 1 - Cudd_CountMinterm(_ddManager, dd, tensor->_rank) / pow(2, tensor->_rank);
	Cudd_RecursiveDeref(_ddManager, dd);
	return sp;
}

/**Function*************************************************************

  Synopsis    [Equivalence checking two tensors.]

  Description []

  SideEffects [The zero LSBs of the two tensor are dropped.]

  SeeAlso     []

***********************************************************************/
bool BDDSystem::eqCheckTwoTensor(Tensor *tensor1, Tensor *tensor2)
{
	if((tensor1->_k % 2) != (tensor2->_k % 2))
		increaseTensorKByOne(tensor2);

	while(isTensorLSBZero(tensor1))
		dropTensorLSB(tensor1);

	while(isTensorLSBZero(tensor2))
		dropTensorLSB(tensor2);

	if(tensor1->_r > tensor2->_r)
		std::swap(tensor1, tensor2);

	for(int i = 0; i < _w; ++i)
	{
		for(int j = 0; j < tensor1->_r; ++j)
		{
			if(tensor1->_allBDD[i][j] != tensor2->_allBDD[i][j])
				return false;
		}

		for(int j = tensor1->_r; j < tensor2->_r; ++j)
		{
			if(tensor2->_allBDD[i][j] != tensor2->_allBDD[i][tensor1->_r-1])
				return false;
		}
	}

	return true;
}
