#include "bddSystem.h"

/**Function*************************************************************

  Synopsis    [Alloc a new tensor.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

Tensor* BDDSystem::newTensor(int r, int rank)
{
	return new Tensor(r, rank);
}

/**Function*************************************************************

  Synopsis    [Delete a tensor.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void BDDSystem::deleteTensor(Tensor* tensor)
{
	for (int i = 0; i < _w; i++)
		for (int j = 0, end_j = tensor->_r; j < end_j; j++)
			Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
}

static void bitVectorPlus1(int length, int *reg)
{
    int one = 1, carry = 0;
    for (int i = 0; i < length; ++i)
    {
        carry = reg[i] & one;
        reg[i] = reg[i] ^ one;

        if (carry == 0)
            break;
    }
}

/**Function*************************************************************

  Synopsis    [Print a tensor.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void BDDSystem::printTensor(Tensor* tensor) const
{
	const double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899;
	const double oneRoot2 = 1 / sqrt(2);
	const double hFactor = pow(oneRoot2, tensor->_k);
	double re = 0, im = 0;
	unsigned long long nEntries = pow(2, tensor->_rank);
	int *assign = new int[tensor->_rank];
	for (int i = 0; i < tensor->_rank; ++i)
		assign[i] = 0;
 
	std::cout << "r: " << tensor->_r << '\n';
	std::cout << "k: " << tensor->_k << '\n';
	std::cout << "rank: " << tensor->_rank << '\n';
	std::cout << '[';
	for (unsigned long long i = 0; i < nEntries; ++i)
	{
		re = 0;
		im = 0;
		
		for (int j = 0; j < _w; j++) // compute every complex value
		{
			long long intValue = 0;
			for (int h = 0; h < tensor->_r; h++) // compute every integer
			{
				DdNode* tmp = Cudd_Eval(_ddManager, tensor->_allBDD[j][h], assign);
				Cudd_Ref(tmp);
				int oneEntry = !(Cudd_IsComplement(tmp));
				Cudd_RecursiveDeref(_ddManager, tmp);
				if (h == tensor->_r - 1)
					intValue -= oneEntry * pow(2, h);
				else
					intValue += oneEntry * pow(2, h);
			}
			/* translate to re and im */
			re += intValue * cos((double) (_w - j - 1)/_w * PI);
			im += intValue * sin((double) (_w - j - 1)/_w * PI);
		}
		re *= hFactor;
		im *= hFactor;

		if ((re == 0)&&(im == 0))
			std::cout << "\"0\"";
		else if (re == 0)
			std::cout << "\"" + std::to_string(im) + "i\"";
		else if (im == 0)
			std::cout << "\"" + std::to_string(re) + "\"";
		else
		{
			if (im < 0)
				std::cout << "\"" + std::to_string(re) + std::to_string(im) + "i\"";
			else
				std::cout << "\"" + std::to_string(re) + "+" + std::to_string(im) + "i\"";
		}
		if (i != nEntries - 1)
			std::cout << ", ";
		bitVectorPlus1(tensor->_rank, assign);
	}
	std::cout << "]\n";
}

/**Function*************************************************************

  Synopsis    [Increase the bitwidth of the given allBDD by 1]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::incBDDsBitWidth(std::vector<std::deque<DdNode*>> &allBDD)
{
	for(int j = 0; j < _w; ++j)
	{
		Cudd_Ref(allBDD[j].back());
		allBDD[j].push_back(allBDD[j].back());
	}
}

/**Function*************************************************************

  Synopsis    [Check if the LSB of the tensor is zero or not.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

bool BDDSystem::isTensorLSBZero(Tensor *tensor) const
{
	for(int i = 0; i < _w; ++i)
	{
		if(tensor->_allBDD[i].front() != Cudd_Not(Cudd_ReadOne(_ddManager)))
		{
			return false;
		}
	}
	return true;
}

/**Function*************************************************************

  Synopsis    [Drop the LSB of the given tensor]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::dropTensorLSB(Tensor *tensor)
{
	for(int j = 0; j < _w; ++j)
	{
		Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[j].front());
		tensor->_allBDD[j].pop_front();
	}
	--tensor->_r;
	tensor->_k -= 2;
}

/**Function*************************************************************

  Synopsis    [Drop the bits of the tensor when the tensor overflows]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/

void BDDSystem::dropTensorBits(Tensor *tensor)
{
	if(isTensorLSBZero(tensor) || _bitWidthMode == BitWidthMode::DropLSB)
		dropTensorLSB(tensor);
}

/**Function*************************************************************

  Synopsis    [Increase the k of a given tensor by 1.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::increaseTensorKByOne(Tensor* tensor)
{
	DdNode *carry, *g, *d, *tmp, *tmp2, *tmp3;
	bool isOverflow = false;

	auto copyAllBDD = tensor->_allBDD;
    for (int i = 0; i < _w; ++i)
         for (int j = 0; j < tensor->_r; ++j)
            Cudd_Ref(copyAllBDD[i][j]);
	
	++tensor->_k;

	for(int i = 0; i < _w; ++i)
	{
		switch(i)
		{
			case 0: carry = Cudd_ReadOne(_ddManager); break;
			case 1: carry = Cudd_Not(Cudd_ReadOne(_ddManager)); break;
			case 2: carry = Cudd_Not(Cudd_ReadOne(_ddManager)); break;
			case 3: carry = Cudd_ReadOne(_ddManager); break;
		}
		Cudd_Ref(carry);

		for(int j = 0; j < tensor->_r; ++j)
		{
			switch(i)
			{
				case 0:
					g = copyAllBDD[1][j];
					d = Cudd_Not(copyAllBDD[3][j]);
					break;
				case 1:
					g = copyAllBDD[2][j];
					d = copyAllBDD[0][j];
					break;
				case 2:
					g = copyAllBDD[1][j];
					d = copyAllBDD[3][j];
					break;
				case 3:
					g = copyAllBDD[2][j];
					d = Cudd_Not(copyAllBDD[0][j]);
					break;
			}

			if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(g, d, carry))
			{
				incBDDsBitWidth(tensor->_allBDD);
				incBDDsBitWidth(copyAllBDD);
				++tensor->_r;
				isOverflow = true;
			}

			tmp = Cudd_bddXor(_ddManager, g, d);
			Cudd_Ref(tmp);
			tmp2 = Cudd_bddXor(_ddManager, tmp, carry);
			Cudd_Ref(tmp2);
			Cudd_RecursiveDeref(_ddManager, tmp);
			Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
			tensor->_allBDD[i][j] = tmp2;

			// Carry
			if (j == tensor->_r - 1)
				Cudd_RecursiveDeref(_ddManager, carry);
			else
			{
				tmp = Cudd_bddAnd(_ddManager, g, d);
				Cudd_Ref(tmp);
				tmp2 = Cudd_bddOr(_ddManager, g, d);
				Cudd_Ref(tmp2);
				tmp3 = Cudd_bddAnd(_ddManager, tmp2, carry);
				Cudd_Ref(tmp3);
				Cudd_RecursiveDeref(_ddManager, tmp2);
				Cudd_RecursiveDeref(_ddManager, carry);
				carry = Cudd_bddOr(_ddManager, tmp, tmp3);
				Cudd_Ref(carry);
				Cudd_RecursiveDeref(_ddManager, tmp);
				Cudd_RecursiveDeref(_ddManager, tmp3);
			}
		}
	}

    for (int i = 0; i < _w; ++i)
        for (int j = 0; j < tensor->_r; ++j)
            Cudd_RecursiveDeref(_ddManager, copyAllBDD[i][j]);
}

/**Function*************************************************************

  Synopsis    [Detect overflow in integer vectors.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
bool BDDSystem::checkAdderOverflow(DdNode *g, DdNode *h, DdNode *carryIn) const
{
    DdNode *tmp, *dd1, *dd2;
    bool overflow;

	// overflow = (g XNOR h) AND (g XOR carryIn)
    dd1 = Cudd_bddXor(_ddManager, g, carryIn);
    Cudd_Ref(dd1);

    dd2 = Cudd_bddXnor(_ddManager, g, h);
    Cudd_Ref(dd2);

    tmp = Cudd_bddAnd(_ddManager, dd1, dd2);
    Cudd_Ref(tmp);
    Cudd_RecursiveDeref(_ddManager, dd1);
    Cudd_RecursiveDeref(_ddManager, dd2);

    if (tmp != Cudd_Not(Cudd_ReadOne(_ddManager)))
        overflow = true;
    else
        overflow = false;
    Cudd_RecursiveDeref(_ddManager, tmp);

    return overflow;
}

/**Function*************************************************************

  Synopsis    [Update max #nodes.]

  Description []

  SideEffects []

  SeeAlso     []

***********************************************************************/
void BDDSystem::updateMaxNodeCount()
{
    _maxNodeCount = std::max(_maxNodeCount, static_cast<unsigned long>(Cudd_ReadNodeCount(_ddManager)));
}
