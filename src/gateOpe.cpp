#include "bddSystem.h"

/*
Apply gate on the quantum data.
*/

void BDDSystem::Toffoli(Tensor* tensor, const std::vector<int> &qubits)
{
    DdNode *term1, *term2, *term3, *g, *tmp;

	const int targ = qubits.back();
	std::vector<int> cont(qubits.begin(), qubits.end()-1);

	assert(targ < tensor->_rank && targ >= 0);
	for(const auto ele: cont)
		assert(ele >= 0 && ele < tensor->_rank && ele != targ);

    g = Cudd_ReadOne(_ddManager);
    Cudd_Ref(g);
    for (int h = cont.size() - 1; h >= 0; h--)
    {
        tmp = Cudd_bddAnd(_ddManager, Cudd_bddIthVar(_ddManager, cont[h]), g);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, g);
        g = tmp;
    }

    for (int i = 0; i < _w; i++)
    {
        for (int j = 0; j < tensor->_r; j++)
        {
            // term1
            term1 = Cudd_ReadOne(_ddManager);
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, tensor->_allBDD[i][j], term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            tmp = Cudd_bddAnd(_ddManager, Cudd_Not(g), term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;

            // term2
            term2 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(term2);

            tmp = Cudd_Cofactor(_ddManager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            // term3
            term3 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(term3);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);

            tmp = Cudd_Cofactor(_ddManager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            // OR
            tensor->_allBDD[i][j] = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(tensor->_allBDD[i][j]);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
            tmp = Cudd_bddOr(_ddManager, term3, tensor->_allBDD[i][j]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = tmp;
        }
    }
    Cudd_RecursiveDeref(_ddManager, g);
}

void BDDSystem::Fredkin(Tensor* tensor, const std::vector<int> &qubits)
{
    DdNode *term1, *term2, *term3, *g, *tmp, *tmp0;

	assert(qubits.size() >= 2);
	const int swapA = qubits[qubits.size()-1];
	const int swapB = qubits[qubits.size()-2];
	std::vector<int> cont(qubits.begin(), qubits.end()-2);

	assert(swapA >= 0 && swapA < tensor->_rank);
	assert(swapB >= 0 && swapB < tensor->_rank);
	assert(swapA != swapB);
	for(const auto ele: cont)
		assert(ele >= 0 && ele < tensor->_rank && ele != swapA && ele != swapB);

    g = Cudd_ReadOne(_ddManager);
    Cudd_Ref(g);
    for (int h = cont.size() - 1; h >= 0; h--)
    {
        tmp = Cudd_bddAnd(_ddManager, Cudd_bddIthVar(_ddManager, cont[h]), g);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, g);
        g = tmp;
    }

    for (int i = 0; i < _w; i++) // F = allBDD[i][j]
    {
        for (int j = 0; j < tensor->_r; j++)
        {
            // term1
            term1 = Cudd_ReadOne(_ddManager);
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, tensor->_allBDD[i][j], term1);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            tmp = Cudd_bddXor(_ddManager, Cudd_bddIthVar(_ddManager, swapA), Cudd_bddIthVar(_ddManager, swapB));
            Cudd_Ref(tmp);
            tmp0 = Cudd_Not(Cudd_bddAnd(_ddManager, g, tmp));
            Cudd_Ref(tmp0);
            Cudd_RecursiveDeref(_ddManager, tmp);
            tmp = Cudd_bddAnd(_ddManager, term1, tmp0);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, tmp0);
            term1 = tmp;

            // term2
            term2 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], g);
            Cudd_Ref(term2);

            tmp = Cudd_Cofactor(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, swapA)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_Cofactor(_ddManager, term2, Cudd_bddIthVar(_ddManager, swapB));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_bddIthVar(_ddManager, swapA));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, swapB)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            // term3
            term3 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], g);
            Cudd_Ref(term3);

            tmp = Cudd_Cofactor(_ddManager, term3, Cudd_bddIthVar(_ddManager, swapA));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_Cofactor(_ddManager, term3, Cudd_Not(Cudd_bddIthVar(_ddManager, swapB)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, g);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, Cudd_Not(Cudd_bddIthVar(_ddManager, swapA)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            tmp = Cudd_bddAnd(_ddManager, term3, Cudd_bddIthVar(_ddManager, swapB));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            term3 = tmp;

            // OR
            tensor->_allBDD[i][j] = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(tensor->_allBDD[i][j]);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
            tmp = Cudd_bddOr(_ddManager, term3, tensor->_allBDD[i][j]);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term3);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = tmp;
        }
    }
    Cudd_RecursiveDeref(_ddManager, g);
}

void BDDSystem::Hadamard(Tensor *tensor, const std::vector<int> &qubits)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back(); 
	
	assert(targ >= 0 && targ < tensor->_rank);

    tensor->_k += 1;

    DdNode *g, *d, *c, *tmp, *term1, *term2;

    int isOverflow = 0;

    for (int i = 0; i < _w; i++) 
    {
        c = Cudd_ReadOne(_ddManager); 
        Cudd_Ref(c);
        tmp = Cudd_bddAnd(_ddManager, c, Cudd_bddIthVar(_ddManager, targ));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, c);
        c = tmp;
        for (int j = 0; j < tensor->_r; j++)
        {
            // G
            g = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(g);

            // D
            term1 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, term1, Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            term2 = Cudd_Not(tensor->_allBDD[i][j]);
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            d = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(d);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

            // Detect overflow
            if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(g, d, c))
			{
				incBDDsBitWidth(tensor->_allBDD); 
				++tensor->_r;
				isOverflow = 1;
			}

            // Sum
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, g, d);
            Cudd_Ref(tensor->_allBDD[i][j]);
            tmp = Cudd_bddXor(_ddManager, tensor->_allBDD[i][j], c);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = tmp;

            // Carry
            if (j == tensor->_r - 1)
            {
                Cudd_RecursiveDeref(_ddManager, c);
                Cudd_RecursiveDeref(_ddManager, g);
                Cudd_RecursiveDeref(_ddManager, d);
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, g, d);
                Cudd_Ref(term1);
                term2 = Cudd_bddOr(_ddManager, g, d);
                Cudd_Ref(term2);
                Cudd_RecursiveDeref(_ddManager, g);
                Cudd_RecursiveDeref(_ddManager, d);
                tmp = Cudd_bddAnd(_ddManager, term2, c);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, term2);
                Cudd_RecursiveDeref(_ddManager, c);
                term2 = tmp;
                c = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(c);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::Rx_pi_2(Tensor *tensor, const std::vector<int> &qubits)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back();
	assert(targ >= 0 && targ < tensor->_rank);

    tensor->_k += 1;

    int nshift = _w / 2;
    int isOverflow = 0;

    DdNode *d, *c, *tmp, *term1, *term2;
	auto copy = tensor->_allBDD;
    for (int i = 0; i < _w; i++)
         for (int j = 0; j < tensor->_r; j++)
            Cudd_Ref(copy[i][j]);

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if (i < nshift)
            c = Cudd_ReadOne(_ddManager);
        else
            c = Cudd_Not(Cudd_ReadOne(_ddManager));
        Cudd_Ref(c);
        for (int j = 0; j < tensor->_r; j++)
        {
            // D
            term1 = Cudd_Cofactor(_ddManager, copy[(i + nshift) % _w][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, term1, Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            term2 = Cudd_Cofactor(_ddManager, copy[(i + nshift) % _w][j], Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            if (i < nshift) d = Cudd_Not(Cudd_bddOr(_ddManager, term1, term2));
            else d = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(d);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
            // Detect overflow
            if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(copy[i][j], d, c))
			{
				incBDDsBitWidth(tensor->_allBDD);
				incBDDsBitWidth(copy);
				++tensor->_r;
				isOverflow = 1;
			}
            // Sum
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, copy[i][j], d);
            Cudd_Ref(tensor->_allBDD[i][j]);
            tmp = Cudd_bddXor(_ddManager, tensor->_allBDD[i][j], c);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = tmp;
            // Carry
            if (j == tensor->_r - 1)
            {
                Cudd_RecursiveDeref(_ddManager, c);
                Cudd_RecursiveDeref(_ddManager, d);
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], d);
                Cudd_Ref(term1);
                term2 = Cudd_bddOr(_ddManager, copy[i][j], d);
                Cudd_Ref(term2);
                Cudd_RecursiveDeref(_ddManager, d);
                tmp = Cudd_bddAnd(_ddManager, term2, c);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, term2);
                Cudd_RecursiveDeref(_ddManager, c);
                term2 = tmp;
                c = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(c);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }

    for (int i = 0; i < _w; i++)
		for (int j = 0; j < tensor->_r; j++)
			Cudd_RecursiveDeref(_ddManager, copy[i][j]);

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::Rx_pi_2_dagger(Tensor *tensor, const std::vector<int> &qubits)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back();

	assert(targ >= 0 && targ < tensor->_rank);

    tensor->_k += 1;

    int nshift = _w / 2;
    int isOverflow = 0;

    DdNode *d, *c, *tmp, *term1, *term2;

	auto copy = tensor->_allBDD;
    for (int i = 0; i < _w; i++)
         for (int j = 0; j < tensor->_r; j++)
            Cudd_Ref(copy[i][j]);

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if (i < nshift)
            c = Cudd_Not(Cudd_ReadOne(_ddManager));
        else
            c = Cudd_ReadOne(_ddManager);
        Cudd_Ref(c);
        for (int j = 0; j < tensor->_r; j++)
        {
            // D
            term1 = Cudd_Cofactor(_ddManager, copy[(i + nshift) % _w][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(term1);
            tmp = Cudd_bddAnd(_ddManager, term1, Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;
            term2 = Cudd_Cofactor(_ddManager, copy[(i + nshift) % _w][j], Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
			if(i < nshift) d = Cudd_bddOr(_ddManager, term1, term2);
			else d = Cudd_Not(Cudd_bddOr(_ddManager, term1, term2));
            Cudd_Ref(d);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
            // Detect overflow
            if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(copy[i][j], d, c))
			{
				incBDDsBitWidth(tensor->_allBDD);
				incBDDsBitWidth(copy);
				++tensor->_r;
				isOverflow = 1;
			}
            // Sum
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, copy[i][j], d);
            Cudd_Ref(tensor->_allBDD[i][j]);
            tmp = Cudd_bddXor(_ddManager, tensor->_allBDD[i][j], c);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = tmp;
            // Carry
            if (j == tensor->_r - 1)
            {
                Cudd_RecursiveDeref(_ddManager, c);
                Cudd_RecursiveDeref(_ddManager, d);
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], d);
                Cudd_Ref(term1);
                term2 = Cudd_bddOr(_ddManager, copy[i][j], d);
                Cudd_Ref(term2);
                Cudd_RecursiveDeref(_ddManager, d);
                tmp = Cudd_bddAnd(_ddManager, term2, c);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, term2);
                Cudd_RecursiveDeref(_ddManager, c);
                term2 = tmp;
                c = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(c);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }
	
	for (int i = 0; i < _w; i++)
        for (int j = 0; j < tensor->_r; j++)
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::Ry_pi_2(Tensor *tensor, const std::vector<int> &qubits, const bool fTranspose)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back();

	assert(targ >= 0 && targ < tensor->_rank);

    tensor->_k += 1;

    int isOverflow = 0;

    DdNode *g, *d, *c, *tmp, *term1, *term2, *var;

    if (fTranspose) var = Cudd_Not(Cudd_bddIthVar(_ddManager, targ));
    else var = Cudd_bddIthVar(_ddManager, targ);

    for (int i = 0; i < _w; i++)
    {
        // Init C
        c = Cudd_ReadOne(_ddManager); 
        Cudd_Ref(c);
        tmp = Cudd_bddAnd(_ddManager, c, Cudd_Not(var));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, c);
        c = tmp;
        for (int j = 0; j < tensor->_r; j++)
        {
            // G
            g = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_Not(var));
            Cudd_Ref(g);
            // D
            term1 = Cudd_bddAnd(_ddManager, tensor->_allBDD[i][j], var);
            Cudd_Ref(term1);
            term2 = Cudd_Not(Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], var));
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(var));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            d = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(d);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

            // Detect overflow
            if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(g, d, c))
			{
				incBDDsBitWidth(tensor->_allBDD);
				++tensor->_r;
				isOverflow = 1;
			}
            // Sum
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, g, d);
            Cudd_Ref(tensor->_allBDD[i][j]);
            tmp = Cudd_bddXor(_ddManager, tensor->_allBDD[i][j], c);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = tmp;
            // Carry
            if (j == tensor->_r - 1)
            {
                Cudd_RecursiveDeref(_ddManager, c);
                Cudd_RecursiveDeref(_ddManager, g);
                Cudd_RecursiveDeref(_ddManager, d);
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, g, d);
                Cudd_Ref(term1);
                term2 = Cudd_bddOr(_ddManager, g, d);
                Cudd_Ref(term2);
                Cudd_RecursiveDeref(_ddManager, g);
                Cudd_RecursiveDeref(_ddManager, d);
                tmp = Cudd_bddAnd(_ddManager, term2, c);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, term2);
                Cudd_RecursiveDeref(_ddManager, c);
                term2 = tmp;
                c = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(c);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::Ry_pi_2_dagger(Tensor *tensor, const std::vector<int> &qubits, const bool fTranspose)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back();

	assert(targ >= 0 && targ < tensor->_rank);

    tensor->_k += 1;

    int isOverflow = 0;

    DdNode *g, *d, *c, *tmp, *term1, *term2, *var;

	if(fTranspose) var = Cudd_bddIthVar(_ddManager, targ);
	else var = Cudd_Not((Cudd_bddIthVar(_ddManager, targ)));

    for (int i = 0; i < _w; i++)
    {
        // Init C
        c = Cudd_ReadOne(_ddManager); 
        Cudd_Ref(c);
        tmp = Cudd_bddAnd(_ddManager, c, Cudd_Not(var));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, c);
        c = tmp;
        for (int j = 0; j < tensor->_r; j++)
        {
            // G
            g = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_Not(var));
            Cudd_Ref(g);
            // D
            term1 = Cudd_bddAnd(_ddManager, tensor->_allBDD[i][j], var);
            Cudd_Ref(term1);
            term2 = Cudd_Not(Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], var));
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(var));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            d = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(d);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

            // Detect overflow
            if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(g, d, c))
			{
				incBDDsBitWidth(tensor->_allBDD);
				++tensor->_r;
				isOverflow = 1;
			}
            // Sum
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, g, d);
            Cudd_Ref(tensor->_allBDD[i][j]);
            tmp = Cudd_bddXor(_ddManager, tensor->_allBDD[i][j], c);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            tensor->_allBDD[i][j] = tmp;
            // Carry
            if (j == tensor->_r - 1)
            {
                Cudd_RecursiveDeref(_ddManager, c);
                Cudd_RecursiveDeref(_ddManager, g);
                Cudd_RecursiveDeref(_ddManager, d);
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, g, d);
                Cudd_Ref(term1);
                term2 = Cudd_bddOr(_ddManager, g, d);
                Cudd_Ref(term2);
                Cudd_RecursiveDeref(_ddManager, g);
                Cudd_RecursiveDeref(_ddManager, d);
                tmp = Cudd_bddAnd(_ddManager, term2, c);
                Cudd_Ref(tmp);
                Cudd_RecursiveDeref(_ddManager, term2);
                Cudd_RecursiveDeref(_ddManager, c);
                term2 = tmp;
                c = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(c);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::Phase_shift(Tensor *tensor, const std::vector<int> &qubits, const int phase)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back();

	assert(targ >= 0 && targ < tensor->_rank);

    int nshift = _w / phase;
    int isOverflow = 0;

    DdNode *g, *c, *tmp, *term1, *term2;

	auto copy = tensor->_allBDD;
    for (int i = 0; i < _w; i++)
         for (int j = 0; j < tensor->_r; j++)
            Cudd_Ref(copy[i][j]);

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if (i >= _w - nshift)
        {
            c = Cudd_bddIthVar(_ddManager, targ);
            Cudd_Ref(c);
        }

        for (int j = 0; j < tensor->_r; j++)
        {
            if (i >= _w - nshift)
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, Cudd_Not(copy[i - (_w - nshift)][j]), Cudd_bddIthVar(_ddManager, targ));
                Cudd_Ref(term2);
                g = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(g);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);

                // Detect overflow
                if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(g, Cudd_Not(Cudd_ReadOne(_ddManager)), c))
				{
					incBDDsBitWidth(tensor->_allBDD);
					incBDDsBitWidth(copy);
					++tensor->_r;
					isOverflow = 1;
				}

                // Plus
                Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
                if (c == Cudd_Not(Cudd_ReadOne(_ddManager)))     // must be constant 0
                    tensor->_allBDD[i][j] = g;
                else
                {
                    // Sum
                    tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, g, c);
                    Cudd_Ref(tensor->_allBDD[i][j]);
                    // Carry
                    if (j == tensor->_r - 1)
                    {
                        Cudd_RecursiveDeref(_ddManager, g);
                        Cudd_RecursiveDeref(_ddManager, c);
                    }
                    else
                    {
                        tmp = Cudd_bddAnd(_ddManager, g, c);
                        Cudd_Ref(tmp);
                        Cudd_RecursiveDeref(_ddManager, g);
                        Cudd_RecursiveDeref(_ddManager, c);
                        c = tmp;
                    }
                }
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, copy[i + nshift][j], Cudd_bddIthVar(_ddManager, targ));
                Cudd_Ref(term2);

                Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
                tensor->_allBDD[i][j] = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(tensor->_allBDD[i][j]);

                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }
	
	for (int i = 0; i < _w; i++)
        for (int j = 0; j < tensor->_r; j++)
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::Phase_shift_dagger(Tensor* tensor, const std::vector<int> &qubits, const int phase)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back();

	assert(targ >= 0 && targ < tensor->_rank);

    int nshift = _w / abs(phase);
    int isOverflow = 0;

    DdNode *g, *c, *tmp, *term1, *term2;

	auto copy = tensor->_allBDD;
    for (int i = 0; i < _w; i++)
         for (int j = 0; j < tensor->_r; j++)
            Cudd_Ref(copy[i][j]);

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if (i < nshift)
        {
            c = Cudd_bddIthVar(_ddManager, targ);
            Cudd_Ref(c);
        }

        for (int j = 0; j < tensor->_r; j++)
        {
            if (i < nshift)
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, Cudd_Not(copy[_w - nshift + i][j]), Cudd_bddIthVar(_ddManager, targ));
                Cudd_Ref(term2);
                g = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(g);
                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);

                // Detect overflow
                if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(g, Cudd_Not(Cudd_ReadOne(_ddManager)), c))
				{
					incBDDsBitWidth(tensor->_allBDD);
					incBDDsBitWidth(copy);
					++tensor->_r;
					isOverflow = 1;
				}

                // Plus
                Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
                if (c == Cudd_Not(Cudd_ReadOne(_ddManager)))     // must be constant 0
                    tensor->_allBDD[i][j] = g;
                else
                {
                    // Sum
                    tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, g, c);
                    Cudd_Ref(tensor->_allBDD[i][j]);
                    // Carry
                    if (j == tensor->_r - 1)
                    {
                        Cudd_RecursiveDeref(_ddManager, g);
                        Cudd_RecursiveDeref(_ddManager, c);
                    }
                    else
                    {
                        tmp = Cudd_bddAnd(_ddManager, g, c);
                        Cudd_Ref(tmp);
                        Cudd_RecursiveDeref(_ddManager, g);
                        Cudd_RecursiveDeref(_ddManager, c);
                        c = tmp;
                    }
                }
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, copy[i - nshift][j], Cudd_bddIthVar(_ddManager, targ));
                Cudd_Ref(term2);

                Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
                tensor->_allBDD[i][j] = Cudd_bddOr(_ddManager, term1, term2);
                Cudd_Ref(tensor->_allBDD[i][j]);

                Cudd_RecursiveDeref(_ddManager, term1);
                Cudd_RecursiveDeref(_ddManager, term2);
            }
        }
    }
	
	for (int i = 0; i < _w; i++)
        for (int j = 0; j < tensor->_r; j++)
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::PauliX(Tensor *tensor, const std::vector<int> &qubits)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back();

	assert(targ >= 0 && targ < tensor->_rank);

    DdNode *tmp, *term1, *term2;

    for (int i = 0; i < _w; i++) // F = tensor->_allBDD[i][j]
    {
        for (int j = 0; j < tensor->_r; j++)
        {
            // term1
            term1 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(term1);

            tmp = Cudd_bddAnd(_ddManager, term1, Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;

            // term2
            term2 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_bddIthVar(_ddManager, targ));
            Cudd_Ref(term2);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(Cudd_bddIthVar(_ddManager, targ)));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            // OR
            tensor->_allBDD[i][j] = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(tensor->_allBDD[i][j]);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
        }
    }
}

void BDDSystem::PauliY(Tensor *tensor, const std::vector<int> &qubits, const bool fTranspose)
{
	assert(qubits.size() == 1);
	const int targ = qubits.back();

	assert(targ >= 0 && targ < tensor->_rank);

    int nshift = _w / 2;

    DdNode *g, *c, *tmp, *term1, *term2, *var;
    int isOverflow = 0;

    if (fTranspose) var = Cudd_Not(Cudd_bddIthVar(_ddManager, targ));
    else var = Cudd_bddIthVar(_ddManager, targ);

    // PauliX(targ);
    for (int i = 0; i < _w; i++) // F = allBDD[i][j]
    {
        for (int j = 0; j < tensor->_r; j++)
        {
            // term1
            term1 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], Cudd_Not(var));
            Cudd_Ref(term1);

            tmp = Cudd_bddAnd(_ddManager, term1, var);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term1);
            term1 = tmp;

            // term2
            term2 = Cudd_Cofactor(_ddManager, tensor->_allBDD[i][j], var);
            Cudd_Ref(term2);
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);

            tmp = Cudd_bddAnd(_ddManager, term2, Cudd_Not(var));
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;

            // OR
            tensor->_allBDD[i][j] = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(tensor->_allBDD[i][j]);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);
        }
    }

	auto copy = tensor->_allBDD;
    for (int i = 0; i < _w; i++)
         for (int j = 0; j < tensor->_r; j++)
            Cudd_Ref(copy[i][j]);

    for (int i = 0; i < _w; i++)
    {
        // Init C
        if (i < nshift)
            c = Cudd_Not(var);
        else
            c = Cudd_bddAnd(_ddManager, Cudd_ReadOne(_ddManager), var);
        Cudd_Ref(c);

        for (int j = 0; j < tensor->_r; j++)
        {
            if (i < nshift)
            {
                term1 = Cudd_bddAnd(_ddManager, copy[i + nshift][j], var);
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, Cudd_Not(copy[i + nshift][j]), Cudd_Not(var));
                Cudd_Ref(term2);
            }
            else
            {
                term1 = Cudd_bddAnd(_ddManager, Cudd_Not(copy[i - nshift][j]), var);
                Cudd_Ref(term1);
                term2 = Cudd_bddAnd(_ddManager, copy[i - nshift][j], Cudd_Not(var));
                Cudd_Ref(term2);
            }
            g = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(g);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

            // Detect overflow
            if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(g, Cudd_Not(Cudd_ReadOne(_ddManager)), c))
			{
				incBDDsBitWidth(tensor->_allBDD); 
				incBDDsBitWidth(copy);
				++tensor->_r;
				isOverflow = 1;
			}

            // Plus 1
            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            if (c == Cudd_Not(Cudd_ReadOne(_ddManager)))
                tensor->_allBDD[i][j] = g;
            else
            {
                // Sum
                tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, g, c);
                Cudd_Ref(tensor->_allBDD[i][j]);
                // Carry
                if (j == tensor->_r - 1)
                {
                    Cudd_RecursiveDeref(_ddManager, g);
                    Cudd_RecursiveDeref(_ddManager, c);
                }
                else
                {
                    tmp = Cudd_bddAnd(_ddManager, g, c);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(_ddManager, g);
                    Cudd_RecursiveDeref(_ddManager, c);
                    c = tmp;
                }
            }
        }
    }
	
	for (int i = 0; i < _w; i++)
        for (int j = 0; j < tensor->_r; j++)
            Cudd_RecursiveDeref(_ddManager, copy[i][j]);

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::PauliZ(Tensor *tensor, const std::vector<int> &qubits)
{
	for(const auto ele: qubits)
		assert(ele >= 0 && ele < tensor->_rank);

    DdNode *c, *tmp, *term1, *term2, *inter, *qubit_and;

    // Init qubit and
    qubit_and = Cudd_ReadOne(_ddManager); 
    Cudd_Ref(qubit_and);
    for (int i = qubits.size() - 1; i >= 0; i--)
    {
        tmp = Cudd_bddAnd(_ddManager, qubit_and, Cudd_bddIthVar(_ddManager, qubits[i]));
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, qubit_and);
        qubit_and = tmp;
    }
    int isOverflow = 0;
    for (int i = 0; i < _w; i++)
    {
        // Init C
        c = Cudd_ReadOne(_ddManager); 
        Cudd_Ref(c);
        tmp = Cudd_bddAnd(_ddManager, c, qubit_and);
        Cudd_Ref(tmp);
        Cudd_RecursiveDeref(_ddManager, c);
        c = tmp;
        for (int j = 0; j < tensor->_r; j++)
        {
            term1 = Cudd_bddAnd(_ddManager, tensor->_allBDD[i][j], Cudd_Not(qubit_and));
            Cudd_Ref(term1);
            term2 = Cudd_Not(tensor->_allBDD[i][j]);
            Cudd_Ref(term2);
            tmp = Cudd_bddAnd(_ddManager, term2, qubit_and);
            Cudd_Ref(tmp);
            Cudd_RecursiveDeref(_ddManager, term2);
            term2 = tmp;
            inter = Cudd_bddOr(_ddManager, term1, term2);
            Cudd_Ref(inter);
            Cudd_RecursiveDeref(_ddManager, term1);
            Cudd_RecursiveDeref(_ddManager, term2);

			// Detect overflow
			if ((j == tensor->_r - 1) && !isOverflow && checkAdderOverflow(inter, Cudd_Not(Cudd_ReadOne(_ddManager)), c))
			{
				incBDDsBitWidth(tensor->_allBDD);
				++tensor->_r;
				isOverflow = 1;
			}

            Cudd_RecursiveDeref(_ddManager, tensor->_allBDD[i][j]);
            // Plus 1
            if (c == Cudd_Not(Cudd_ReadOne(_ddManager)))
                tensor->_allBDD[i][j] = inter;
            else
            {
                // Sum
                tensor->_allBDD[i][j] = Cudd_bddXor(_ddManager, inter, c);
                Cudd_Ref(tensor->_allBDD[i][j]);
                // Carry
                if (j == tensor->_r - 1)
				{
                    Cudd_RecursiveDeref(_ddManager, inter);
					Cudd_RecursiveDeref(_ddManager, c);
				}
                else
                {
                    tmp = Cudd_bddAnd(_ddManager, inter, c);
                    Cudd_Ref(tmp);
                    Cudd_RecursiveDeref(_ddManager, c);
                    Cudd_RecursiveDeref(_ddManager, inter);
                    c = tmp;
                }
            }
        }
    }
    Cudd_RecursiveDeref(_ddManager, qubit_and);

	if(isOverflow)
		dropTensorBits(tensor);
}

void BDDSystem::applyGate(const Gate* gate, Tensor *tensor, bool fTranspose)
{
	const GateType gateType = gate->getType();
	const std::vector<int> &qubits = gate->getQubits();

	if (gateType == GateType::X) PauliX(tensor, qubits);
	else if (gateType == GateType::Y) PauliY(tensor, qubits, fTranspose);
	else if (gateType == GateType::Z) PauliZ(tensor, qubits);
	else if (gateType == GateType::H) Hadamard(tensor, qubits);
	else if (gateType == GateType::S) Phase_shift(tensor, qubits, 2);
	else if (gateType == GateType::SDG) Phase_shift_dagger(tensor, qubits, -2);
	else if (gateType == GateType::T) Phase_shift(tensor, qubits, 4);
	else if (gateType == GateType::TDG) Phase_shift_dagger(tensor, qubits, -4);
	else if (gateType == GateType::RX_PI_2) Rx_pi_2(tensor, qubits);
	else if (gateType == GateType::RX_PI_2_DG) Rx_pi_2_dagger(tensor, qubits);
	else if (gateType == GateType::RY_PI_2) Ry_pi_2(tensor, qubits, fTranspose);
	else if (gateType == GateType::RY_PI_2_DG) Ry_pi_2_dagger(tensor, qubits, fTranspose);
	else if (gateType == GateType::CX) Toffoli(tensor, qubits);
	else if (gateType == GateType::CZ) PauliZ(tensor, qubits);
	else if (gateType == GateType::SWAP) Fredkin(tensor, qubits);
	else if (gateType == GateType::CSWAP) Fredkin(tensor, qubits);
	else if (gateType == GateType::CCX) Toffoli(tensor, qubits);

	updateMaxNodeCount();
}
