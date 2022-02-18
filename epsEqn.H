/*******************************************************
Calculates the evolution of aggregation intensity over time.

Wall shear stress tau_w
*******************************************************/
Info<<"Calculating wall shear stress\n"<<endl;

tau_w = nu*rho*mag(fvc::grad(U));     //Wall shear stress in interior cells

/*******************************************************



*******************************************************/

//Thrombus susceptability equation
const fvPatchList& patches = mesh.boundary();

forAll(patches, patchI)
{
    const fvPatch& curPatch = patches[patchI];
    if(isType<wallFvPatch>(curPatch))
    {
	forAll(curPatch,faceI)
	{
	    label faceCellI = curPatch.faceCells()[faceI];
	    check[faceCellI] = 1.0;
	    if(tau_w[faceCellI] >= tau_high_wall.value())
	    {
	        P_TSP[faceCellI] = 0.0;
	    }
	    else if (tau_w[faceCellI] <= tau_low_wall.value())
	    {
	        P_TSP[faceCellI] = 1.0;
	    }
	    else
	    {
	        P_TSP[faceCellI] = (tau_w[faceCellI]/tau_high_thr.value())*(Foam::exp((tau_w[faceCellI]-tau_low_thr.value())/(tau_high_thr.value() - tau_low_thr.value()))-1.0);
	        P_TSP[faceCellI] /= (Foam::constant::mathematical::e - 1.0);
	        P_TSP[faceCellI] = 1.0 - P_TSP[faceCellI];
	    }

	    //Thrombus breakdown rate
	    if(tau_w[faceCellI] >= tau_bd_thr.value())
	    {
		beta_epsilon[faceCellI] = B.value();
	    }
	    else
	    {
		beta_epsilon[faceCellI] = 0.0;
	    }
	}
    }
}

/******************************************************************************
The following block of code calculates the probability of thrombosis within/on
the surface of the thrombus. 
******************************************************************************/

forAll(eps, cellI)
{ 
    if (check[cellI] == 1.0)
        continue;
    
    if (eps[cellI] >= epsilon_t.value())
    {
        if (tau_w[cellI] >= tau_high_thr.value())
        {
            P_TSP[cellI] = 0.0;
        }
        else if (tau_w[cellI] <= tau_low_thr.value())
        {
            P_TSP[cellI] = 1.0;
        }
        else
        {
            P_TSP[cellI] = (tau_w[cellI]/tau_high_thr.value())*(Foam::exp((tau_w[cellI] - tau_low_thr.value())/(tau_high_thr.value() - tau_low_thr.value())) - 1.0);
            P_TSP[cellI] /= (Foam::constant::mathematical::e - 1.0);
            P_TSP[cellI] = 1.0 - P_TSP[cellI];
        }
    }
    else    
    {
        const labelList& adjCells = mesh.cellCells()[cellI];
        forAll(adjCells, i)
        {
            if (eps[adjCells[i]] >= epsilon_t.value())
            {
                if (tau_w[cellI] >= tau_high_thr.value())
                {
                    P_TSP[cellI] = 0.0;
                }
                else if (tau_w[cellI] <= tau_low_thr.value())
                {
                    P_TSP[cellI] = 1.0;
                }
                else
                {
                    P_TSP[cellI] = (tau_w[cellI]/tau_high_thr.value())*(Foam::exp((tau_w[cellI] - tau_low_thr.value())/(tau_high_thr.value() - tau_low_thr.value())) - 1.0);
                    P_TSP[cellI] /= (Foam::constant::mathematical::e - 1.0);
                    P_TSP[cellI] = 1.0 - P_TSP[cellI];
                }
                break;
            }
        }
    }
    
    //Breakdown    
    if (tau_w[cellI] >= tau_bd_thr.value())
    {
        beta_epsilon[cellI] = B.value();
    }
    else
    {
        beta_epsilon[cellI] = 0.0;
    }
}

/******************************************************************************
Solve the epsilon equation
******************************************************************************/

Info<< "Solving eps equation\n" << endl;
solve(fvm::ddt(eps) == alpha_epsilon*P_TSP*(AP*AP) - fvm::Sp(beta_epsilon, eps));
Info<< "Max/min eps: " << max(eps).value() << "/" << min(eps).value() << nl << endl;






