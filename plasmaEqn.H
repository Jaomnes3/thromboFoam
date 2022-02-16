/*************************************************************
The activation of platelets due to the activity of chemicals is quantified. 
The exact order was decided in reference to Yingming(2020):
Omega >> Kpa >> Gamma
**************************************************************/

Info<< "\nCaulculating platelet activation rate kpa\n" << endl;

forAll(Omega,[i])
{
    Omega[i] = ADP[i]/2.0 + TxA2[i]/0.6 + TRB[i]/0.1;
}

forAll(k_pa,[i])
{
    if(Omega[i] < 1.0)
	k_pa[i] = 0.0;
    else
	k_pa[i] = Omega[i]/t_act.value();
}

forAll(Gamma,[i])
{
    Gamma[i] = k_1j.value()*H.value()*AT[i]/
               ( alpha.value()*K_AT.value()*K_T.value()
                +alpha.value()*K_AT.value()*beta.value()*TRB[i]
		+AT[i]*beta.value()*TRB[i]);
}

/************************************************************
Calculation of shear stress by multiplying dynamic viscosity and 
the norm of the velocity gradient.

Copied from clotFoam by Saurev(2016).

Interesting observations regarding OpenFoam coding:
1. The += operator is used for positive updates.
	a += b >> a = a + b
2. The xy separation of the tauTensor using the ~.xx() notation.
************************************************************/

Info<<"n\Calculating mechanical activation\n" << endl;

volTensorField tauTensor = nu*rho*fvc::grad(U);
tauTensor += tauTensor.T();

forAll(tau,cellI)
{
    tau[cellI] = tauTensor[cellI].xx()*tauTensor[cellI].xx()
                +tauTensor[cellI].yy()*tauTensor[cellI].yy()
		+tauTensor[cellI].zz()*tauTensor[cellI].zz()
		-tauTensor[cellI].xx()*tauTensor[cellI].yy()
		-tauTensor[cellI].xx()*tauTensor[cellI].zz()
		-tauTensor[cellI].yy()*tauTensor[cellI].zz()
	  +3.0*( tauTensor[cellI].xy()*tauTensor[cellI].xy()
		+tauTensor[cellI].xz()*tauTensor[cellI].xz()
		+tauTensor[cellI].yz()*tauTensor[cellI].yz());
    
    tau[cellI] = Foam::sqrt(tau[cellI]/3.0);
}

/************************************************************
Chemical reactions for 7 species equations.
************************************************************/

//Resting platelets
Info<< "\nSolving resting platelet equation" << endl;
solve(fvm::ddt(RP) + fvm::div(phi,RP) == fvm::laplacian(D_RP,RP)-k_pa*RP);
Info<< "Max/min resting platelets: " << max(RP).value() << "/" << min(RP).value() << nl << endl;

//Activated platelets
Info<< "\nSolving activated platelet equation" << endl;
solve(fvm::ddt(AP) + fvm::div(phi,AP) == fvm::laplacian(D_AP,AP)+k_pa*RP);
Info<< "Max/min resting platelets: " << max(AP).value() << "/" << min(AP).value() << nl << endl;

//Adenosine Diphosphate
Info<< "\nSolving ADP equation" << endl;
solve(fvm::ddt(ADP) + fvm::div(phi,ADP) == fvm::laplacian(D_ADP,ADP)-lambda_j*k_pa*RP-k_1j*ADP);
Info<< "Max/min ADP: " << max(ADP).value() << "/" << min(ADP).value() << nl << endl;

//Thromboxane A_2
Info<< "\nSolving TxA_2 equation" << endl;
solve(fvm::ddt(TxA_2) + fvm::div(phi,TxA_2) == fvm::laplacian(D_TxA_2,TxA_2)+ s_pj*AP - k_1j*TxA_2);
Info<< "Max/min thromboxane: " << max(TxA_2).value() << "/" << min(TxA_2).value() << nl << endl;

//Prothrombin
Info<< "\nSolving PT equation" << endl;
solve(fvm::ddt(PT) + fvm::di(phi,PT) == fvm::laplacian(D_PT,PT) - beta*PT* (phi_at*AP + phi_rt*RP));
Info<< "Max/min prothrombin: " << max(PT).value() << "/" << min(PT).value() << nl << endl;

//Thrombin
Info<< "\nSolving TRB equation" << endl;
solve(fvm::ddt(TRB) + fvm::div(phi,TRB) == fvm::laplacian(D_TRB,TRB) -Gamma*TRB + PT*(phi_at*AP + phi_rt*RP));
Info<< "Max/min thrombin: " << max(TRB).value() << "/" << min(TRB).value() << nl << endl;

//Antithrombin III
Info<< "\nSolving AT equation" << endl;
solve(fvm::ddt(AT) + fvm::div(phi,AT) == fvm::laplacian(D_AT,AT) -Gamma*beta*TRB);
Info<< "Max/min antithrombin: " << max(AT).value() << "/" << min(AT).value() << nl << endl;


