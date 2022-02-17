/*************************************************
List of all volume fields and scalars.
Note that parameters with discrete files in the "0" directories
don't specify a value or units, as they are specified elsewhere.

Naming scheme:
1. subscript >> A_B
2. "," >> ignore
3. Greek letters >> full name of Greek letter
4. Uppercase Greek letters >> capital name of Greek letter
**************************************************/


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field Resting Platelets\n" << endl;
volScalarField RP
(
    IOobject
    (
        "RP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Activated Platelets\n" << endl;
volScalarField AP
(
    IOobject
    (
        "AP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Adenosine Diphosphate\n" << endl;
volScalarField ADP
(
    IOobject
    (
        "ADP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Thromboxane\n" << endl;
volScalarField TxA_2
(
    IOobject
    (
        "TxA_2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Prothrombin\n" << endl;
volScalarField PT
(
    IOobject
    (
        "PT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Thrombin\n" << endl;
volScalarField TRB
(
    IOobject
    (
        "TRB",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Antithrombin III\n" << endl;
volScalarField AT
(
    IOobject
    (
        "AT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Omega\n" << endl;
volScalarField Omega
(
    IOobject
    (
        "Omega",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Omega", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

Info<< "Reading field epsilon\n" << endl;
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field k_pa\n" << endl;
volScalarField k_pa
(
    IOobject
    (
        "k_pa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("k_pa", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field Gamma\n" << endl;
volScalarField Gamma
(
    IOobject
    (
        "Gamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Gamma", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

/*
Info<< "Reading field tau\n" << endl;
volScalarField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tau", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);
*/


Info<< "Reading field P_TSP\n" << endl;
volScalarField P_TSP
(
    IOobject
    (
        "P_TSP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("P_TSP", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

Info<< "Reading field beta_epsilon\n" << endl;
volScalarField beta_epsilon
(
    IOobject
    (
        "beta_epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("beta_epsilon", dimensionSet(0,0,-1,0,0,0,0), 0.0)
);

Info<< "Reading field wall shear stress\n" << endl;
volScalarField tau_w
(
    IOobject
    (
        "tau_w",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tau_w", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

Info<< "Reading field check\n" << endl;
volScalarField check
(
    IOobject
    (
        "check",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("check", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


dimensionedScalar nu
(
    "nu",
    transportProperties.lookup("nu")
);

dimensionedScalar D_AP
(
    "D_AP",
    transportProperties.lookup("D_AP")
);

dimensionedScalar D_RP
(
    "D_RP",
    transportProperties.lookup("D_RP")
);

dimensionedScalar D_ADP
(
    "D_ADP",
    transportProperties.lookup("D_ADP")
);

dimensionedScalar D_TxA_2
(
    "D_TxA_2",
    transportProperties.lookup("D_TxA_2")
);

dimensionedScalar D_PT
(
    "D_PT",
    transportProperties.lookup("D_PT")
);

dimensionedScalar D_TRB
(
    "D_TRB",
    transportProperties.lookup("D_TRB")
);

dimensionedScalar D_AT
(
    "D_AT",
    transportProperties.lookup("D_AT")
);

dimensionedScalar lambda_j
(
    "lambda_j",
    transportProperties.lookup("lambda_j")
);

dimensionedScalar k_1j
(
    "k_1j",
    transportProperties.lookup("k_1j")
);

dimensionedScalar s_pj
(
    "s_pj",
    transportProperties.lookup("s_pj")
);

dimensionedScalar beta
(
    "beta",
    transportProperties.lookup("beta")
);

dimensionedScalar phi_at
(
    "phi_at",
    transportProperties.lookup("phi_at")
);

dimensionedScalar phi_rt
(
    "phi_rt",
    transportProperties.lookup("phi_rt")
);

dimensionedScalar t_act
(
    "t_act",
    transportProperties.lookup("t_act")
);

dimensionedScalar k_1T
(
    "k_1T",
    transportProperties.lookup("k_1T")
);

dimensionedScalar H
(
    "H",
    transportProperties.lookup("H")
);

dimensionedScalar alpha
(
    "alpha",
    transportProperties.lookup("alpha")
);

dimensionedScalar K_AT
(
    "K_AT",
    transportProperties.lookup("K_AT")
);

dimensionedScalar K_T
(
    "K_T",
    transportProperties.lookup("K_T")
);

dimensionedScalar k
(
    "k",
    transportProperties.lookup("k")
);

dimensionedScalar alpha_epsilon
(
    "alpha_epsilon",
    transportProperties.lookup("alpha_epsilon")
);

dimensionedScalar B
(
    "B",
    transportProperties.lookup("B")
);

dimensionedScalar tau_low_wall
(
    "tau_low_wall",
    transportProperties.lookup("tau_low_wall")
);

dimensionedScalar tau_high_wall
(
    "tau_high_wall",
    transportProperties.lookup("tau_high_wall")
);

dimensionedScalar tau_bd_wall
(
    "tau_bd_wall",
    transportProperties.lookup("tau_bd_wall")
);

dimensionedScalar tau_low_thr
(
    "tau_low_thr",
    transportProperties.lookup("tau_low_thr")
);

dimensionedScalar tau_high_thr
(
    "tau_high_thr",
    transportProperties.lookup("tau_high_thr")
);

dimensionedScalar tau_bd_thr
(
    "tau_bd_thr",
    transportProperties.lookup("tau_bd_thr")
);

dimensionedScalar epsilon_t
(
    "epsilon_t",
    transportProperties.lookup("epsilon_t")
);

dimensionedScalar rho
(
    "rho",
    transportProperties.lookup("rho")
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::momentumTransportModel> turbulence
(
    incompressible::momentumTransportModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
