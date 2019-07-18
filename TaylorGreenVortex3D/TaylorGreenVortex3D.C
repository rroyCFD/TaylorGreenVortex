/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "TaylorGreenVortex3D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TaylorGreenVortex3D, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::TaylorGreenVortex3D::setInitialFieldsAsAnalytical
(
    volVectorField UIn, surfaceScalarField phiIn, volScalarField pIn
)
{
    tmp<volScalarField> x_c = mesh_.C().component(vector::X);
    tmp<volScalarField> y_c = mesh_.C().component(vector::Y);
    tmp<volScalarField> z_c = mesh_.C().component(vector::Z);

    UIn =  Uinit_*(
        vector(1,0,0) * sin(x_c.ref()/L_) * cos(y_c.ref()/L_)* cos(z_c.ref()/L_)
      - vector(0,1,0) * cos(x_c.ref()/L_) * sin(y_c.ref()/L_)* cos(z_c.ref()/L_)
      + vector(0,0,1) * scalar(0.));

    phiIn = fvc::interpolate(U_)& mesh_.Sf();

    pIn = sqr(Uinit_)/16 *
        (cos(2.*x_c.ref()/L_) + cos(2.*y_c.ref()/L_))*(cos(2*z_c.ref()/L_) + 2);

    x_c.clear();
    y_c.clear();
    z_c.clear();
}


void Foam::TaylorGreenVortex3D::setPropertiesOutput()
{
    // create output file
    fileName outputDir;
    autoPtr<fileName> outFilePath_;

    if(Pstream::parRun() && Pstream::master())
    {
        outputDir = runTime_.path()/"../postProcessing";

        outFilePath_.reset( new fileName
        (
            runTime_.path()/"../postProcessing"/
            ("TaylorGreenVortexProperties.dat_"+runTime_.timeName())
        )
        );

    }
    else
    {
        outputDir = runTime_.path()/"postProcessing";

        outFilePath_.reset( new fileName
        (
            runTime_.path()/"postProcessing"/
            ("TaylorGreenVortexProperties.dat_"+runTime_.timeName())
        )
        );
    }

    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    globalPropertiesFile_.reset(new OFstream(*outFilePath_));

    // write header of the log file
    globalPropertiesFile_()
        << "time" << tab << "Ek(t)" << tab << "epsilon(t)" << tab << endl;

    Info << "TGV3D: Global properties are written in\n"
         << globalPropertiesFile_().name() << endl;
}


void Foam::TaylorGreenVortex3D::calcGlobalProperties()
{

    Info << "Calculating kinetic energy and dissipation rate" << endl;

    // Assuming laminar 3D Taylor-Green case
    // volScalarField dissipation =
    //     (U_ &
    //         (
    //             fvc::laplacian(nu_, U_)
    //           + fvc::div(nu_*dev(T(fvc::grad(U_))))
    //         )
    //     );


    // Generic turbulence model case
    tmp<volScalarField> nuEff
    (
        mesh_.lookupObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        ).nuEff()
    );

    volScalarField dissipation =
        (U_ &
            (
                fvc::laplacian(nuEff.ref(), U_)
              + fvc::div(nuEff.ref()*dev(T(fvc::grad(U_))))
            )
        );

    nuEff.clear();

    // dissipation rate weighted average
    epsilon_ = dissipation.weightedAverage(mesh_.V())/(pow(Uinit_,3)/L_);

    // Kinetic energy weighted average
    Ek_ = 0.5*magSqr(U_)().weightedAverage(mesh_.V())/ pow(Uinit_,2);
}

void Foam::TaylorGreenVortex3D::write()
{
    Info << "writing to log file" << endl;
    if(Pstream::master())
    {
        globalPropertiesFile_()
            << runTime_.value()/tc_.value() << tab
            << Ek_.value() << tab << epsilon_.value() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TaylorGreenVortex3D::TaylorGreenVortex3D
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& p
)
:
    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    // Set the pointer to the velocity, flux and pressure-correction field
    U_(U),
    phi_(phi),
    p_(p),

    TGV3DProperties_
    (
        IOobject
        (
            "TaylorGreenVortexProperties",
            runTime_.constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    Uinit_( "Uinit", dimVelocity,TGV3DProperties_),

    L_("L", dimLength, TGV3DProperties_),

    tc_("tc", dimTime, L_.value()/Uinit_.value()),

    epsilon_ ("epsilon", dimless, 0.0),

    Ek_ ("Ek", dimless, 0.0)

{
    Info << "TaylorGreenVortex3D constructor" << endl;

    // Set global properties output file
    setPropertiesOutput();

    // Get molucular viscosity
    const dictionary& transportProperties_ =
        mesh_.lookupObject<dictionary>("transportProperties");

    nu_.reset(new dimensionedScalar ("nu", dimViscosity, transportProperties_));

    // calculate Global Ek and epsilon values
    calcGlobalProperties();

    // write to log files
    write();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TaylorGreenVortex3D::~TaylorGreenVortex3D()
{
    Info << "TaylorGreenVortex3D Destructor" << endl;
}

// ************************************************************************* //
