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

#include "TaylorGreenVortex2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TaylorGreenVortex2D, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void::Foam::TaylorGreenVortex2D::setPropertiesOutput()
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
        << "time" << tab << "Ek" << tab << "epsilon" << endl;

    globalPropertiesFile_().precision(8); // set precision

    Info << "TGV2D: Global properties are written in\n"
        << globalPropertiesFile_().name() << endl;
}

void Foam::TaylorGreenVortex2D::calcGlobalProperties()
{
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
                fvc::laplacian(nuEff, U_)
              + fvc::div(nuEff*dev(T(fvc::grad(U_))))
            )
        );

    // dissipation rate weighted average
    epsilon_ = dissipation.weightedAverage(mesh_.V())/(pow(Uinit_,3)/L_);

    // Kinetic energy weighted average
    Ek_ = 0.5*magSqr(U_)().weightedAverage(mesh_.V())/ pow(Uinit_,2);


    // write to log file
    if(Pstream::master())
    {
        globalPropertiesFile_()
            << runTime_.value()/tc_.value() << tab
            << Ek_.value() << tab << epsilon_.value()  << tab

            << endl;
    }
}


void::Foam::TaylorGreenVortex2D::setInitialFields
(
    volVectorField UIn,
    surfaceScalarField phiIn,
    volScalarField pIn
)
{
    Ua_ = UIn;
    phia_ = phiIn;
    pa_ = pIn;
}


void::Foam::TaylorGreenVortex2D::setInitialFieldsAsAnalytical()
{
    tmp<volScalarField> x_c = mesh_.C().component(vector::X);
    tmp<volScalarField> y_c = mesh_.C().component(vector::Y);
    tmp<volScalarField> z_c = mesh_.C().component(vector::Z);

    Ua_ =  Uinit_*( vector(1,0,0) * sin(x_c/L_) * cos(y_c/L_)
                  - vector(0,1,0) * cos(x_c/L_) * sin(y_c/L_)
                  + vector(0,0,1) * scalar(0.));

    phia_ = fvc::interpolate(Ua_)& mesh_.Sf();

    pa_ = sqr(Uinit_)/4 * (cos(2.*x_c/L_) + cos(2.*y_c/L_));
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TaylorGreenVortex2D::TaylorGreenVortex2D
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

    Ua_
    (
        IOobject
        (
            "UStart",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_
    ),

    phia_
    (
        IOobject
        (
            "phiStart",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi_
    ),

    pa_
    (
        IOobject
        (
            "pStart",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p_
    ),

    TGV2DProperties_
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

    Uinit_( "Uinit", dimVelocity,TGV2DProperties_),

    L_("L", dimLength, TGV2DProperties_),

    tc_("tc", dimTime, L_.value()/Uinit_.value()),

    epsilon_ ("epsilon", dimless, 0.0),

    Ek_ ("Ek", dimless, 0.0)



{
    Info << "TaylorGreenVortex2D constructor" << endl;


    setPropertiesOutput();

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TaylorGreenVortex2D::~TaylorGreenVortex2D()
{
    Info << "TaylorGreenVortex2D Destructor" << endl;
}


// ************************************************************************* //
