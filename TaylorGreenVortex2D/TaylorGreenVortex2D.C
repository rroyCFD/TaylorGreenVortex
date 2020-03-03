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
void Foam::TaylorGreenVortex2D::setAnalyticalFields()
{
    tmp<volScalarField> x_c = mesh_.C().component(vector::X);
    tmp<volScalarField> y_c = mesh_.C().component(vector::Y);

    Ua_ =  Uinit_*( vector(1,0,0) * sin(x_c.ref()/L_) * cos(y_c.ref()/L_)
                  - vector(0,1,0) * cos(x_c.ref()/L_) * sin(y_c.ref()/L_)
                  + vector(0,0,1) * scalar(0.));

    phia_ = fvc::interpolate(Ua_)& mesh_.Sf();

    pa_ = sqr(Uinit_)/4 * (cos(2.*x_c.ref()/L_) + cos(2.*y_c.ref()/L_));

    x_c.clear();
    y_c.clear();

    calcAnalytical_ =false;

    return;
}


void Foam::TaylorGreenVortex2D::getNu()
{
    // Get molucular viscosity
    const dictionary& transportProperties_ =
        mesh_.lookupObject<dictionary>("transportProperties");

    nu_.reset(new dimensionedScalar ("nu", dimViscosity, transportProperties_));
}


void Foam::TaylorGreenVortex2D::calcErrorNorms()
{
    // Velocity error norms
    Info << "Calculating Uerror" << endl;
    Uerror_ = Ua_*Foam::exp(-2.0*nu_().value()*runTime_.value()) - U_;

    ULinfErr_.value() = gMax(mag(Uerror_)());   // infinity norm
    UL2err_ = sqrt(magSqr(Uerror_)().weightedAverage(mesh_.V())); // 2nd norm


    // Pressure error norms
    Info << "Calculating perror" << endl;
    if(pRefOn_)
    {
        tmp<volScalarField> pAna_ =
            pa_*Foam::exp(-4.0*nu_().value()*runTime_.value());

        dimensionedScalar pRefValAna(
            "", p_.dimensions(), getRefCellValue(pAna_.ref(), pRefCell_));
        dimensionedScalar pRefValSim(
            "", p_.dimensions(), getRefCellValue(p_, pRefCell_));

        perror_ = pAna_.ref() - p_ - pRefValAna + pRefValSim; //

        pAna_.clear();
    }
    else
    {
        perror_ = pa_*Foam::exp(-4.0*nu_().value()*runTime_.value()) - p_;
    }

    pLinfErr_.value() = gMax(mag(perror_)()); // infinity norm
    pL2err_ = sqrt(magSqr(perror_)().weightedAverage(mesh_.V()));  // 2nd norm

    return;
}


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //


void Foam::TaylorGreenVortex2D::setInitialFieldsAsAnalytical(
    volVectorField& UIn,
    surfaceScalarField& phiIn,
    volScalarField& pIn
)
{
    if(calcAnalytical_)
    {
        setAnalyticalFields();
    }

    Info << "Time " << runTime_.timeName()
         << ": making U and p analytical; phi is interpolated." << endl;

    if(runTime_.startTime().value() == 0)
    {
       pIn  = pa_;
       UIn   = Ua_;
       phiIn = phia_;
    }
    else
    {
        pIn   = pa_*Foam::exp(-4.0*nu_().value()*runTime_.value());
        UIn  = Ua_*Foam::exp(-2.0*nu_().value()*runTime_.value());
        phiIn = phia_*Foam::exp(-2.0*nu_().value()*runTime_.value());
    }

    UIn.write();
    phiIn.write();
    pIn.write();

    // Calculate and write global properties
        // calculate Global Ek and epsilon values
        calcGlobalProperties();

        // write to log files
        writeGlobalProperties();

    return;
}


void Foam::TaylorGreenVortex2D::setPropertiesOutput()
{
    if(Pstream::parRun() && !(Pstream::master()))
    {
        return;
    }
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
        << "time"      << tab
        << "Ek(t)"     << tab << "epsilon(t)" << tab
        << "L2(U,t)"   << tab << "L2(p,t)"    << tab
        << "Linf(U,t)" << tab << "Linf(p,t)"
        << endl;

    globalPropertiesFile_().precision(12); // set precision

    Info << "TGV2D: Global properties are written in\n"
         << globalPropertiesFile_().name() << endl;

    return;
}



void Foam::TaylorGreenVortex2D::calcGlobalProperties()
{

    Info << "Calculating kinetic energy and dissipation rate" << endl;

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


    // calculate error fields and error norms
    calcErrorNorms();

    return;
}



void Foam::TaylorGreenVortex2D::writeGlobalProperties()
{
    // Write error fields
    if(runTime_.outputTime())
    {
        perror_.write();
        Uerror_.write();
    }

    // return if not the master directory
    if(Pstream::parRun() && !(Pstream::master()))
    {
        return;
    }

    Info << "writing to log file" << endl;
    if(Pstream::master())
    {
        globalPropertiesFile_()
            << runTime_.value()/tc_.value() << tab
            << Ek_.value()       << tab << epsilon_.value() << tab
            << UL2err_.value()   << tab << pL2err_.value()  << tab
            << ULinfErr_.value() << tab << pLinfErr_.value()
            << endl;
    }

    return;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TaylorGreenVortex2D::TaylorGreenVortex2D
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& p,
    const label& pRefCell
)
:
    // construct base class
    TaylorGreenVortex (U),

    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    // Set the pointer to the velocity, flux and pressure-correction field
    U_(U),
    phi_(phi),
    p_(p),
    pRefCell_(pRefCell),
    pRefOn_(p.needReference()),

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

    calcAnalytical_(true),

    Uerror_
    (
        IOobject
        (
            "Uerror",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE // AUTO_WRITE
        ),
        (Ua_ - U_)
    ),

    perror_
    (
        IOobject
        (
            "perror",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE // AUTO_WRITE
        ),
        (pa_ - p_)
    ),

    ULinfErr_(max(mag(Uerror_))),
    UL2err_(sqrt(sum(magSqr(Uerror_)))),

    pLinfErr_(max(mag(perror_))),
    pL2err_(sqrt(sum(magSqr(perror_)))),


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
    // Set analytical fields
    if(calcAnalytical_)
    {
        setAnalyticalFields();
    }

    // Set global properties output file
    setPropertiesOutput();

    // get molucular viscosity
    getNu();

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TaylorGreenVortex2D::~TaylorGreenVortex2D()
{
    Info << "TaylorGreenVortex2D Destructor" << endl;
}

// ************************************************************************* //
