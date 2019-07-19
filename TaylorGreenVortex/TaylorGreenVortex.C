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

#include "TaylorGreenVortex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(TaylorGreenVortex, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::TaylorGreenVortex::setInitialFieldsAsAnalytical()
{
   if(dimension_ == "2D")
    {
        TaylorGreenVortex2D::setInitialFieldsAsAnalytical();
    }
    else
    {
        TaylorGreenVortex3D::setInitialFieldsAsAnalytical();
    }

    return;
}


void Foam::TaylorGreenVortex::setupProperties()
{
   if(dimension_ == "2D")
    {
        TaylorGreenVortex2D::setupProperties();
    }
    else
    {
        TaylorGreenVortex3D::setupProperties();
    }

    return;
}

void Foam::TaylorGreenVortex::calcProperties()
{
    if(dimension_ == "2D")
    {
        // calculate error fields and error norms
        TaylorGreenVortex2D::calcError();

        // calculate Global Ek and epsilon values
        TaylorGreenVortex2D::calcGlobalProperties();
    }
    else
    {
        // calculate Global Ek and epsilon values
        TaylorGreenVortex3D::calcGlobalProperties();
    }

    return;
}


void Foam::TaylorGreenVortex::writeProperties()
{
    // write to log files
    if(dimension_ == "2D")
    {
        TaylorGreenVortex2D::write();
    }
    else
    {
        TaylorGreenVortex3D::write();
    }

    return;
}


void Foam::TaylorGreenVortex::setPrecision(label nPrecision)
{
    if(dimension_ == "2D")
    {
        TaylorGreenVortex2D::setPrecision(nPrecision);
    }
    else
    {
        TaylorGreenVortex3D::setPrecision(nPrecision);
    }

    return;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TaylorGreenVortex::TaylorGreenVortex
(
    // const volVectorField& U,
    // const surfaceScalarField& phi,
    // const volScalarField& p,
    // const label& pRefCell

    volVectorField& U,
    surfaceScalarField& phi,
    volScalarField& p,
    label& pRefCell
)
:
    // constructor of parent classes
    TaylorGreenVortex2D (U, phi, p, pRefCell),
    TaylorGreenVortex3D (U, phi, p),

    // Set the pointer to runTime
    runTime_(U.time()),

    // Set the pointer to the mesh
    mesh_(U.mesh()),

    TGVProperties_
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

    // dimension_(TGVProperties_.lookupOrDefault<word>("dimension", "1D"))
    dimension_(TGVProperties_.lookup("dimension"))
{
    Info << "TaylorGreenVortex constructor" << endl;

    if (dimension_ == "2D")
    {
        TaylorGreenVortex2D::getNu();
        Info << "Taylor-Green Vortex dimension is " << dimension_ << endl;
    }
    else if (dimension_ == "3D")
    {
        TaylorGreenVortex3D::getNu();
        Info << "Taylor-Green Vortex dimension is " << dimension_ << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Taylor-Green Vortex dimension " << dimension_ << " is invalid!"
            << "\n Options are 2D or 3D." << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TaylorGreenVortex::~TaylorGreenVortex()
{
    Info << "TaylorGreenVortex Destructor" << endl;
}

// ************************************************************************* //
