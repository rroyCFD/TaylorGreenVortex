/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::TaylorGreenVortex2D

Description
    Coriolis force object

SourceFiles
    TaylorGreenVortex2D.C

\*---------------------------------------------------------------------------*/

#ifndef TaylorGreenVortex2D_H
#define TaylorGreenVortex2D_H

#include "fvCFD.H"
#include "IOdictionary.H"
#include "turbulenceModel.H"
#include "OFstream.H"
#include "vector.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class TaylorGreenVortex2D Declaration
\*---------------------------------------------------------------------------*/

class TaylorGreenVortex2D
{
    // Private data

        //- Runtime pointer
        const Time& runTime_;

        //- Mesh pointer
        const fvMesh& mesh_;

        // Pointer to solver fields
        //- Extrapolated velocity field pointer
        const volVectorField& U_;

        // Extrapolated surface volume flux pointer
        const surfaceScalarField& phi_;

        // pressure field p pointer
        const volScalarField& p_;

        const label& pRefCell_;

        const bool pRefOn_;

        // analytical velocity, flux and pressure field
        volVectorField Ua_;
        surfaceScalarField phia_;
        volScalarField pa_;

        // error fields and error norms
        volVectorField Uerror_;
        volScalarField perror_;

        dimensionedScalar ULinfErr_;
        dimensionedScalar UL2err_;

        dimensionedScalar pLinfErr_;
        dimensionedScalar pL2err_;

        // Case characteristic properties
        IOdictionary TGV2DProperties_;
        dimensionedScalar Uinit_;    // initial Velocity
        dimensionedScalar L_;        //  characteristic length scale
        dimensionedScalar tc_;       // characteristic time scale
        dimensionedScalar epsilon_;  // dissipation rate
        dimensionedScalar Ek_;       // kinetic energy

        autoPtr<dimensionedScalar> nu_;

        //- Output file pointer
        autoPtr<OFstream> globalPropertiesFile_;

    // Private member functions


public:

    //- Declare name of the class and its debug switch
    ClassName("TaylorGreenVortex2D");

    // Constructors
    TaylorGreenVortex2D
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& p,
        const label& pRefCell
    );

    // Destructor
    virtual ~TaylorGreenVortex2D();

    // Public Member functions
        volVectorField getUAnalytical()
        {
            return Ua_;
        }

        surfaceScalarField getPhiAnalytical()
        {
            return phia_;
        }

        volScalarField getPAnalytical()
        {
            return pa_;
        }

        dimensionedScalar epsilon()
        {
            return epsilon_;
        }

        dimensionedScalar Ek()
        {
            return Ek_;
        }


        void calcGlobalProperties();

        void setInitialFields(
            volVectorField UIn, surfaceScalarField phiIn, volScalarField pIn
        );

        void setInitialFieldsAsAnalytical();

        void setPropertiesOutput();

        void calcError();

        void write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
