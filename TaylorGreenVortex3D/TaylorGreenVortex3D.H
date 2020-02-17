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
    Foam::TaylorGreenVortex3D

Description
    Taylor-Green Vortex 3D object

SourceFiles
    TaylorGreenVortex3D.C

\*---------------------------------------------------------------------------*/

#ifndef TaylorGreenVortex3D_H
#define TaylorGreenVortex3D_H

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
                           Class TaylorGreenVortex3D Declaration
\*---------------------------------------------------------------------------*/

class TaylorGreenVortex3D
{
    // Private data

        //- Runtime pointer
        const Time& runTime_;

        //- Mesh pointer
        const fvMesh& mesh_;

        // Pointer to solver fields
        //- Extrapolated velocity field pointer
        // const volVectorField& U_;
        volVectorField& U_;

        // Extrapolated surface volume flux pointer
        // const surfaceScalarField& phi_;
        surfaceScalarField& phi_;

        // pressure field p pointer
        // const volScalarField& p_;
        volScalarField& p_;

        // Case characteristic properties
        IOdictionary TGV3DProperties_;
        dimensionedScalar Uinit_;    // initial Velocity
        dimensionedScalar L_;        //  characteristic length scale
        dimensionedScalar tc_;       // characteristic time scale
        dimensionedScalar epsilon_;  // dissipation rate
        dimensionedScalar viscEpsilon_;  // viscous dissipation rate
        dimensionedScalar lapTurbEpsilon_ ; // laplacian cmpt of turb dissip
        dimensionedScalar divGradTurbEpsilon_;// divGrad cmpt of turb dissip
        dimensionedScalar turbEpsilon_;  // total turbulent dissipation rate

        dimensionedScalar Ek_;       // kinetic energy

        dimensionedScalar nu_;
        autoPtr<dimensionedScalar> nuPtr_;

        //- Output file pointer
        autoPtr<OFstream> globalPropertiesFile_;

    // Private member functions


public:

    //- Declare name of the class and its debug switch
    ClassName("TaylorGreenVortex3D");

    // Constructors
    TaylorGreenVortex3D
    (
        // const volVectorField& U,
        // const surfaceScalarField& phi,
        // const volScalarField& p
        volVectorField& U,
        surfaceScalarField& phi,
        volScalarField& p
    );

    // Destructor
    virtual ~TaylorGreenVortex3D();

    // Public Member functions
        inline dimensionedScalar epsilon()
        {
            return epsilon_;
        }

        inline dimensionedScalar Ek()
        {
            return Ek_;
        }

        // call-up the initial member functions in case of parent class
        void setupProperties();

        void setInitialFieldsAsAnalytical(
            volVectorField& UIn, surfaceScalarField& phiIn, volScalarField& pIn
        );
        void setInitialFieldsAsAnalytical();

        void getNu();

        void setPropertiesOutput();

        void calcGlobalProperties();

        void write();

        void setPrecision(label nPrecision)
        {
            globalPropertiesFile_().precision(nPrecision);
            Info << "Output precision: "
                 << globalPropertiesFile_().precision() << endl;

            return;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
