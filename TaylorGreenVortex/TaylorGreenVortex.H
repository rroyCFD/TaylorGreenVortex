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
    Foam::TaylorGreenVortex

Description
    Coriolis force object

SourceFiles
    TaylorGreenVortex.C

\*---------------------------------------------------------------------------*/

#ifndef TaylorGreenVortex_H
#define TaylorGreenVortex_H

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
                           Class TaylorGreenVortex Declaration
\*---------------------------------------------------------------------------*/

class TaylorGreenVortex
{
    // Private data

        //- Runtime pointer
        const Time& runTime_;

        //- Mesh pointer
        const fvMesh& mesh_;

    // Private member functions


public:

    //- Declare name of the class and its debug switch
    ClassName("TaylorGreenVortex");

    // Constructors
    TaylorGreenVortex
    (
        const volVectorField& U
    );

    // Destructor
    virtual ~TaylorGreenVortex();

    // Public Member functions
        virtual void setInitialFieldsAsAnalytical(
            volVectorField& UIn,
            surfaceScalarField& phiIn,
            volScalarField& pIn
        ) = 0;

        virtual void setPropertiesOutput() = 0;

        virtual void calcGlobalProperties() = 0;

        virtual void writeGlobalProperties() = 0;

        virtual void setPrecision(label nPrecision) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
