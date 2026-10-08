/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::radiation::wsggAbsorptionEmission

Description
    Base class for weighted-sum-of-gray-gases-model (WSGGM) for gas phase absorption/emission.
    This class provides an interface for getWeights function for any wsgg model (Bordbar, Dorigon etc.)


SourceFiles
    wsggAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef wsggAbsorptionEmission_H
#define wsggAbsorptionEmission_H

#include "absorptionEmissionModel.H"
#include "fluidThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wsggAbsorptionEmission
: public absorptionEmissionModel
{
private:
    // Private data

    //- Absorption model dictionary
    dictionary coeffsDict_;

    //- Pointer to the thermo model
    const basicThermo* thermoPtr_;

    // double Mr_;

    label nBands_;

public:
    //- Runtime type information
    TypeName("wsggAbsorptionEmission");


    // Constructors

    //- Construct from components
    wsggAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggAbsorptionEmission();

    //- Member function
    virtual Foam::scalarList getWeights(
        Foam::scalar T,
        Foam::scalar Xr) const = 0;

    //- Initialize thermo class
    void setThermoPtr(
        const basicThermo& thermo);
    // Public getter for thermoPtr
    const basicThermo* thermoPtr() const
    {
        return thermoPtr_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
