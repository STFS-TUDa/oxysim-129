/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::radiation::wsggBinaryAbsorptionEmission

Group
    grpRadiationAbsorptionEmissionSubModels

Description
    Radiation coefficient based on two absorption models to be used in wsgg.
    First model is a gas phase wsgg model. Second model is a cloud model, using the same wsgg coeffs from the first model.

SourceFiles
    wsggBinaryAbsorptionEmission.C

\*---------------------------------------------------------------------------*/

#ifndef radiation_wsggBinaryAbsorptionEmission_H
#define radiation_wsggBinaryAbsorptionEmission_H

#include "wsggAbsorptionEmission.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class wsggBinaryAbsorptionEmission Declaration
\*---------------------------------------------------------------------------*/

class wsggBinaryAbsorptionEmission
: public absorptionEmissionModel
{
    // Private data

    //- Coefficients dictionary
    dictionary coeffsDict_;

    //- Pointer to the thermo model
    const basicThermo* thermoPtr_;

public:
    //- First absorption model from the gas phase. Needs to be wsgg type
    autoPtr<wsggAbsorptionEmission> model1_;

    //- Second absorption model
    autoPtr<absorptionEmissionModel> model2_;

    //- nBands of overall model
    label nBands_;

    // public:

    //- Runtime type information
    TypeName("wsggBinaryAbsorptionEmission");


    // Constructors

    //- Construct from components
    wsggBinaryAbsorptionEmission(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~wsggBinaryAbsorptionEmission();


    // Member Functions

    // Access

    // Absorption coefficient

    //- Absorption coefficient for continuous phase
    virtual tmp<volScalarField> aCont(const label bandI = 0) const;

    //- Absorption coefficient for dispersed phase
    virtual tmp<volScalarField> aDisp(const label bandI = 0) const;


    // Emission coefficient

    //- Emission coefficient for continuous phase
    virtual tmp<volScalarField> eCont(const label bandI = 0) const;

    //- Emission coefficient for dispersed phase
    virtual tmp<volScalarField> eDisp(const label bandI = 0) const;


    // Emission contribution

    //- Emission contribution for continuous phase
    virtual tmp<volScalarField> ECont(const label bandI = 0) const;

    //- Emission contribution for continuous phase
    virtual tmp<volScalarField> EDisp(const label bandI = 0) const;

    inline label nBands() const
    {
        return nBands_;
    }

    // Access
    void correct(
        Foam::volScalarField& a,
        PtrList<Foam::volScalarField>& aLambda) const;

    //- Set pointer to thermo class
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
