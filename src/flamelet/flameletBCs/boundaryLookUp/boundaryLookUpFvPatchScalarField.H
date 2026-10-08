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
    Foam::boundaryLookUpFvPatchScalarField

Description
    A fixed boundary condition for enthalpy

SourceFiles
    boundaryLookUpFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryLookUpFvPatchScalarField_H
#define boundaryLookUpFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class boundaryLookUpFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
/** @defgroup boundaryConditionsFlamelet
 *  @ingroup flamelet
 *  @{
 */
class boundaryLookUpFvPatchScalarField
: public fixedValueFvPatchScalarField
{

    const word lookupType_;

public:
    //- Runtime type information
    TypeName("boundaryLookup");

    // Constructors

    //- Construct from patch and internal field
    boundaryLookUpFvPatchScalarField(
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&);

    //- Construct from patch, internal field and dictionary
    boundaryLookUpFvPatchScalarField(
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&);

    //- Construct by mapping given boundaryLookUpFvPatchScalarField
    // onto a new patch
    boundaryLookUpFvPatchScalarField(
        const boundaryLookUpFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&);

    //- Construct as copy
    boundaryLookUpFvPatchScalarField(
        const boundaryLookUpFvPatchScalarField&);

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const
    {
        return tmp<fvPatchScalarField>(
            new boundaryLookUpFvPatchScalarField(*this));
    }

    //- Construct as copy setting internal field reference
    boundaryLookUpFvPatchScalarField(
        const boundaryLookUpFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone(
        const DimensionedField<scalar, volMesh>& iF) const
    {
        return tmp<fvPatchScalarField>(
            new boundaryLookUpFvPatchScalarField(*this, iF));
    }


    // Member functions

    // Evaluation functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;

    int getPatchID();

    void init(const fvPatch& p);

    const word getLookupType() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**@} */
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
