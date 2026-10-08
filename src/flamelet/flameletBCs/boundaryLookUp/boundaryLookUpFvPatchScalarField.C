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

\*---------------------------------------------------------------------------*/

#include "boundaryLookUpFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryLookUpFvPatchScalarField::boundaryLookUpFvPatchScalarField(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF)
: fixedValueFvPatchScalarField(p, iF),
  lookupType_("standardBoundaryLookup")
{
}


Foam::boundaryLookUpFvPatchScalarField::boundaryLookUpFvPatchScalarField(
    const boundaryLookUpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper)
: fixedValueFvPatchScalarField(ptf, p, iF, mapper),
  lookupType_(ptf.lookupType_)
{
}


Foam::boundaryLookUpFvPatchScalarField::boundaryLookUpFvPatchScalarField(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict)
: fixedValueFvPatchScalarField(p, iF, dict),
  lookupType_(dict.lookupOrDefault<word>("lookupType", "standardBoundaryLookup"))
{
}


Foam::boundaryLookUpFvPatchScalarField::boundaryLookUpFvPatchScalarField(
    const boundaryLookUpFvPatchScalarField& ptf)
: fixedValueFvPatchScalarField(ptf),
  lookupType_(ptf.lookupType_)
{
}

Foam::boundaryLookUpFvPatchScalarField::boundaryLookUpFvPatchScalarField(
    const boundaryLookUpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF)
: fixedValueFvPatchScalarField(ptf, iF),
  lookupType_(ptf.lookupType_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryLookUpFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::boundaryLookUpFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    this->writeEntry("value", os);
    os.writeKeyword("lookupType") << lookupType_ << token::END_STATEMENT << nl;
}


const Foam::word Foam::boundaryLookUpFvPatchScalarField::getLookupType() const
{
    return lookupType_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField(
    fvPatchScalarField,
    boundaryLookUpFvPatchScalarField);
}

// ************************************************************************* //
