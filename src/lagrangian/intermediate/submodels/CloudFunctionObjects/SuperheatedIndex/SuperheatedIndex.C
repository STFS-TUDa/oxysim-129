/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "SuperheatedIndex.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SuperheatedIndex<CloudType>::SuperheatedIndex(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName)
: CloudFunctionObject<CloudType>(dict, owner, modelName, typeName)
{
}

// gp=gas props
template<class CloudType>
Foam::SuperheatedIndex<CloudType>::SuperheatedIndex(
    const SuperheatedIndex<CloudType>& rp)
: CloudFunctionObject<CloudType>(rp)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SuperheatedIndex<CloudType>::postEvolve(
    const typename parcelType::trackingData& td)
{
    const auto& c = this->owner();


    if (!c.template foundObject<IOField<scalar>>("Rp"))
    {
        IOField<scalar>* RpPtr =
            new IOField<scalar>(
                IOobject(
                    "Rp",
                    c.time().timeName(),
                    c,
                    IOobject::NO_READ));

        RpPtr->store();
    }


    auto& Rp = c.template lookupObjectRef<IOField<scalar>>("Rp");
    Rp.setSize(c.size());

    const auto& thermo = c.db().template lookupObject<SLGThermo>("SLGThermo");
    const auto& liquids = thermo.liquids();

    const auto& UInterp = td.UInterp();
    const auto& pInterp = td.pInterp();
    const auto& rhoInterp = td.rhoInterp();
    const auto& TInterp = td.TInterp();
    const auto& muInterp = td.muInterp();

    label parceli = 0;
    forAllConstIters(c, parcelIter)
    {
        const parcelType& p = parcelIter();

        const auto& coords = p.coordinates();
        const auto& tetIs = p.currentTetIndices();

        const vector Uc(UInterp.interpolate(coords, tetIs));

        const scalar pc =
            max(
                pInterp.interpolate(coords, tetIs),
                c.constProps().pMin());

        const scalar rhoc(rhoInterp.interpolate(coords, tetIs));
        const scalar muc(muInterp.interpolate(coords, tetIs));
        const scalar Tc(TInterp.interpolate(coords, tetIs));

        const scalarField X(liquids.X(p.YLiquid()));
        // saturation pressure for liquid species lid [Pa]
        const scalar pSat = liquids.pv(pc, p.T(), X);

        Rp[parceli++] = 50000.0 / pSat;
    }


    if (c.size() && c.time().writeTime())
    {
        Rp.write();
    }
}


// ************************************************************************* //
